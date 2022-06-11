import re
import os
import glob
from collections import namedtuple
from subprocess import Popen, PIPE
from contextlib import contextmanager


import numpy as np
from scipy.interpolate import make_interp_spline


def get_file_size(path):
    """Return file size"""
    if os.path.isfile(path):
        size = os.stat(path).st_size
        for unit in ['Bytes','KB','MB','GB','TB']:
            if size < 1024.0:
                return "%3.2f %s" % (size,unit)
            size /= 1024.0
    else:
        return ''

@contextmanager
def set_dir(path):
    "work in some directory and come back"
    current = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(current)

def interpolate_data(x,y,n=10,k=3):
    """
    - Returns interpolated xnew,ynew. If two points are same, it will add 0.1*min(dx>0) to compensate it.
    - **Parameters**
        - x: 1D array of size p,
        - y: ndarray of size p*q*r,....
        - n: Number of points to add between two given points.
        - k: Polynomial order to interpolate.

    - Only axis 0 will be interpolated. If you want general interploation, use `from scipy.interpolate import make_interp_spline, BSpline`

    - **General Usage**: K(p),E(p,q) input from bandstructure.
        - `Knew,Enew= interpolate_data(K,E,n=10,k=3)`. cubic interploation
    """
    #Add very small values at simliar points to make interpolation work.
    ind=[i for i in range(0,len(x)) if x[i-1]==x[i]] #Duplicate indices
    xa=np.unique(x)
    dx=0.1*np.min(xa[1:]-xa[:-1])
    if(ind):
        for pt in ind:
            x[pt:]=x[pt:]-x[pt]+x[pt-1]+dx
    # Now Apply interpolation
    xnew=[np.linspace(x[i],x[i+1],n) for i in range(len(x)-1)]
    xnew=np.reshape(xnew,(-1))
    spl = make_interp_spline(x, y, k=k) #BSpline object
    ynew = spl(xnew)
    return xnew,ynew

def ps2py(ps_command='Get-ChildItem', exec_type='-Command', path_to_ps='powershell.exe'):
    """
    - Captures powershell output in python.
    - **Parameters**
        - ps_command: enclose ps_command in ' ' or " ".
        - exec_type : type of execution, default '-Command', could be '-File'.
        - path_to_ps: path to powerhell.exe if not added to PATH variables.
    """
    try: # Works on Linux and Windows if PS version > 5.
        cmd = ['pwsh', '-ExecutionPolicy', 'Bypass', exec_type, ps_command]
        proc = Popen(cmd, stdout=PIPE, stderr=PIPE)
    except FileNotFoundError:
        try: # Works only on Windows.
            cmd = ['powershell', '-ExecutionPolicy', 'Bypass', exec_type, ps_command]
            proc = Popen(cmd, stdout=PIPE, stderr=PIPE)
        except FileNotFoundError:
            # Works in case nothing above works and you know where is executable.
            cmd = [path_to_ps, '-ExecutionPolicy', 'Bypass', exec_type, ps_command]
            proc = Popen(cmd, stdout=PIPE, stderr=PIPE)

    out=[]; #save to out.
    while True:
        line = proc.stdout.readline()
        if line!=b'':
            line=line.strip()
            u_line=line.decode("utf-8")
            out.append(u_line)
        else:
            break
    out=[item for item in out if item!=''] #filter out empty lines
    return out

def ps2std(ps_command='Get-ChildItem', exec_type='-Command', path_to_ps='powershell.exe'):
    """
    - Prints powershell output in python std.
    - **Parameters**
        - ps_command: enclose ps_command in ' ' or " ".
        - exec_type: type of execution, default '-Command', could be '-File'.
        - path_to_ps: path to powerhell.exe if not added to PATH variables.
    """
    out = ps2py(path_to_ps=path_to_ps,exec_type=exec_type,ps_command=ps_command)
    for item in out:
        print(item)
    return None


def get_child_items(path = os.getcwd(),depth=None,recursive=True,include=None,exclude=None,filesOnly=False,dirsOnly= False):
    """
    - Returns selected directories/files recursively from a parent directory.
    - **Parameters**
        - path    : path to a parent directory, default is `"."`
        - depth   : int, subdirectories depth to get recursively, default is None to list all down.
        - recursive : If False, only list current directory items, if True,list all items recursively down the file system.
        - include: Default is None and includes everything. String of patterns separated by | to keep, could be a regular expression.
        - exclude: Default is None and removes nothing. String of patterns separated by | to drop,could be a regular expression.
        - filesOnly : Boolean, if True, returns only files.
        - dirsOnly  : Boolean, if True, returns only directories.
    - **Returns**
        - GLOB : Tuple (children,parent), children is list of selected directories/files and parent is given path. Access by index of by `get_child_items().{children,path}`.
    """
    path = os.path.abspath(path) # important
    pattern = path + '**/**' # Default pattern
    if depth != None and type(depth) == int:
        pattern = path + '/'.join(['*' for i in range(depth+1)])
        if glob.glob(pattern) == []: #If given depth is more, fall back.
            pattern = path + '**/**' # Fallback to default pattern if more depth to cover all.
    glob_files = glob.iglob(pattern, recursive=recursive)
    if dirsOnly == True:
        glob_files = filter(lambda f: os.path.isdir(f),glob_files)
    if filesOnly == True:
        glob_files = filter(lambda f: os.path.isfile(f),glob_files)
    list_dirs=[]
    for g_f in glob_files:
        list_dirs.append(os.path.relpath(g_f,path))
    # Include check
    if include:
        list_dirs = [l for l in list_dirs if re.search(include,l)]
    # Exclude check
    if exclude:
        list_dirs = [l for l in list_dirs if not re.search(exclude,l)]
    # Keep only unique
    req_dirs = list(np.unique(list_dirs))
    out_files = namedtuple('GLOB',['children','parent'])
    return out_files(req_dirs,os.path.abspath(path))

def prevent_overwrite(path):
    """Prevents overwiting as file/directory by adding numbers in given file/directory path."""
    if os.path.exists(path):
        name, ext = os.path.splitext(path)
        # Check existing files
        i = 0
        _path = name + '-{}' + ext
        while os.path.isfile(_path.format(i)):
            i +=1
        out_path = _path.format(i)
        print(f"Found existing path: {path!r}\nConverting to: {out_path!r}")
        return out_path
    return path

class color:
     def bg(text,r,g,b):
          """Provide r,g,b component in range 0-255"""
          return f"\033[48;2;{r};{g};{b}m{text}\033[00m"
     def fg(text,r,g,b):
          """Provide r,g,b component in range 0-255"""
          return f"\033[38;2;{r};{g};{b}m{text}\033[00m"
     # Usual Colos
     r  = lambda text: f"\033[0;91m {text}\033[00m"
     rb = lambda text: f"\033[1;91m {text}\033[00m"
     g  = lambda text: f"\033[0;92m {text}\033[00m"
     gb = lambda text: f"\033[1;92m {text}\033[00m"
     b  = lambda text: f"\033[0;34m {text}\033[00m"
     bb = lambda text: f"\033[1;34m {text}\033[00m"
     y  = lambda text: f"\033[0;93m {text}\033[00m"
     yb = lambda text: f"\033[1;93m {text}\033[00m"
     m  = lambda text: f"\033[0;95m {text}\033[00m"
     mb = lambda text: f"\033[1;95m {text}\033[00m"
     c  = lambda text: f"\033[0;96m {text}\033[00m"
     cb = lambda text: f"\033[1;96m {text}\033[00m"

def nav_links(current_index=0,
            doc_url = r"https://massgh.github.io/pivotpy/",
            items   = ["Index",
                       "Example",
                       "StaticPlots",
                       "InteractivePlots",
                       "SpinProjectedSurfaces",
                       "StructureIO",
                       "Widgets",
                       "MainAPI",
                       ],
            horizontal = False,
            out_string = False):
    from IPython.display import Markdown
    links   = [doc_url+item if not 'Index' in item else doc_url for item in items]
    style = """<style>a{text-decoration: none !important;color:lightkblue;font-weight:bold;}
                a:focus,a:active,a:hover{color:hotpink !important;}</style>\n"""
    md_str = style
    for i,(link,item) in enumerate(zip(links,items)):
        if current_index == i: item = "{}●".format(item)
        if not horizontal:
            md_str += "> [&nbsp;`▶` {}&nbsp;]({})  \n".format(item,link)
        else:
            md_str += "> [&nbsp;`▶` {}&nbsp;]({})\n".format(item,link)
    if out_string:
        return md_str
    return Markdown(md_str)

def transform_color(arr,s=1,c=1,b=0,mixing_matrix=None):
    """
    - Color transformation such as brightness, contrast, saturation and mixing of an input color array. `c = -1` would invert color,keeping everything else same.
    - **Parameters**
        - arr: input array, a single RGB/RGBA color or an array with inner most dimension equal to 3 or 4. e.g. [[[0,1,0,1],[0,0,1,1]]].
        - c  : contrast, default is 1. Can be a float in [-1,1].
        - s  : saturation, default is 1. Can be a float in [-1,1]. If s = 0, you get a gray scale image.
        - b  : brightness, default is 0. Can be a float in [-1,1] or list of three brightnesses for RGB components.
        - mixing_matrix: A 3x3 matrix to mix RGB values, such as `pp.color_matrix`.

    [Recoloring](https://docs.microsoft.com/en-us/windows/win32/gdiplus/-gdiplus-recoloring-use?redirectedfrom=MSDN)
    [Rainmeter](https://docs.rainmeter.net/tips/colormatrix-guide/)
    """
    arr = np.array(arr) # Must
    t = (1-c)/2 # For fixing gray scale when contrast is 0.
    whiteness = np.array(b)+t # need to clip to 1 and 0 after adding to color.
    sr = (1-s)*0.2125 #red saturation from red luminosity
    sg = (1-s)*0.7154 #green saturation from green luminosity
    sb = (1-s)*0.0721 #blue saturation from blue luminosity
    # trans_matrix is multiplied from left, or multiply its transpose from right.
    # trans_matrix*color is not normalized but value --> value - int(value) to keep in [0,1].
    trans_matrix = np.array([
        [c*(sr+s), c*sg,      c*sb],
        [c*sr,   c*(sg+s),    c*sb],
        [c*sr,     c*sg,  c*(sb+s)]])
    if np.ndim(arr) == 1:
        new_color = np.dot(trans_matrix,arr)
    else:
        new_color = np.dot(arr[...,:3],trans_matrix.T)
    if mixing_matrix is not None and np.size(mixing_matrix)==9:
        new_color = np.dot(new_color,np.transpose(mixing_matrix))
    new_color[new_color > 1] = new_color[new_color > 1] - new_color[new_color > 1].astype(int)
    new_color = np.clip(new_color + whiteness,a_max=1,a_min=0)
    if np.shape(arr)[-1]==4:
        axis = len(np.shape(arr))-1 #Add back Alpha value if present
        new_color = np.concatenate([new_color,arr[...,3:]],axis=axis)
    return new_color