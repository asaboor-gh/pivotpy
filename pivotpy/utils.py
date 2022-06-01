import re
import os
import glob
from collections import namedtuple
from subprocess import Popen, PIPE
from contextlib import contextmanager

from importlib.machinery import SourceFileLoader

import numpy as np
from scipy.interpolate import make_interp_spline

from . import serializer

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

def load_ps_exported(path= './vasprun.xml',kseg_inds =[], shift_kpath = 0, path_to_ps='pwsh', skipk = None, max_filled = 10, max_empty = 10, keep_files = True):
    """
    - Returns a full dictionary of all objects from `vasprun.xml` file exported using powershell.
    - **Parameters**
        - path       : Path to `vasprun.xml` file. Default is `'./vasprun.xml'`.
        - skipk      : Default is None. Automatically detects kpoints to skip.
        - path_to_ps : Path to `powershell.exe`. Automatically picks on Windows and Linux if added to PATH.
        - kseg_inds : List of indices of kpoints where path is broken.
        - shift_kpath: Default 0. Can be used to merge multiple calculations side by side.
        - keep_files : Could be use to clean exported text files. Default is True.
        - max_filled : Number of filled bands below and including VBM. Default is 10.
        - max_empty  : Number of empty bands above VBM. Default is 10.
    
    **Returns**: `pivotpy.serializer.VasprunData` object.
    """
    from . import parser # This should be inside function, not outside, otherwise circular import.
    that_loc, file_name = os.path.split(os.path.abspath(path)) # abspath is important to split.
    with set_dir(that_loc):
        # Goes there and work
        i = 0
        required_files = ['Bands.txt','tDOS.txt','pDOS.txt','Projection.txt','SysInfo.py']
        for _file in required_files:
            if os.path.isfile(_file):
               i = i + 1
        if i < 5:
            if skipk != None:
                ps2std(path_to_ps=path_to_ps,ps_command='Import-Module Vasp2Visual; Export-VR -InputFile {} -MaxFilled {} -MaxEmpty {} -SkipK {}'.format(path,max_filled,max_empty,skipk))
            else:
                ps2std(path_to_ps=path_to_ps,ps_command='Import-Module Vasp2Visual; Export-VR -InputFile {} -MaxFilled {} -MaxEmpty {}'.format(path,max_filled,max_empty))

        # get info from same directory
        iscartesian = parser._is_poscar_cartesian()
        kpoints_info = parser.get_kpoints_info()

        # Enable loading SysInfo.py file as source.
        _vars = SourceFileLoader("SysInfo", "./SysInfo.py").load_module()

        SYSTEM            = _vars.SYSTEM
        NKPTS             = _vars.NKPTS
        NBANDS            = _vars.NBANDS
        NFILLED           = _vars.NFILLED
        TypeION           = _vars.TypeION
        NION              = _vars.NION
        NELECT            = _vars.NELECT
        nField_Projection = _vars.nField_Projection
        E_Fermi           = _vars.E_Fermi
        ISPIN             = _vars.ISPIN
        ElemIndex         = _vars.ElemIndex
        ElemName          = _vars.ElemName
        poscar            = {'SYSTEM': SYSTEM,
                            'volume':_vars.volume,
                            'basis' : np.array(_vars.basis),
                            'rec_basis': np.array(_vars.rec_basis),
                            'cartesian':iscartesian,
                            'positions': np.array(_vars.positions)
                            }
        fields            = _vars.fields
        incar             = _vars.INCAR

        # Elements Labels
        elem_labels = []
        for i, name in enumerate(ElemName):
            for ind in range(ElemIndex[i],ElemIndex[i+1]):
                elem_labels.append(f"{name} {str(ind - ElemIndex[i] + 1)}")
        poscar.update({'labels': elem_labels})
        # Unique Elements Ranges
        unique_d = {}
        for i,e in enumerate(ElemName):
            unique_d.update({e:range(ElemIndex[i],ElemIndex[i+1])})
        poscar.update({'unique': unique_d})

        # Load Data
        bands= np.loadtxt('Bands.txt').reshape((-1,NBANDS+4)) #Must be read in 2D even if one row only.
        start = int(open('Bands.txt').readline().split()[4][1:])
        pro_bands= np.loadtxt('Projection.txt').reshape((-1,NBANDS*nField_Projection))
        pro_dos = np.loadtxt('pDOS.txt')
        dos= np.loadtxt('tDOS.txt')

        # Keep or delete only if python generates files (i < 5 case.)
        if(keep_files==False and i==5):
            for file in required_files:
                os.remove(file)
        # Returns back

    # Work now!
    sys_info = {'SYSTEM': SYSTEM,'NION': NION,'NELECT':NELECT,'TypeION': TypeION,'ElemName': ElemName,
                'E_Fermi': E_Fermi,'fields':fields, 'incar': incar,'ElemIndex': ElemIndex,'ISPIN': ISPIN,
                'kpts_info': kpoints_info}
    dim_info = {'kpoints': '(NKPTS,3)','kpath': '(NKPTS,1)','bands': '⇅(NKPTS,NBANDS)','dos': '⇅(grid_size,3)',
                'pro_dos': '⇅(NION,grid_size,en+pro_fields)','pro_bands': '⇅(NION,NKPTS,NBANDS,pro_fields)'}

    bands_dic,tdos_dic,pdos_dic,pro_dic,kpath={},{},{},{},[]
    if(ISPIN==1):
        kpath   = bands[:,3]
        kpoints = bands[:,:3]
        evals   = bands[:,4:]
        bands_dic = {'E_Fermi': E_Fermi, 'ISPIN': ISPIN, 'NBANDS': NBANDS, 'evals': evals, 'indices': range(start,start+NBANDS)}
        tdos_dic  = {'E_Fermi': E_Fermi, 'ISPIN': ISPIN,'tdos': dos}
        pdos      = pro_dos.reshape(NION,-1,nField_Projection+1)
        pdos_dic  = {'labels': fields,'pros': pdos}
        pros      = pro_bands.reshape(NION,NKPTS,NBANDS,-1)
        pro_dic   = {'labels': fields,'pros': pros}
    if(ISPIN==2):
        # Bands
        kpath   = bands[:NKPTS,3]
        kpoints = bands[:NKPTS,:3]
        SpinUp  = bands[:NKPTS,4:]
        SpinDown= bands[NKPTS:,4:]
        evals   = {'SpinUp':SpinUp,'SpinDown': SpinDown}
        bands_dic = {'E_Fermi': E_Fermi, 'ISPIN': ISPIN, 'NBANDS': NBANDS, 'evals': evals,'indices': range(start,start+NBANDS)}
        # tDOS
        dlen    = int(np.shape(dos)[0]/2)
        SpinUp  = dos[:dlen,:]
        SpinDown= dos[dlen:,:]
        tdos    = {'SpinUp':SpinUp,'SpinDown': SpinDown}
        tdos_dic= {'E_Fermi': E_Fermi, 'ISPIN': ISPIN,'tdos': tdos}

        # pDOS
        plen    = int(np.shape(pro_dos)[0]/2)
        SpinUp  = pro_dos[:plen,:].reshape(NION,-1,nField_Projection+1)
        SpinDown= pro_dos[plen:,:].reshape(NION,-1,nField_Projection+1)
        pdos    = {'SpinUp':SpinUp,'SpinDown': SpinDown}
        pdos_dic= {'labels': fields,'pros': pdos}

        # projections
        pblen  = int(np.shape(pro_bands)[0]/2)
        SpinUp  = pro_bands[:pblen,:].reshape(NION,NKPTS,NBANDS,-1)
        SpinDown= pro_bands[pblen:,:].reshape(NION,NKPTS,NBANDS,-1)
        pros    = {'SpinUp': SpinUp,'SpinDown': SpinDown}
        pro_dic = {'labels': fields,'pros': pros}
    # If broken path, then join points.
    kpath = parser.join_ksegments(kpath,kseg_inds)
    kpath=[k+shift_kpath for k in kpath.copy()] # Shift kpath
    full_dic = {'sys_info': sys_info,'dim_info': dim_info,'kpoints': kpoints,'kpath':kpath,               'bands':bands_dic,'tdos':tdos_dic,'pro_bands': pro_dic ,'pro_dos': pdos_dic,
               'poscar':serializer.PoscarData(poscar)}
    return serializer.VasprunData(full_dic)

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