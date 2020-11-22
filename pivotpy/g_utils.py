# AUTOGENERATED! DO NOT EDIT! File to edit: Utilities.ipynb (unless otherwise specified).

__all__ = ['get_file_size', 'interpolate_data', 'ps_to_py', 'ps_to_std', 'select_dirs', 'select_files',
           'get_child_items', 'invert_color', 'printr', 'printg', 'printb', 'printy', 'printm', 'printc',
           'EncodeFromNumpy', 'DecodeToNumpy', 'link_to_class', 'nav_links']

# Cell
def get_file_size(path):
    import os
    if os.path.isfile(path):
        size = os.stat(path).st_size
        for unit in ['Bytes','KB','MB','GB','TB']:
            if size < 1024.0:
                return "%3.2f %s" % (size,unit)
            size /= 1024.0
    else:
        return ''

# Cell
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
    import numpy as np
    #Add very small values at simliar points to make interpolation work.
    ind=[i for i in range(0,len(x)) if x[i-1]==x[i]] #Duplicate indices
    xa=np.unique(x)
    dx=0.1*np.min(xa[1:]-xa[:-1])
    if(ind):
        for pt in ind:
            x[pt:]=x[pt:]-x[pt]+x[pt-1]+dx
    # Now Apply interpolation
    from scipy.interpolate import make_interp_spline, BSpline
    xnew=[np.linspace(x[i],x[i+1],n) for i in range(len(x)-1)]
    xnew=np.reshape(xnew,(-1))
    spl = make_interp_spline(x, y, k=k) #BSpline object
    ynew = spl(xnew)
    return xnew,ynew

# Cell
def ps_to_py(ps_command='Get-ChildItem', exec_type='-Command', path_to_ps='powershell.exe'):
    """
    - Captures powershell output in python.
    - **Parameters**
        - ps_command: enclose ps_command in ' ' or " ".
        - exec_type : type of execution, default '-Command', could be '-File'.
        - path_to_ps: path to powerhell.exe if not added to PATH variables.
    """
    from subprocess import Popen, PIPE
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

# Cell
def ps_to_std(ps_command='Get-ChildItem', exec_type='-Command', path_to_ps='powershell.exe'):
    """
    - Prints powershell output in python std.
    - **Parameters**
        - ps_command: enclose ps_command in ' ' or " ".
        - exec_type: type of execution, default '-Command', could be '-File'.
        - path_to_ps: path to powerhell.exe if not added to PATH variables.
    """
    out=ps_to_py(path_to_ps=path_to_ps,exec_type=exec_type,ps_command=ps_command)
    for item in out:
        print(item)
    return None

# Cell
import os
import glob
#Selection of required project directories.
def select_dirs(path = os.getcwd(),include=[],exclude=[]):
    """
    - Returns selected directories recursively from a parent directory.
    - **Parameters**
        - path    : path to a parent directory, default is `"."`
        - include : list of keywords to include directories, avoid wildcards.
        - exclude : list of keywords to exclude directories, avoid wildcards.
    - **Returns**
        - Tuple of two elements, list of selcted directories and given path.
    """
    print('Use command `get_child_items()` instead for more flexibility.')
    list_dirs=[]; req_dirs=[];
    for filename in glob.iglob(path + '**/**', recursive=True):
        if os.path.isdir(filename):
            list_dirs.append(filename)
    for item in list_dirs:
        for check in include:
            if(check in item):
                if(path != os.getcwd()):
                    req_dirs.append(item.replace("\\","/"))
                if(path == os.getcwd()):
                    req_dirs.append('.'+(item.split(os.getcwd())[-1]).replace("\\","/"))
    for item in req_dirs.copy():
        for ex in exclude:
            if ex in item:
                req_dirs.remove(item)
    return (req_dirs,path.replace("\\","/"))
#Selction of files in selected directories.
def select_files(path=os.getcwd(),include=[],exclude=[]):
    """
    - Returns selected files from a given directory.
    - **Parameters**
        - path    : path to a parent directory, default is `"."`
        - include : list of keywords to include files, avoid wildcards.
        - exclude : list of keywords to exclude files, avoid wildcards.
    - **Returns**
        - Tuple of two elements, list of selcted files and given path.
    """
    print('Use command `get_child_items()` instead for more flexibility.')
    req_files=[]
    all_files=os.listdir(path)
    for file in all_files:
        for check in include:
                    if(check in file):
                        req_files.append(file)
    for item in req_files.copy():
        for ex in exclude:
            if ex in item:
                req_files.remove(item)
    return (req_files,path.replace("\\","/"))

# Cell
def get_child_items(path = os.getcwd(),depth=None,recursive=True,include=[],exclude=[],filesOnly=False,dirsOnly= False):
    """
    - Returns selected directories/files recursively from a parent directory.
    - **Parameters**
        - path    : path to a parent directory, default is `"."`
        - depth   : int, subdirectories depth to get recursively, default is None to list all down.
        - recursive : If False, only list current directory items, if True,list all items recursively down the file system.
        - include : list or str of keywords to include directories/files, avoid wildcards.
        - exclude : list or str of keywords to exclude directories/files, avoid wildcards.
        - filesOnly : Boolean, if True, returns only files.
        - dirsOnly  : Boolean, if True, returns only directories.
    - **Returns**
        - GLOB : Tuple (children,parent), children is list of selected directories/files and parent is given path. Access by index of by `get_child_items().{children,path}`.
    """
    import os
    import glob
    import numpy as np
    from collections import namedtuple
    if include != None and type(include) == str:
        include = [include,]
    if exclude != None and type(exclude) == str:
        exclude = [exclude,]
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
    req_dirs=[]
    if include != []:
        for check in include:
            req_dirs.extend(list(filter(lambda f: check in f ,list_dirs)))
    elif include == []:
        req_dirs = list_dirs
    # Exclude check
    to_exclude = []
    if exclude != []:
        for ex in exclude:
            to_exclude.extend(list(filter(lambda f: ex in f ,req_dirs)))
        req_dirs = [r_d for r_d in req_dirs if r_d not in to_exclude]
    # Keep only unique
    req_dirs = list(np.unique(req_dirs))
    out_files = namedtuple('GLOB',['children','parent'])
    return out_files(req_dirs,os.path.abspath(path))

# Cell
def invert_color(color=(1,1,1)):
    """
    - Returns opposite of given complementary color.
    - Input: Tuple (r,g,b).
    """
    r = min(color)+max(color)
    return tuple(r-c for c in color)

# Cell
def printr(s): print("\033[91m {}\033[00m" .format(s))
def printg(s): print("\033[92m {}\033[00m" .format(s))
def printb(s): print("\033[34m {}\033[00m" .format(s))
def printy(s): print("\033[93m {}\033[00m" .format(s))
def printm(s): print("\033[95m {}\033[00m" .format(s))
def printc(s): print("\033[96m {}\033[00m" .format(s))

# Cell
import json
class EncodeFromNumpy(json.JSONEncoder):
    """
    - Serializes python/Numpy objects via customizing json encoder.
    - **Usage**
        - `json.dumps(python_dict, cls=EncodeFromNumpy)` to get json string.
        - `json.dump(*args, cls=EncodeFromNumpy)` to create a file.json.
    """
    def default(self, obj):
        import numpy
        if isinstance(obj, numpy.ndarray):
            return {
                "_kind_": "ndarray",
                "_value_": obj.tolist()
            }
        if isinstance(obj, numpy.integer):
            return int(obj)
        elif isinstance(obj, numpy.floating):
            return float(obj)
        elif isinstance(obj,range):
            value = list(obj)
            return {
                "_kind_" : "range",
                "_value_" : [value[0],value[-1]+1]
            }
        return super(EncodeFromNumpy, self).default(obj)



class DecodeToNumpy(json.JSONDecoder):
    """
    - Deserilizes JSON object to Python/Numpy's objects.
    - **Usage**
        - `json.loads(json_string,cls=DecodeToNumpy)` from string, use `json.load()` for file.
    """
    def __init__(self, *args, **kwargs):
        json.JSONDecoder.__init__(self, object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, obj):
        import numpy
        if '_kind_' not in obj:
            return obj
        kind = obj['_kind_']
        if kind == 'ndarray':
            return numpy.array(obj['_value_'])
        elif kind == 'range':
            value = obj['_value_']
            return range(value[0],value[-1])
        return obj

# Cell
def link_to_class(cls):
    """
    - Binds wrapper of a function to class as attribute that does exactly the same as function. Also function returned from wrapper can be used normally as well.
    - **Parameters**
        - cls : A class object to which function is attached.
    """
    from functools import wraps
    def decorator(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            return func(*args, **kwargs)
        setattr(cls, func.__name__, wrapper)
        return func
    return decorator

# Cell
def nav_links(current_index=0,
            doc_url = r"https://massgh.github.io/pivotpy/",
            items   = ["Index",
                       "XmlElementTree",
                       "StaticPlots",
                       "InteractivePlots",
                       "Utilities",
                       "StructureIO",
                       "Widgets"
                       ]):
    from IPython.display import Markdown,HTML
    links   = [doc_url+item if not 'Index' in item else doc_url for item in items]
    style = """<style>
                .mydiv {background:#eaf0f0; padding:2px;display:inline-block;border:1px solid #93b2b2;}
                a{text-decoration: none;}
                a:focus,a:active.a:hover{color:hotpink;}
                a:visited{opacity:0.5;color:green;}
                </style>\n"""
    md_str = style
    for i,(link,item) in enumerate(zip(links,items)):
        if current_index == i: item = "●{}".format(item)
        md_str += "<div class='mydiv'><b><a href='{}'>  {}  </a></b></div>".format(link,item)
    return HTML(f"<div>{md_str}</div>")