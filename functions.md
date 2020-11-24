## Functions Reference


```python
from nbdev import show_doc
import pivotpy as pp
all_ = pp.__all__
show_doc(pp.Dict2Data)
show_doc(pp.Dict2Data.to_dict)
show_doc(pp.Dict2Data.to_json)
show_doc(pp.Dict2Data.to_pickle)
_ = [show_doc(eval('pp.{}'.format(f))) for f in all_ if f not in ['Dict2Data','savefig','show']]
```


<h2 id="Dict2Data" class="doc_header"><code>class</code> <code>Dict2Data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L10" class="source_link" style="float:right">[source]</a></h2>

> <code>Dict2Data</code>(**`d`**) :: `dict`

- Returns a Data object with dictionary keys as attributes of Data accessible by dot notation.
- **Parmeters**
    - dict : Python dictionary (nested as well) containing any python data types.
- **Methods**
    - to_dict()  : Converts a Data object to dictionary if it could be made a dictionary, otherwise throws relevant error.
    - to_json()  : Converts to json str or save to file if `outfil` given. Accepts `indent` as parameter.
    - to_pickle(): Converts to bytes str or save to file if `outfile` given.
- **Example**
    > x = Dict2Data({'A':1,'B':{'C':2}})
    > x
    > Data(
    >     A = 1
    >     B = Data(
    >         C = 2
    >         )
    >     )
    > x.B.to_dict()
    > {'C': 2}



<h4 id="Dict2Data.to_dict" class="doc_header"><code>Dict2Data.to_dict</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L41" class="source_link" style="float:right">[source]</a></h4>

> <code>Dict2Data.to_dict</code>()

- Converts a [`Dict2Data`](/pivotpy/XmlElementTree.html#Dict2Data) object (root or nested level) to a dictionary.



<h4 id="Dict2Data.to_json" class="doc_header"><code>Dict2Data.to_json</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L52" class="source_link" style="float:right">[source]</a></h4>

> <code>Dict2Data.to_json</code>(**`outfile`**=*`None`*, **`indent`**=*`1`*)

- Dumps a [`Dict2Data`](/pivotpy/XmlElementTree.html#Dict2Data) object (root or nested level) to json.
- **Parameters**
    - outfile : Default is None and returns string. If given, writes to file.
    - indent  : Json indent. Default is 1.



<h4 id="Dict2Data.to_pickle" class="doc_header"><code>Dict2Data.to_pickle</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L62" class="source_link" style="float:right">[source]</a></h4>

> <code>Dict2Data.to_pickle</code>(**`outfile`**=*`None`*)

- Dumps a [`Dict2Data`](/pivotpy/XmlElementTree.html#Dict2Data) object (root or nested level) to pickle.
- **Parameters**
    - outfile : Default is None and returns string. If given, writes to file.



<h4 id="read_asxml" class="doc_header"><code>read_asxml</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L90" class="source_link" style="float:right">[source]</a></h4>

> <code>read_asxml</code>(**`path`**=*`None`*, **`suppress_warning`**=*`False`*)

- Reads a big vasprun.xml file into memory once and then apply commands.
If current folder contains `vasprun.xml` file, it automatically picks it.

- **Parameters**
    - path             : Path/To/vasprun.xml
    - suppress_warning : False by defualt. Warns about memory usage for large files > 100 MB.
- **Returns**
    - xml_data : Xml object to use in other functions



<h4 id="exclude_kpts" class="doc_header"><code>exclude_kpts</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L141" class="source_link" style="float:right">[source]</a></h4>

> <code>exclude_kpts</code>(**`xml_data`**=*`None`*)

- Returns number of kpoints to exclude used from IBZKPT.
- **Parameters**
    - xml_data : From [`read_asxml`](/pivotpy/XmlElementTree.html#read_asxml) function
- **Returns**
    - int      : Number of kpoints to exclude.



<h4 id="get_ispin" class="doc_header"><code>get_ispin</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L162" class="source_link" style="float:right">[source]</a></h4>

> <code>get_ispin</code>(**`xml_data`**=*`None`*)

- Returns value of ISPIN.
- **Parameters**
    - xml_data : From [`read_asxml`](/pivotpy/XmlElementTree.html#read_asxml) function
- **Returns**
    - int      : Value of ISPIN.



<h4 id="get_summary" class="doc_header"><code>get_summary</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L179" class="source_link" style="float:right">[source]</a></h4>

> <code>get_summary</code>(**`xml_data`**=*`None`*)

- Returns overview of system parameters.
- **Parameters**
    - xml_data : From [`read_asxml`](/pivotpy/XmlElementTree.html#read_asxml) function
- **Returns**
    - Data     : pivotpy.Dict2Data with attibutes accessible via dot notation.



<h4 id="get_kpts" class="doc_header"><code>get_kpts</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L217" class="source_link" style="float:right">[source]</a></h4>

> <code>get_kpts</code>(**`xml_data`**=*`None`*, **`skipk`**=*`0`*, **`joinPathAt`**=*`[]`*)

- Returns kpoints and calculated kpath.
- **Parameters**
    - xml_data   : From [`read_asxml`](/pivotpy/XmlElementTree.html#read_asxml) function
    - skipk      : Number of initil kpoints to skip
    - joinPathAt : List of indices of kpoints where path is broken
- **Returns**
    - Data     : pivotpy.Dict2Data with attibutes `kpath` and `kpoints`



<h4 id="get_tdos" class="doc_header"><code>get_tdos</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L251" class="source_link" style="float:right">[source]</a></h4>

> <code>get_tdos</code>(**`xml_data`**=*`None`*, **`spin_set`**=*`1`*, **`elim`**=*`[]`*)

- Returns total dos for a spin_set (default 1) and energy limit. If spin-polarized calculations, gives SpinUp and SpinDown keys as well.
- **Parameters**
    - xml_data : From [`read_asxml`](/pivotpy/XmlElementTree.html#read_asxml) function
    - spin_set : int, default is 1.and
    - elim     : List [min,max] of energy, default empty.
- **Returns**
    - Data     : pivotpy.Dict2Data with attibutes E_Fermi, ISPIN,tdos.



<h4 id="get_evals" class="doc_header"><code>get_evals</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L304" class="source_link" style="float:right">[source]</a></h4>

> <code>get_evals</code>(**`xml_data`**=*`None`*, **`skipk`**=*`None`*, **`elim`**=*`[]`*)

- Returns eigenvalues as numpy array. If spin-polarized calculations, gives SpinUp and SpinDown keys as well.
- **Parameters**
    - xml_data : From [`read_asxml`](/pivotpy/XmlElementTree.html#read_asxml) function
    - skipk    : Number of initil kpoints to skip.
    - elim     : List [min,max] of energy, default empty.
- **Returns**
    - Data     : pivotpy.Dict2Data with attibutes evals and related parameters.



<h4 id="get_bands_pro_set" class="doc_header"><code>get_bands_pro_set</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L357" class="source_link" style="float:right">[source]</a></h4>

> <code>get_bands_pro_set</code>(**`xml_data`**=*`None`*, **`spin_set`**=*`1`*, **`skipk`**=*`0`*, **`bands_range`**=*`None`*)

- Returns bands projection of a spin_set(default 1) as numpy array. If spin-polarized calculations, gives SpinUp and SpinDown keys as well.
- **Parameters**
    - xml_data    : From [`read_asxml`](/pivotpy/XmlElementTree.html#read_asxml) function
    - skipk       : Number of initil kpoints to skip (Default 0).
    - spin_set    : Spin set to get, default is 1.
    - bands_range : If elim used in [`get_evals`](/pivotpy/XmlElementTree.html#get_evals),that will return bands_range to use here..
- **Returns**
    - Data     : pivotpy.Dict2Data with attibutes of bands projections and related parameters.



<h4 id="get_dos_pro_set" class="doc_header"><code>get_dos_pro_set</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L410" class="source_link" style="float:right">[source]</a></h4>

> <code>get_dos_pro_set</code>(**`xml_data`**=*`None`*, **`spin_set`**=*`1`*, **`dos_range`**=*`None`*)

- Returns dos projection of a spin_set(default 1) as numpy array. If spin-polarized calculations, gives SpinUp and SpinDown keys as well.
- **Parameters**
    - xml_data    : From [`read_asxml`](/pivotpy/XmlElementTree.html#read_asxml) function
    - spin_set    : Spin set to get, default 1.
    - dos_range   : If elim used in [`get_tdos`](/pivotpy/XmlElementTree.html#get_tdos),that will return dos_range to use here..
- **Returns**
    - Data     : pivotpy.Dict2Data with attibutes of dos projections and related parameters.



<h4 id="get_structure" class="doc_header"><code>get_structure</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L453" class="source_link" style="float:right">[source]</a></h4>

> <code>get_structure</code>(**`xml_data`**=*`None`*)

- Returns structure's volume,basis,positions and rec-basis.
- **Parameters**
    - xml_data : From [`read_asxml`](/pivotpy/XmlElementTree.html#read_asxml) function.
- **Returns**
    - Data     : pivotpy.Dict2Data with attibutes volume,basis,positions and rec_basis.



<h4 id="export_vasprun" class="doc_header"><code>export_vasprun</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L481" class="source_link" style="float:right">[source]</a></h4>

> <code>export_vasprun</code>(**`path`**=*`None`*, **`skipk`**=*`None`*, **`elim`**=*`[]`*, **`joinPathAt`**=*`[]`*, **`shift_kpath`**=*`0`*)

- Returns a full dictionary of all objects from `vasprun.xml` file. It first try to load the data exported by powershell's `Export-VR(Vasprun)`, which is very fast for large files. It is recommended to export large files in powershell first.
- **Parameters**
    - path       : Path to `vasprun.xml` file. Default is `'./vasprun.xml'`.
    - skipk      : Default is None. Automatically detects kpoints to skip.
    - elim       : List [min,max] of energy interval. Default is [], covers all bands.
    - joinPathAt : List of indices of kpoints where path is broken.
    - shift_kpath: Default 0. Can be used to merge multiple calculations on single axes side by side.
- **Returns**
    - Data : Data accessible via dot notation containing nested Data objects:
        - sys_info  : System Information
        - dim_info  : Contains information about dimensions of returned objects.
        - kpoints   : numpy array of kpoints with excluded IBZKPT points
        - kpath     : 1D numpy array directly accessible for plot.
        - bands     : Data containing bands.
        - tdos      : Data containing total dos.
        - pro_bands : Data containing bands projections.
        - pro_dos   : Data containing dos projections.
        - poscar    : Data containing basis,positions, rec_basis and volume.



<h4 id="load_export" class="doc_header"><code>load_export</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L566" class="source_link" style="float:right">[source]</a></h4>

> <code>load_export</code>(**`path`**=*`'./vasprun.xml'`*, **`joinPathAt`**=*`[]`*, **`shift_kpath`**=*`0`*, **`path_to_ps`**=*`'pwsh'`*, **`skipk`**=*`None`*, **`max_filled`**=*`10`*, **`max_empty`**=*`10`*, **`keep_files`**=*`True`*)

- Returns a full dictionary of all objects from `vasprun.xml` file exported using powershell.
- **Parameters**
    - path       : Path to `vasprun.xml` file. Default is `'./vasprun.xml'`.
    - skipk      : Default is None. Automatically detects kpoints to skip.
    - path_to_ps : Path to `powershell.exe`. Automatically picks on Windows and Linux if added to PATH.
    - joinPathAt : List of indices of kpoints where path is broken.
    - shift_kpath: Default 0. Can be used to merge multiple calculations side by side.
    - keep_files : Could be use to clean exported text files. Default is True.
    - max_filled : Number of filled bands below and including VBM. Default is 10.
    - max_empty  : Number of empty bands above VBM. Default is 10.
- **Returns**
    - Data : Data accessible via dot notation containing nested Data objects:
        - sys_info  : System Information
        - dim_info  : Contains information about dimensions of returned objects.
        - kpoints   : numpy array of kpoints with excluded IBZKPT points
        - kpath     : 1D numpy array directly accessible for plot.
        - bands     : Data containing bands.
        - tdos      : Data containing total dos.
        - pro_bands : Data containing bands projections.
        - pro_dos   : Data containing dos projections.
        - poscar    : Data containing basis,positions, rec_basis and volume.



<h4 id="dump_dict" class="doc_header"><code>dump_dict</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L720" class="source_link" style="float:right">[source]</a></h4>

> <code>dump_dict</code>(**`dict_data`**=*`None`*, **`dump_to`**=*`'pickle'`*, **`outfile`**=*`None`*, **`indent`**=*`1`*)

- Dump an [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun) or [`load_export`](/pivotpy/XmlElementTree.html#load_export)'s `Data` object or any dictionary to json or pickle string/file. It convert [`Dict2Data`](/pivotpy/XmlElementTree.html#Dict2Data) to dictionary before serializing to json/pickle, so json/pickle.loads() of converted Data would be a simple dictionary, pass that to [`Dict2Data`](/pivotpy/XmlElementTree.html#Dict2Data) to again make accessible via dot notation.
- **Parameters**
    - dict_data : Any dictionary/Dict2Data object containg numpy arrays, including [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun) or [`load_export`](/pivotpy/XmlElementTree.html#load_export) output.
    - dump_to  : Defualt is `pickle` or `json`.
    - outfile  : Defualt is None and return string. File name does not require extension.
    - indent   : Defualt is 1. Only works for json.



<h4 id="load_from_dump" class="doc_header"><code>load_from_dump</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L753" class="source_link" style="float:right">[source]</a></h4>

> <code>load_from_dump</code>(**`file_or_str`**, **`keep_as_dict`**=*`False`*)

- Loads a json/pickle dumped file or string by auto detecting it.
- **Parameters**
    - file_or_str : Filename of pickl/json or their string.
    - keep_as_dict: Defualt is False and return `Data` object. If True, returns dictionary.



<h4 id="get_file_size" class="doc_header"><code>get_file_size</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L8" class="source_link" style="float:right">[source]</a></h4>

> <code>get_file_size</code>(**`path`**)





<h4 id="interpolate_data" class="doc_header"><code>interpolate_data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L20" class="source_link" style="float:right">[source]</a></h4>

> <code>interpolate_data</code>(**`x`**, **`y`**, **`n`**=*`10`*, **`k`**=*`3`*)

- Returns interpolated xnew,ynew. If two points are same, it will add 0.1*min(dx>0) to compensate it.
- **Parameters**
    - x: 1D array of size p,
    - y: ndarray of size p*q*r,....
    - n: Number of points to add between two given points.
    - k: Polynomial order to interpolate.

- Only axis 0 will be interpolated. If you want general interploation, use `from scipy.interpolate import make_interp_spline, BSpline`

- **General Usage**: K(p),E(p,q) input from bandstructure.
    - `Knew,Enew= interpolate_data(K,E,n=10,k=3)`. cubic interploation



<h4 id="ps_to_py" class="doc_header"><code>ps_to_py</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L51" class="source_link" style="float:right">[source]</a></h4>

> <code>ps_to_py</code>(**`ps_command`**=*`'Get-ChildItem'`*, **`exec_type`**=*`'-Command'`*, **`path_to_ps`**=*`'powershell.exe'`*)

- Captures powershell output in python.
- **Parameters**
    - ps_command: enclose ps_command in ' ' or " ".
    - exec_type : type of execution, default '-Command', could be '-File'.
    - path_to_ps: path to powerhell.exe if not added to PATH variables.



<h4 id="ps_to_std" class="doc_header"><code>ps_to_std</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L85" class="source_link" style="float:right">[source]</a></h4>

> <code>ps_to_std</code>(**`ps_command`**=*`'Get-ChildItem'`*, **`exec_type`**=*`'-Command'`*, **`path_to_ps`**=*`'powershell.exe'`*)

- Prints powershell output in python std.
- **Parameters**
    - ps_command: enclose ps_command in ' ' or " ".
    - exec_type: type of execution, default '-Command', could be '-File'.
    - path_to_ps: path to powerhell.exe if not added to PATH variables.



<h4 id="select_dirs" class="doc_header"><code>select_dirs</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L102" class="source_link" style="float:right">[source]</a></h4>

> <code>select_dirs</code>(**`path`**=*`'e:\\Research\\pivotpy'`*, **`include`**=*`[]`*, **`exclude`**=*`[]`*)

- Returns selected directories recursively from a parent directory.
- **Parameters**
    - path    : path to a parent directory, default is `"."`
    - include : list of keywords to include directories, avoid wildcards.
    - exclude : list of keywords to exclude directories, avoid wildcards.
- **Returns**
    - Tuple of two elements, list of selcted directories and given path.



<h4 id="select_files" class="doc_header"><code>select_files</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L130" class="source_link" style="float:right">[source]</a></h4>

> <code>select_files</code>(**`path`**=*`'e:\\Research\\pivotpy'`*, **`include`**=*`[]`*, **`exclude`**=*`[]`*)

- Returns selected files from a given directory.
- **Parameters**
    - path    : path to a parent directory, default is `"."`
    - include : list of keywords to include files, avoid wildcards.
    - exclude : list of keywords to exclude files, avoid wildcards.
- **Returns**
    - Tuple of two elements, list of selcted files and given path.



<h4 id="get_child_items" class="doc_header"><code>get_child_items</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L154" class="source_link" style="float:right">[source]</a></h4>

> <code>get_child_items</code>(**`path`**=*`'e:\\Research\\pivotpy'`*, **`depth`**=*`None`*, **`recursive`**=*`True`*, **`include`**=*`[]`*, **`exclude`**=*`[]`*, **`filesOnly`**=*`False`*, **`dirsOnly`**=*`False`*)

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



<h4 id="invert_color" class="doc_header"><code>invert_color</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L209" class="source_link" style="float:right">[source]</a></h4>

> <code>invert_color</code>(**`color`**=*`(1, 1, 1)`*)

- Returns opposite of given complementary color.
- Input: Tuple (r,g,b).



<h4 id="printr" class="doc_header"><code>printr</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L218" class="source_link" style="float:right">[source]</a></h4>

> <code>printr</code>(**`s`**)





<h4 id="printg" class="doc_header"><code>printg</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L219" class="source_link" style="float:right">[source]</a></h4>

> <code>printg</code>(**`s`**)





<h4 id="printb" class="doc_header"><code>printb</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L220" class="source_link" style="float:right">[source]</a></h4>

> <code>printb</code>(**`s`**)





<h4 id="printy" class="doc_header"><code>printy</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L221" class="source_link" style="float:right">[source]</a></h4>

> <code>printy</code>(**`s`**)





<h4 id="printm" class="doc_header"><code>printm</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L222" class="source_link" style="float:right">[source]</a></h4>

> <code>printm</code>(**`s`**)





<h4 id="printc" class="doc_header"><code>printc</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L223" class="source_link" style="float:right">[source]</a></h4>

> <code>printc</code>(**`s`**)





<h2 id="EncodeFromNumpy" class="doc_header"><code>class</code> <code>EncodeFromNumpy</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L227" class="source_link" style="float:right">[source]</a></h2>

> <code>EncodeFromNumpy</code>(**`skipkeys`**=*`False`*, **`ensure_ascii`**=*`True`*, **`check_circular`**=*`True`*, **`allow_nan`**=*`True`*, **`sort_keys`**=*`False`*, **`indent`**=*`None`*, **`separators`**=*`None`*, **`default`**=*`None`*) :: `JSONEncoder`

- Serializes python/Numpy objects via customizing json encoder.
- **Usage**
    - `json.dumps(python_dict, cls=EncodeFromNumpy)` to get json string.
    - `json.dump(*args, cls=EncodeFromNumpy)` to create a file.json.



<h2 id="DecodeToNumpy" class="doc_header"><code>class</code> <code>DecodeToNumpy</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L255" class="source_link" style="float:right">[source]</a></h2>

> <code>DecodeToNumpy</code>(**\*`args`**, **\*\*`kwargs`**) :: `JSONDecoder`

- Deserilizes JSON object to Python/Numpy's objects.
- **Usage**
    - `json.loads(json_string,cls=DecodeToNumpy)` from string, use `json.load()` for file.



<h4 id="link_to_class" class="doc_header"><code>link_to_class</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L277" class="source_link" style="float:right">[source]</a></h4>

> <code>link_to_class</code>()

- Binds wrapper of a function to class as attribute that does exactly the same as function. Also function returned from wrapper can be used normally as well.
- **Parameters**
    - cls : A class object to which function is attached.



<h4 id="nav_links" class="doc_header"><code>nav_links</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L293" class="source_link" style="float:right">[source]</a></h4>

> <code>nav_links</code>(**`current_index`**=*`0`*, **`doc_url`**=*`'https://massgh.github.io/pivotpy/'`*, **`items`**=*`['Index', 'XmlElementTree', 'StaticPlots', 'InteractivePlots', 'Utilities', 'StructureIO', 'Widgets']`*, **`horizontal`**=*`False`*, **`out_string`**=*`False`*)





<h4 id="plot_bands" class="doc_header"><code>plot_bands</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L8" class="source_link" style="float:right">[source]</a></h4>

> <code>plot_bands</code>(**`ax`**=*`None`*, **`kpath`**=*`None`*, **`bands`**=*`None`*, **`showlegend`**=*`False`*, **`E_Fermi`**=*`None`*, **`color1`**=*`(0, 0, 0.8)`*, **`style1`**=*`'solid'`*, **`lw1`**=*`0.7`*, **`color2`**=*`(0.8, 0, 0)`*, **`style2`**=*`'dashed'`*, **`lw2`**=*`0.7`*)

- Returns axes object and plot on which all matplotlib allowed actions could be performed.
- **Parameters**
    - ax         : Matplotlib axes object, if not given, one is created.
    - kpath      : 1D array from [`get_kpts`](/pivotpy/XmlElementTree.html#get_kpts)().kpath or [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun)().kpath.
    - bands      : Dictionary Object from [`get_evals`](/pivotpy/XmlElementTree.html#get_evals) or [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun)().bands.
    - showlegend : Boolean, default is False, if true, gives legend for spin-polarized calculations.
    - E_Fermi    : If not given, automatically picked from bands object.
    - **kwargs   : lines color,width and style to distinguish spin Up and Down.
- **Returns**
    - ax : matplotlib axes object with plotted bands.



<h4 id="modify_axes" class="doc_header"><code>modify_axes</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L60" class="source_link" style="float:right">[source]</a></h4>

> <code>modify_axes</code>(**`ax`**=*`None`*, **`xticks`**=*`[]`*, **`xt_labels`**=*`[]`*, **`xlim`**=*`[]`*, **`yticks`**=*`[]`*, **`yt_labels`**=*`[]`*, **`ylim`**=*`[]`*, **`xlabel`**=*`None`*, **`ylabel`**=*`None`*, **`vlines`**=*`True`*, **`zeroline`**=*`True`*)

- Returns None, applies given settings on axes. Prefered to use before other plotting.
- **Parameters**
    - ax  : Matplotlib axes object.
    - (x,y)ticks : List of positions on (x,y axes).
    - (xt,yt)_labels : List of labels on (x,y) ticks points.
    - (x,y)lim : [min, max] of (x,y) axes.
    - (x,y)label : axes labels.
    - vlines : If true, drawn when `ylim` is not empty.
    - zeroline : If True, drawn when `xlim` is not empty.



<h4 id="quick_bplot" class="doc_header"><code>quick_bplot</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L105" class="source_link" style="float:right">[source]</a></h4>

> <code>quick_bplot</code>(**`path_evr`**=*`None`*, **`ax`**=*`None`*, **`skipk`**=*`None`*, **`joinPathAt`**=*`[]`*, **`elim`**=*`[]`*, **`xt_indices`**=*`[]`*, **`xt_labels`**=*`[]`*, **`E_Fermi`**=*`None`*, **`figsize`**=*`(3.4, 2.6)`*, **`txt`**=*`None`*, **`xytxt`**=*`[0.05, 0.9]`*, **`ctxt`**=*`'black'`*)

- Returns axes object and plot on which all matplotlib allowed actions could be performed.
- **Parameters**
    - path_evr   : path/to/vasprun.xml or output of [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun). Auto picks in CWD.
    - ax         : Matplotlib axes object, if not given, one is created.
    - skipk      : Number of kpoints to skip, default will be from IBZKPT.
    - joinPathAt : Points where kpath is broken.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun).
    - xt_indices : High symmetry kpoints indices.abs
    - xt_labels  : High Symmetry kpoints labels.
    - **kwargs   : figsize=(3.4,2.6). Text,its position and color.
- **Returns**
    - ax : matplotlib axes object with plotted bands.



<h4 id="add_text" class="doc_header"><code>add_text</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L168" class="source_link" style="float:right">[source]</a></h4>

> <code>add_text</code>(**`ax`**=*`None`*, **`xs`**=*`0.05`*, **`ys`**=*`0.9`*, **`txts`**=*`'[List]'`*, **`colors`**=*`'r'`*)

- Adds text entries on axes, given single string or list.
- **Parameters**
    - xs    : List of x coordinates relative to axes or single coordinate.
    - ys    : List of y coordinates relative to axes or single coordinate.
    - txts  : List of strings or one string.
    - colors: List of x colors of txts or one color.



<h4 id="add_legend" class="doc_header"><code>add_legend</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L193" class="source_link" style="float:right">[source]</a></h4>

> <code>add_legend</code>(**`ax`**=*`None`*, **`colors`**=*`[]`*, **`labels`**=*`[]`*, **`styles`**=*`'solid'`*, **`widths`**=*`0.7`*, **`anchor`**=*`(0, 1)`*, **`ncol`**=*`3`*, **`loc`**=*`'lower left'`*, **`fontsize`**=*`'small'`*, **`frameon`**=*`False`*, **\*\*`legend_kwargs`**)

- Adds custom legeneds on a given axes,returns None.
- **Parameters**
    - ax       : Matplotlib axes.
    - colors   : List of colors.
    - labels   : List of labels.
    - styles   : str or list of line styles.
    - widths   : str or list of line widths.
    - **kwargs : Matplotlib's legend arguments.



<h4 id="add_colorbar" class="doc_header"><code>add_colorbar</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L232" class="source_link" style="float:right">[source]</a></h4>

> <code>add_colorbar</code>(**`ax`**=*`None`*, **`colors`**=*`[]`*, **`n`**=*`20`*, **`ticks`**=*`[20, 60, 100]`*, **`ticklabels`**=*`['s', 'p', 'd']`*, **`linewidth`**=*`None`*, **`vertical`**=*`False`*, **`fontsize`**=*`8`*)

- Plots colorbar on a given axes. This axes should be only for colorbar. Returns None.
- **Parameters**
    - ax         : Matplotlib axes object.
    - colors     : List of colors in colorbar, if not given, RGB colorbar is added.
    - vertical   : Boolean, default is Fasle.
    - n          : int, number of points between colors. Default 20.
    - ticks      : List of tick points to show on colorbar.
    - ticklabels : List of labels for ticks.
    - linewidth  : To tweek in order to make smooth gradient.
    - fontsize   : Default 8. Adjustable according to plot space.



<h4 id="create_rgb_lines" class="doc_header"><code>create_rgb_lines</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L297" class="source_link" style="float:right">[source]</a></h4>

> <code>create_rgb_lines</code>(**`ax`**=*`None`*, **`kpath`**=*`None`*, **`evals_set`**=*`None`*, **`pros_set`**=*`None`*, **`ions`**=*`[0]`*, **`orbs`**=*`[[0], [], []]`*, **`labels`**=*`['', '', '']`*, **`uni_width`**=*`False`*, **`max_width`**=*`2.5`*, **`uni_color`**=*`False`*, **`color`**=*`'red'`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*, **`scale_color`**=*`False`*)

- Plot on a given axes or returns line collection, lines and colors if axes is None, which can be added to an ax only onces,by using `ax.add_collection(collection)` and then `ax.autoscale_view()` will make it visible.
- **Parameters**
    - ax       : Matplotlib axes object, if not given, linecollection is returned.
    - kapath   : [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun)().kpath or [`get_kpts`](/pivotpy/XmlElementTree.html#get_kpts)().kpath.
    - evals_set: [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun)().bands.evals or [`get_evals`](/pivotpy/XmlElementTree.html#get_evals)().evals. If calculations are spin-polarized, it will be `...evals.SpinUp/SpinDown` for both. You need to create collections twice for SpinUp and SpinDown separately.
    - pros_set : `export_vasprun().pro_bands.pros` or [`get_bands_pro_set`](/pivotpy/XmlElementTree.html#get_bands_pro_set)().pros. If calculations are spin-polarized, it will be `...pros.SpinUp/SpinDown` for both. You need to create collections twice for SpinUp and SpinDown separately.
    - ions     : List of ions to project on, could be `range(start,stop,step)` as well, remember that `stop` is not included in python. so `range(0,2)` will generate 0 and 1 indices.
    - orbs     : List of three lists of orbitals indices. `[[red],[green],[blue]]`, you can create any color by this combination. For example, to get `s-orbital in yellow color`, you will use `[[0],[0],[]]`. Do not remove empty list from there, it will not effect your orbital selection.
    - uni_width: If True, will keep equal `width=max_width/2` of lines.
    - max_width: Default is 5. Orbitals' projections are added and Normalized to this thickness.
    - uni_color: If True, will not change color in a band from point to point,width is reduced.
    - color    : (str,rgb,rgba), if `uni_color=True`, color will be applied to line.
    - interpolate: Deafult is false, if True, it will add n points between nearest kpoints.
    - n        : int, default is 5. Adds n points between nearest kpoints.
    - k        : int, order of interpolation, defualt is 3. `n > k` should be hold.
    - scale_color: If True, colors are scaled to 1 at each points.
- **Returns**
    - line collection : Matplotlib line collection object.
    - line patches    : An (N,2,2) dimensional arry.
    - colors          : An (N,4) or (N,3) dimensional list.(Not scaled.)
- **Exception**
    - If `uni_color` and `uni_width` are True together, this leads to simple plot. No collections will be created. Use `bands_plot()` instead.



<h4 id="quick_rgb_lines" class="doc_header"><code>quick_rgb_lines</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L398" class="source_link" style="float:right">[source]</a></h4>

> <code>quick_rgb_lines</code>(**`path_evr`**=*`None`*, **`ax`**=*`None`*, **`skipk`**=*`None`*, **`joinPathAt`**=*`[]`*, **`elim`**=*`[]`*, **`elements`**=*`[[0], [], []]`*, **`orbs`**=*`[[0], [], []]`*, **`labels`**=*`['Elem0-s', '', '']`*, **`max_width`**=*`2.5`*, **`xt_indices`**=*`[0, -1]`*, **`xt_labels`**=*`['$\\Gamma$', 'M']`*, **`E_Fermi`**=*`None`*, **`figsize`**=*`(3.4, 2.6)`*, **`txt`**=*`None`*, **`xytxt`**=*`[0.05, 0.9]`*, **`ctxt`**=*`'black'`*, **`uni_width`**=*`False`*, **`interpolate`**=*`False`*, **`spin`**=*`'both'`*, **`n`**=*`5`*, **`k`**=*`3`*, **`scale_color`**=*`True`*, **`colorbar`**=*`True`*)

- Returns axes object and plot on which all matplotlib allowed actions could be performed. In this function,orbs,labels,elements all have list of length 3. Inside list, sublists or strings could be any length but should be there even if empty.
- **Parameters**
    - path_evr   : path/to/vasprun.xml or output of [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun). Auto picks in CWD.
    - ax         : Matplotlib axes object, if not given, one is created.
    - skipk      : Number of kpoints to skip, default will be from IBZKPT.
    - joinPathAt : Points where kpath is broken.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun).
    - xt_indices : High symmetry kpoints indices.abs
    - xt_labels  : High Symmetry kpoints labels.
    - elements   : List [[0],[],[]] by default and plots s orbital of first ion..
    - orbs       : List [[r],[g],[b]] of indices of orbitals, could be empty, but shape should be same.
    - labels     : List [str,str,str] of projection labels. empty string should exist to maintain shape. Auto adds `↑`,`↓` for ISPIN=2.
    - max_width  : Width to scale whole projections. if `uni_width=True, width=max_width/2`.
    - figsize    : Tuple (width,height) in inches. Default (3.4.2.6) is article column's width.
    - txt        : Text on figure, if None, SYSTEM's name is printed.
    - xytxt      : [x_coord,y_coord] of text relative to axes.
    - ctxt       : color of text.
    - uni_width  : If True, width of bands kept uniform.
    - uni_color  : If True, color of bands kept same.
    - color      : (str,rgb,rgba), if `uni_color=True`, color is applied.
    - spin       : Plot spin-polarized for spin {'up','down','both'}. Default is both.
    - interpolate: Default is False, if True, bands are interpolated.
    - n          : int, number of points, default is 5.
    - k          : int, order of interpolation 0,1,2,3. Defualt 3. `n > k` should be hold.
    - scale_color: Boolean. Default True, colors are scaled to 1 at each point.
    - colorbar   : Default is True. Displays a vertical RGB colorbar.
- **Returns**
    - ax : matplotlib axes object with plotted projected bands.



<h4 id="quick_color_lines" class="doc_header"><code>quick_color_lines</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L586" class="source_link" style="float:right">[source]</a></h4>

> <code>quick_color_lines</code>(**`path_evr`**=*`None`*, **`axes`**=*`None`*, **`skipk`**=*`None`*, **`joinPathAt`**=*`[]`*, **`elim`**=*`[]`*, **`elements`**=*`[[0]]`*, **`orbs`**=*`[[0]]`*, **`labels`**=*`['s']`*, **`color_map`**=*`'gist_rainbow'`*, **`max_width`**=*`2.5`*, **`xt_indices`**=*`[0, -1]`*, **`xt_labels`**=*`['$\\Gamma$', 'M']`*, **`E_Fermi`**=*`None`*, **`showlegend`**=*`True`*, **`figsize`**=*`(3.4, 2.6)`*, **`txt`**=*`None`*, **`xytxt`**=*`[0.05, 0.85]`*, **`ctxt`**=*`'black'`*, **`spin`**=*`'both'`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*, **`legend_kwargs`**=*`{'ncol': 4, 'anchor': (0, 0.85), 'handletextpad': 0.5, 'handlelength': 1, 'fontsize': 'small', 'frameon': True}`*, **\*\*`subplots_adjust_kwargs`**)

- Returns axes object and plot on which all matplotlib allowed actions could be performed. If given, axes,elements,orbs colors, and labels must have same length. If not given, zeroth ion is plotted with s-orbital.
- **Parameters**
    - path_evr   : Path/to/vasprun.xml or output of [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun). Auto picks in CWD.
    - axes       : Matplotlib axes object with one or many axes, if not given, auto created.
    - skipk      : Number of kpoints to skip, default will be from IBZKPT.
    - joinPathAt : Points where kpath is broken.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun).
    - xt_indices : High symmetry kpoints indices.abs
    - xt_labels  : High Symmetry kpoints labels.
    - elements   : List [[0],], by defualt and plot first ion's projections.
    - orbs       : List [[0],] lists of indices of orbitals, could be empty.
    - labels     : List [str,] of orbitals labels. len(labels)==len(orbs) must hold.  Auto adds `↑`,`↓` for ISPIN=2.
    - color_map  : Matplotlib's standard color maps. Default is 'gist_ranibow'.
    - showlegend : True by defualt.
    - max_width  : Width to scale whole projections. if `uni_width=True, width=max_width/2`.
    - figsize    : Tuple (width,height) in inches. Default (3.4.2.6) is article column's width.
    - txt        : Text on figure, if None, SYSTEM's name is printed.
    - xytxt      : [x_coord,y_coord] of text relative to axes.
    - ctxt       : color of text.
    - spin       : Plot spin-polarized for spin {'up','down','both'}. Default is both.
    - interpolate: Default is False, if True, bands are interpolated.
    - n          : int, number of points, default is 5.
    - k          : int, order of interpolation 0,1,2,3. Defualt 3. `n > k` should be hold.
    - legend_kwargs: Dictionary to contain legend arguments to fix.
    - **subplots_adjust_kwargs : plt.subplots_adjust parameters.
- **Returns**
    - ax : matplotlib axes object with plotted projected bands.



<h4 id="init_figure" class="doc_header"><code>init_figure</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L771" class="source_link" style="float:right">[source]</a></h4>

> <code>init_figure</code>(**`figsize`**=*`(3.4, 2.6)`*, **`nrows`**=*`1`*, **`ncols`**=*`1`*, **`widths`**=*`[]`*, **`heights`**=*`[]`*, **`axes_off`**=*`[]`*, **`sharex`**=*`False`*, **`sharey`**=*`False`*, **\*\*`subplots_adjust_kwargs`**)

- Returns all axes of initialized figure, based on plt.subplots().
- **Parameters**
    - figsize   : Tuple (width, height). Default is (3.4,2.6).
    - nrows     : Default 1.
    - ncols     : Default 1.
    - widths    : List with len(widths)==nrows, to set width ratios of subplots.
    - heights   : List with len(heights)==ncols, to set height ratios of subplots.
    - share(x,y): Share axes between plots, this removes shared ticks automatically.
    - axes_off  : Turn off axes visibility, If `nrows = ncols = 1, set True/False`, If anyone of `nrows or ncols > 1`, provide list of axes indices to turn off. If both `nrows and ncols > 1`, provide list of tuples (x_index,y_index) of axes.
    - **subplots_adjust_kwargs : These are same as `plt.subplots_adjust()`'s arguements.



<h4 id="select_pdos" class="doc_header"><code>select_pdos</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L835" class="source_link" style="float:right">[source]</a></h4>

> <code>select_pdos</code>(**`tdos`**=*`None`*, **`pdos_set`**=*`None`*, **`ions`**=*`[0]`*, **`orbs`**=*`[0]`*, **`E_Fermi`**=*`0`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*)

- Returns (interpolated/orginal) enrgy(N,), tdos(N,), and pdos(N,) of selected ions/orbitals.
- **Parameters**
    - tdos     : [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun)().tdos or [`get_tdos`](/pivotpy/XmlElementTree.html#get_tdos)().tdos. If calculations are spin-polarized, it will be `..tdos.SpinUp/SpinDown` for both. You need to apply this function twice for SpinUp and SpinDown separately.
    - pdos_set : `export_vasprun().pro_dos.pros` or [`get_dos_pro_set`](/pivotpy/XmlElementTree.html#get_dos_pro_set)().pros. If calculations are spin-polarized, it will be `...pros.SpinUp/SpinDown` for both.
    - ions     : List of ions to project on, could be `range(start,stop,step)` as well, remember that `stop` is not included in python. so `range(0,2)` will generate 0 and 1 indices.
    - orbs     : List of orbitals indices to pick.
    - E_Fermi  : Here it is zero. Needs to be input.
    - interpolate: Deafult is false, if True, it will add n points between nearest points.
    - n        : int, default is 5. Adds n points between nearest kpoints.
    - k        : int, order of interpolation, defualt is 3. `n > k` should be hold.



<h4 id="collect_dos" class="doc_header"><code>collect_dos</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L875" class="source_link" style="float:right">[source]</a></h4>

> <code>collect_dos</code>(**`path_evr`**=*`None`*, **`elim`**=*`[]`*, **`elements`**=*`[[0]]`*, **`orbs`**=*`[[0]]`*, **`labels`**=*`['s']`*, **`E_Fermi`**=*`None`*, **`spin`**=*`'both'`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*)

- Returns lists of energy,tdos, pdos and labels. If given,elements,orbs and labels must have same length. If not given, zeroth ions is collected with s-orbital.
- **Parameters**)
    - path_evr   : Path/to/vasprun.xml or output of [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun). Auto picks in CWD.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun).
    - elements   : List [[0],], by defualt and plot first ion's projections.
    - orbs       : List [[0],] lists of indices of orbitals, could be empty.
    - labels     : List [str,] of orbitals labels. len(labels)==len(orbs) must hold.  Auto adds `↑`,`↓` for ISPIN=2.
    - spin       : Plot spin-polarized for spin {'up','down','both'}. Default is both.
    - interpolate: Default is False, if True, bands are interpolated.
    - n          : int, number of points, default is 5.
    - k          : int, order of interpolation 0,1,2,3. Defualt 3. `n > k` should be hold.
- **Returns**
    - Energy : (N,1) size.
    - tdos   : (N,1) size or [(N,1),(N,1)] if spin polarized.
    - pdos   : [(N,1),(N,1),...], spin polarized is auto-fixed.
    - labels : ['label1,'label2',...] spin polarized is auto-fixed.
    - vr     : Exported vasprun.



<h4 id="quick_dos_lines" class="doc_header"><code>quick_dos_lines</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L991" class="source_link" style="float:right">[source]</a></h4>

> <code>quick_dos_lines</code>(**`path_evr`**=*`None`*, **`ax`**=*`None`*, **`elim`**=*`[]`*, **`include_dos`**=*`'both'`*, **`elements`**=*`[[0]]`*, **`orbs`**=*`[[0]]`*, **`labels`**=*`['s']`*, **`color_map`**=*`'gist_rainbow'`*, **`tdos_color`**=*`(0.8, 0.95, 0.8)`*, **`linewidth`**=*`0.5`*, **`fill_area`**=*`True`*, **`vertical`**=*`False`*, **`E_Fermi`**=*`None`*, **`figsize`**=*`(3.4, 2.6)`*, **`txt`**=*`None`*, **`xytxt`**=*`[0.05, 0.85]`*, **`ctxt`**=*`'black'`*, **`spin`**=*`'both'`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*, **`showlegend`**=*`True`*, **`legend_kwargs`**=*`{'ncol': 4, 'anchor': (0, 1), 'handletextpad': 0.5, 'handlelength': 1, 'fontsize': 'small', 'frameon': True}`*)

- Returns ax object (if ax!=False) and plot on which all matplotlib allowed actions could be performed, returns lists of energy,tdos and pdos and labels. If given,elements,orbs colors, and labels must have same length. If not given, zeroth ions is plotted with s-orbital.
- **Parameters**)
    - path_evr   : Path/to/vasprun.xml or output of [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun). Auto picks in CWD.
    - ax         : Matplotlib axes object, if None, one is created. If False, data lists are returned.
    - include_dos: One of {'both','tdos','pdos'}.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun).
    - elements   : List [[0],], by defualt and plot first ion's projections.
    - orbs       : List [[0],] lists of indices of orbitals, could be empty.
    - labels     : List [str,] of orbitals labels. len(labels)==len(orbs) must hold.  Auto adds `↑`,`↓` for ISPIN=2.
    - color_map  : Matplotlib's standard color maps. Default is 'gist_ranibow'. Use 'RGB' if want to compare with [`quick_rgb_lines`](/pivotpy/StaticPlots.html#quick_rgb_lines) with 3 projection inputs (len(orbs)==3).
    - fill_area  : Default is True and plots filled area for dos. If False, plots lines only.
    - vertical   : False, If True, plots along y-axis.
    - showlegend : True by defualt.
    - figsize    : Tuple (width,height) in inches. Default (3.4.2.6) is article column's width.
    - txt        : Text on figure, if None, SYSTEM's name is printed.
    - xytxt      : [x_coord,y_coord] of text relative to axes.
    - ctxt       : color of text.
    - spin       : Plot spin-polarized for spin {'up','down','both'}. Default is both.
    - interpolate: Default is False, if True, bands are interpolated.
    - n          : int, number of points, default is 5.
    - k          : int, order of interpolation 0,1,2,3. Defualt 3. `n > k` should be hold.
    - legend_kwargs: Dictionary to contain legend arguments to fix.
- **Returns**
    - ax         : Matplotlib axes.



<h4 id="plt_to_html" class="doc_header"><code>plt_to_html</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L1145" class="source_link" style="float:right">[source]</a></h4>

> <code>plt_to_html</code>(**`plt_fig`**=*`None`*, **`dpi`**=*`300`*, **`dash_html`**=*`None`*)

- Returns base64 encoded Image to display in notebook or HTML <img/> or plotly's dash_html_components.Img object.
- **Parameters**
    - plt_fig  : Matplotlib's figure instance, auto picks as well.
    - dpi      : PNG images's DPI, default is 300.
    - dash_html: Default is None which results in an image display in jupyter notebook.
        - If True, returns html.Img object for plotly's dash.
        - If False, returns html <img/> object to embed in HTML DOM.



<h4 id="get_rgb_data" class="doc_header"><code>get_rgb_data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L7" class="source_link" style="float:right">[source]</a></h4>

> <code>get_rgb_data</code>(**`kpath`**=*`None`*, **`evals_set`**=*`None`*, **`pros_set`**=*`None`*, **`elements`**=*`[[0], [], []]`*, **`orbs`**=*`[[0], [], []]`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*, **`scale_color`**=*`False`*)

- Returns a formatted RGB colored data to pass into [`rgb_to_plotly`](/pivotpy/InteractivePlots.html#rgb_to_plotly) function. Two arguments, `elements` and `orbs` should be in one-to-one correspondence. Returned item has transpose data shape, so that main iteration is over bands.
- **Parameters**
    - ax       : Matplotlib axes object, if not given, linecollection is returned.
    - kapath   : [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun)().kpath or [`get_kpts`](/pivotpy/XmlElementTree.html#get_kpts)().kpath.
    - evals_set: [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun)().bands.evals or [`get_evals`](/pivotpy/XmlElementTree.html#get_evals)().evals. If calculations are spin-polarized, it will be `...evals.SpinUp/SpinDown` for both. You need to apply twice for SpinUp and SpinDown separately.
    - pros_set : `export_vasprun().pro_bands.pros` or [`get_bands_pro_set`](/pivotpy/XmlElementTree.html#get_bands_pro_set)().pros. If calculations are spin-polarized, it will be `...pros.SpinUp/SpinDown` for both. You need to create collections twice for SpinUp and SpinDown separately.
    - elements : List three lists of ions to project on, each element could be `range(start,stop,step)` as well, remember that `stop` is not included in python. so `range(0,2)` will generate 0 and 1 indices.
    - orbs     : List of three lists of orbitals indices. `[[red],[green],[blue]]`, you can create any color by this combination. For example, to get `s-orbital in yellow color`, you will use `[[0],[0],[]]`. Do not remove empty list from there, it will not effect your orbital selection.
    - interpolate: Deafult is false, if True, it will add n points between nearest kpoints.
    - n        : int, default is 5. Adds n points between nearest kpoints.
    - k        : int, order of interpolation, defualt is 3. `n > k` should be hold.
    - scale_color: If True, colors are scaled to 1 at each points.
- **Returns**
    - kpath    : List of NKPTS. (interpolated if given.)
    - evals    : An (NBAND,NKPTS) numpy arry.
    - colors   : An (NBANDS,NKPTS,3) numpy array.



<h4 id="flip_even_patches" class="doc_header"><code>flip_even_patches</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L74" class="source_link" style="float:right">[source]</a></h4>

> <code>flip_even_patches</code>(**`array_1d`**, **`patch_length`**)

- When you reshape bands data to 1D array, you may need to draw lines which do not link ends of plot, for that, it is required to flip patches, so that next band start from where 1st end and so one.
- **Parameters**
    - array_1d     : Numpy 1d array or list.
    - patch_length : length of xaxis patches, e.g NKPTS.
- **Returns**
    - 1D list
- ** Example**
    > k=[1,2,3,1,2,3]
    > flip_even_patches(k,3)
    > [1,2,3,3,2,1]



<h4 id="rgb_to_plotly" class="doc_header"><code>rgb_to_plotly</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L98" class="source_link" style="float:right">[source]</a></h4>

> <code>rgb_to_plotly</code>(**`rgb_data`**=*`None`*, **`mode`**=*`'markers'`*, **`max_width`**=*`5`*, **`showlegend`**=*`False`*, **`name`**=*`''`*, **`labels`**=*`['s', 'p', 'd']`*, **`symbol`**=*`0`*)

- Returns data object of plotly's figure using [`get_rgb_data`](/pivotpy/InteractivePlots.html#get_rgb_data). Returned data could be fed to a plolty's figure.
- ** Parameters**
    - rgb_data    : output of [`get_rgb_data`](/pivotpy/InteractivePlots.html#get_rgb_data).
    - mode        : Three plotting modes are available:
        - 'markers' : Plot whole data as a single scatter object. Its too fast.
        - 'bands'   : Plot data such that each band is accessible via legend.
        - 'lines'   : A replica of `matplotlib LineCollection` object. It plots at each point separately, slower than other two modes.
    - max_width  : Line/Scatter thickness is scaled to `max_width`.
    - name       : Name to be shown on hover text or legend.
    - labels     : Optional, show red green blue colors corresponding orbitals.
    - showlegend : Optional, only suitbale if spin up/down or 'bands' mode is ON.
    - symbol     : Plotly's marker symbol. 0 for circle, 5/6 for Up/Down.



<h4 id="plotly_to_html" class="doc_header"><code>plotly_to_html</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L163" class="source_link" style="float:right">[source]</a></h4>

> <code>plotly_to_html</code>(**`fig`**, **`filename`**=*`None`*, **`out_string`**=*`False`*, **`modebar`**=*`True`*)

- Writes plotly's figure as HTML file or display in IPython which is accessible when online. It is different than plotly's `fig.to_html` as it is minimal in memory. If you need to have offline working file, just use `fig.write_html('file.html')` which will be larger in size.
- **Parameters**
    - fig      : A plotly's figure object.
    - filename : Name of file to save fig. Defualt is None and show plot in Colab/Online or return hrml string.
    - out_string: If True, returns HTML string, if False displays graph if possible.



<h4 id="plotly_rgb_lines" class="doc_header"><code>plotly_rgb_lines</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L215" class="source_link" style="float:right">[source]</a></h4>

> <code>plotly_rgb_lines</code>(**`path_evr`**=*`None`*, **`elements`**=*`[[], [], []]`*, **`orbs`**=*`[[], [], []]`*, **`labels`**=*`['', '', '']`*, **`mode`**=*`'markers'`*, **`elim`**=*`[]`*, **`E_Fermi`**=*`None`*, **`skipk`**=*`None`*, **`joinPathAt`**=*`[]`*, **`max_width`**=*`6`*, **`title`**=*`None`*, **`xt_indices`**=*`[0, -1]`*, **`xt_labels`**=*`['Γ', 'M']`*, **`figsize`**=*`None`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*)

- Returns plotly's figure object, takes care of spin-polarized calculations automatically. `elements`,`orbs` and `labels` are required to be one-to-one lists of size 3 where each item in list could be another list or integer.
- **Parameters**
    - path_ever  : Path/to/vasprun.xml or xml output of [`read_asxml`](/pivotpy/XmlElementTree.html#read_asxml).
    - elements   : List of size 3 of list of indices of ions. If not given, picks all ions for each orbital.
    - orbs       : List of size 3 of list of orbital indices, if not gievn, s,p,d plotted.
    - labels  : List of labels for projection.
    - mode       : Three plotting modes are available:
        - 'markers' : Plot whole data as a single scatter object. Its too fast.
        - 'bands'   : Plot data such that each band is accessible via legend.
        - 'lines'   : A replica of `matplotlib LineCollection` object. It plots at each point separately, slower than other two modes.
    - **kwargs      : interpolate, ticks, figsize,elim,joinPathAt,max_width,title etc.



<h4 id="plotly_dos_lines" class="doc_header"><code>plotly_dos_lines</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L369" class="source_link" style="float:right">[source]</a></h4>

> <code>plotly_dos_lines</code>(**`path_evr`**=*`None`*, **`elim`**=*`[]`*, **`elements`**=*`[[0]]`*, **`orbs`**=*`[[0]]`*, **`labels`**=*`['s']`*, **`color_map`**=*`'gist_rainbow'`*, **`tdos_color`**=*`(0.5, 0.95, 0)`*, **`linewidth`**=*`2`*, **`fill_area`**=*`True`*, **`vertical`**=*`False`*, **`E_Fermi`**=*`None`*, **`figsize`**=*`None`*, **`spin`**=*`'both'`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*, **`title`**=*`None`*)

- Returns ax object (if ax!=False) and plot on which all matplotlib allowed actions could be performed, returns lists of energy,tdos and pdos and labels. If given,elements,orbs colors, and labels must have same length. If not given, zeroth ions is plotted with s-orbital.
- **Parameters**)
    - path_evr   : Path/to/vasprun.xml or output of [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun). Auto picks in CWD.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from [`export_vasprun`](/pivotpy/XmlElementTree.html#export_vasprun).
    - elements   : List [[0,],] of ions indices, by defualt plot first ion's projections.
    - orbs       : List [[0,],] lists of indices of orbitals, could be empty.
    - labels     : List [str,] of orbitals labels. len(labels)==len(orbs) must hold.
    - color_map  : Matplotlib's standard color maps. Default is 'gist_ranibow'. Use 'RGB' if want to compare with [`plotly_rgb_lines`](/pivotpy/InteractivePlots.html#plotly_rgb_lines) with 3 projection inputs (len(orbs)==3).
    - fill_area  : Default is True and plots filled area for dos. If False, plots lines only.
    - vertical   : False, If True, plots along y-axis.
    - figsize    : Tuple in pixels (width,height).
    - interpolate: Default is False, if True, bands are interpolated.
    - n          : int, number of points, default is 5.
    - k          : int, order of interpolation 0,1,2,3. Defualt 3. `n > k` should be hold.
    - legend_kwargs: Dictionary to contain legend arguments to fix.
- **Returns**
    - fig        : Plotly's figure object.



<h4 id="save_mp_API" class="doc_header"><code>save_mp_API</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L8" class="source_link" style="float:right">[source]</a></h4>

> <code>save_mp_API</code>(**`api_key`**)

- Save materials project api key for autoload in functions.



<h4 id="load_mp_data" class="doc_header"><code>load_mp_data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L29" class="source_link" style="float:right">[source]</a></h4>

> <code>load_mp_data</code>(**`formula`**, **`api_key`**=*`None`*, **`mp_id`**=*`None`*, **`max_sites`**=*`None`*)

- Returns fetched data using request api of python form materials project website.
- **Parameters**
    - formula  : Material formula such as 'NaCl'.
    - api_key  : API key for your account from material project site. Auto picks if you already used [`save_mp_API`](/pivotpy/StructureIO.html#save_mp_API) function.
    - mp_id    : Optional, you can specify material ID to filter results.
    -max_sites : Option, you can set maximum number of sites to load fastly as it will not fetch all large data sets.



<h4 id="get_crystal" class="doc_header"><code>get_crystal</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L72" class="source_link" style="float:right">[source]</a></h4>

> <code>get_crystal</code>(**`formula`**, **`api_key`**=*`None`*, **`mp_id`**=*`None`*, **`max_sites`**=*`None`*)

- Returns crystal information dictionary including cif data format.
- **Parameters**
    - formula  : Material formula such as 'NaCl'.
    - api_key  : API key for your account from material project site. Auto picks if you already used [`save_mp_API`](/pivotpy/StructureIO.html#save_mp_API) function.
    - mp_id    : Optional, you can specify material ID to filter results.
    -max_sites : Option, you can set maximum number of sites to load fastly as it will not fetch all large data sets.



<h4 id="get_poscar" class="doc_header"><code>get_poscar</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L96" class="source_link" style="float:right">[source]</a></h4>

> <code>get_poscar</code>(**`formula`**, **`api_key`**=*`None`*, **`mp_id`**=*`None`*, **`max_sites`**=*`None`*)

- Returns poscar information dictionary including cif data format.
- **Parameters**
    - formula  : Material formula such as 'NaCl'.
    - api_key  : API key for your account from material project site. Auto picks if you already used [`save_mp_API`](/pivotpy/StructureIO.html#save_mp_API) function.
    - mp_id    : Optional, you can specify material ID to filter results.
    -max_sites : Option, you can set maximum number of sites to load fastly as it will not fetch all large data sets.
- **Usage**
    - `get_poscar('GaAs',api_key,**kwargs)`. Same result is returned from `Get-POSCAR` command in PowerShell terminal if Vasp2Visual module is installed.



<h4 id="get_kpath" class="doc_header"><code>get_kpath</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L168" class="source_link" style="float:right">[source]</a></h4>

> <code>get_kpath</code>(**`hsk_list`**=*`[]`*, **`labels`**=*`[]`*, **`n`**=*`5`*, **`weight`**=*`None`*, **`ibzkpt`**=*`None`*, **`outfile`**=*`None`*)

- Generate list of kpoints along high symmetry path. Options are write to file or return KPOINTS list. It generates uniformly spaced point with input `n` as just a scale factor of number of points per unit length. You can also specify custom number of kpoints in an interval by putting number of kpoints as 4th entry in left kpoint.
- **Parameters**
    - hsk_list : N x 3 list of N high symmetry points, if broken path then [[N x 3],[M x 3],...]. Optionally you can put a 4 values point where 4th entry will decide number of kpoints in current interval. Make sure that points in a connected path patch are at least two i.e. `[[x1,y1,z1],[x2,y2,z2]]` or `[[x1,y1,z1,N],[x2,y2,z2]]`.
    - n        ; int, number per unit length, this makes uniform steps based on distance between points.
    - weight : Float, if None, auto generates weights.
    - gamma  : If True, shifts mesh at gamma point.
    - ibzkpt : Path to ibzkpt file, required for HSE calculations.
    - labels : Hight symmetry points labels. Good for keeping record of lables and points indices for later use.                - Note: If you do not want to label a point, label it as 'skip' at its index and it will be removed.
    - outfile: Path/to/file to write kpoints.
- **Attributes**
    - If `outfile = None`, a tuple is returned which consists of:
        - nkpts   : get_kmesh().nkpts.
        - kpoints : get_kmesh().kpoints.
        - weights : get_kmesh().weights.



<h4 id="get_kmesh" class="doc_header"><code>get_kmesh</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L240" class="source_link" style="float:right">[source]</a></h4>

> <code>get_kmesh</code>(**`n_xyz`**=*`[5, 5, 5]`*, **`weight`**=*`None`*, **`gamma`**=*`True`*, **`ibzkpt`**=*`None`*, **`poscar`**=*`None`*, **`outfile`**=*`None`*, **`plot`**=*`False`*)

- Generate uniform mesh of kpoints. Options are write to file, plot or return KPOINTS list.
- **Parameters**
    - n_xyz  : List of [nx ny nz] or integer. If integere given, kmesh is autoscaled.
    - weight : Float, if None, auto generates weights.
    - gamma  : Default True, shifts mesh at gamma point.
    - ibzkpt : Path to ibzkpt file, required for HSE calculations.
    - poscar : POSCAR file or real space lattice vectors, if None, cubic symmetry is used and it is fast.
    - outfile: Path/to/file to write kpoints.
    - plot   : If True, returns interactive plot. You can look at mesh before you start calculation.
- **Attributes**
    - If `plot = False`, a tuple is returned which consists of:
        - nkpts   : get_kmesh().nkpts.
        - kpoints : get_kmesh().kpoints.
        - weight : get_kmesh().weight, its one float number, provided or calculated.



<h4 id="intersect_3p_p_3v" class="doc_header"><code>intersect_3p_p_3v</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L360" class="source_link" style="float:right">[source]</a></h4>

> <code>intersect_3p_p_3v</code>(**`a`**, **`b`**, **`c`**)

- Returns intersection point of 3 planes in 3D.
- **Parameters**
    - a,b,c : three vectors in 3D, whose perpendicular planes will be intersected.



<h4 id="centroid" class="doc_header"><code>centroid</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L373" class="source_link" style="float:right">[source]</a></h4>

> <code>centroid</code>(**`points`**)

- Returns centroid of a list of 3D points.
- **Parameters**
    - points: List[List(len=3)]



<h4 id="order" class="doc_header"><code>order</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L390" class="source_link" style="float:right">[source]</a></h4>

> <code>order</code>(**`points`**)

- Returns counterclockwise ordered vertices of a plane in 3D. Append first vertex at end to make loop.
- **Parameters**
    - points: List[List(len=3)]



<h4 id="in_vol_sector" class="doc_header"><code>in_vol_sector</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L438" class="source_link" style="float:right">[source]</a></h4>

> <code>in_vol_sector</code>(**`test_point`**, **`p1`**, **`p2`**, **`p3`**)

- Returns True if test_point lies inside/on the overlapping planes of three vectors.
- **Parameters**
    - p1,p2,p3: Three vectors points in 3D.



<h4 id="out_bz_plane" class="doc_header"><code>out_bz_plane</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L460" class="source_link" style="float:right">[source]</a></h4>

> <code>out_bz_plane</code>(**`test_point`**, **`plane`**)

- Returns True if test_point is between plane and origin. Could be used to sample BZ mesh in place of ConvexHull.
- **Parameters**
    - test_points: 3D point.
    - plane      : List of at least three coplanar points.



<h4 id="to_xy" class="doc_header"><code>to_xy</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L479" class="source_link" style="float:right">[source]</a></h4>

> <code>to_xy</code>(**`v`**)

- Rotate a 3D vector v in xy-plane.
- **Parameters**
    - v: Ponit in 3D.



<h4 id="rad_angle" class="doc_header"><code>rad_angle</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L495" class="source_link" style="float:right">[source]</a></h4>

> <code>rad_angle</code>(**`v1`**, **`v2`**)

- Returns interier angle between two vectors.
- **Parameters**
    - v1,v2 : Two vectors/points in 3D.



<h4 id="arctan_full" class="doc_header"><code>arctan_full</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L509" class="source_link" style="float:right">[source]</a></h4>

> <code>arctan_full</code>(**`perp`**, **`base`**)

- Returns full angle from x-axis counter clockwise.
- **Parameters**
    - perp: Perpendicular componet of vector including sign.
    - base: Base compoent of vector including sign.



<h4 id="get_bz" class="doc_header"><code>get_bz</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L541" class="source_link" style="float:right">[source]</a></h4>

> <code>get_bz</code>(**`poscar`**=*`None`*, **`loop`**=*`True`*, **`digits`**=*`8`*)

- Return required information to construct first Brillouin zone in form of tuple (basis, normals, vertices, faces).
- **Parameters**
    - poscar : POSCAR file or list of 3 vectors in 3D aslist[list,list,list].
    - loop   : If True, joins the last vertex of a BZ plane to starting vertex in order to complete loop.
    - digits : int, rounding off decimal places, no effect on intermediate calculations, just for pretty final results
- **Attributes**
    - basis   : get_bz().basis, recprocal lattice basis.
    - normals : get_bz().normals, all vertors that are perpendicular BZ faces/planes.
    - vertices: get_bz().vertices, all vertices of BZ, could be input into ConvexHull for constructing 3D BZ.
    - faces   : get_bz().faces, vertices arranged into faces, could be input to Poly3DCollection of matplotlib for creating BZ from faces' patches.
    - specials : get_bz().specials, Dictionary of high symmetry KPOINTS with keys as points relative to basis and values are corresponding positions in recirprocal coordinates space.



<h4 id="plot_bz" class="doc_header"><code>plot_bz</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L636" class="source_link" style="float:right">[source]</a></h4>

> <code>plot_bz</code>(**`poscar_or_bz`**=*`None`*, **`fill`**=*`True`*, **`color`**=*`'rgba(168,204,216,0.4)'`*, **`background`**=*`'rgb(255,255,255)'`*)

- Plots interactive figure showing axes,BZ surface, special points and basis, each of which could be hidden or shown.
- **Parameters**
    - pocar_or_bz: POSCAR or 3 basis vectors' list forming POSCAR. Auto picks in working directory. Output of get_bz() also works.
    - fill       : True by defult, determines whether to fill surface of BZ or not.
    - color      : color to fill surface 'rgba((168,204,216,0.4)` by default.
    - background : Plot background color, default is 'rgb(255,255,255)'.
- **Returns**
    - fig   : plotly.graph_object's Figure instance.



<h4 id="widget-toggle-button {
                    color: black !important;
                    min-width: max-content !important;
                    background-color: #c3d4d4;
                    border-radius: 5px !important;
                }
            </style>" class="doc_header"><code>widget-toggle-button {
                    color: black !important;
                    min-width: max-content !important;
                    background-color: #c3d4d4;
                    border-radius: 5px !important;
                }
            </style></code><a href="" class="source_link" style="float:right">[source]</a></h4>

str(object='') -> str
str(bytes_or_buffer[, encoding[, errors]]) -> str

Create a new string object from the given object. If encoding or
errors is specified, then the object must expose a data buffer
that will be decoded using the given encoding and error handler.
Otherwise, returns the result of object.__str__() (if defined)
or repr(object).
encoding defaults to sys.getdefaultencoding().
errors defaults to 'strict'.



<h4 id="widget-toggle-button {
                    color: 	whitesmoke !important;
                    min-width: max-content !important;
                    background-color: #3f5959;
                    border-radius: 5px !important;
                }

            </style>" class="doc_header"><code>widget-toggle-button {
                    color: 	whitesmoke !important;
                    min-width: max-content !important;
                    background-color: #3f5959;
                    border-radius: 5px !important;
                }

            </style></code><a href="" class="source_link" style="float:right">[source]</a></h4>

str(object='') -> str
str(bytes_or_buffer[, encoding[, errors]]) -> str

Create a new string object from the given object. If encoding or
errors is specified, then the object must expose a data buffer
that will be decoded using the given encoding and error handler.
Otherwise, returns the result of object.__str__() (if defined)
or repr(object).
encoding defaults to sys.getdefaultencoding().
errors defaults to 'strict'.



<h4 id="get_files_gui" class="doc_header"><code>get_files_gui</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L214" class="source_link" style="float:right">[source]</a></h4>

> <code>get_files_gui</code>(**`auto_fill`**=*`'vasprun.xml'`*, **`html_style`**=*`None`*, **`height`**=*`320`*)

- Creates a GUI interface for files/folders filtering.
- **Parmeters**
    - auto_fill  : Default is `vasprun.xml`, any file/folder.
    - html_style : None,[`dark_style`](/pivotpy/Widgets.html#dark_style) or [`light_style`](/pivotpy/Widgets.html#light_style).
    - height     : Height of Grid box.
- **Returns**
    - Tuple(GUI_gridbox,Files_Dropdown). Access second one by item itself.



<h4 id="get_input_gui" class="doc_header"><code>get_input_gui</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L311" class="source_link" style="float:right">[source]</a></h4>

> <code>get_input_gui</code>(**`rgb`**=*`True`*, **`sys_info`**=*`None`*, **`html_style`**=*`None`*, **`height`**=*`400`*)

- Creates a GUI interface for input/selection of orbitals/ions projection.
- **Parmeters**
    - rgb        : Default is `True` and generates input for `plotly(quick)_rgb_lines`, if `False` creates input for `quick(plotly)_dos(color)_lines`
    - html_style : None,[`dark_style`](/pivotpy/Widgets.html#dark_style) or [`light_style`](/pivotpy/Widgets.html#light_style).
    - height     : Height of Grid box.
- **Returns**
    - Tuple(GUI_gridbox,json_in_HTML). Access second one by item.value.



<h4 id="read_data" class="doc_header"><code>read_data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L475" class="source_link" style="float:right">[source]</a></h4>

> <code>read_data</code>(**`tabel_w`**, **`poscar`**=*`None`*, **`sys_info`**=*`None`*)





<h4 id="click_data" class="doc_header"><code>click_data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L488" class="source_link" style="float:right">[source]</a></h4>

> <code>click_data</code>(**`sel_en_w`**, **`fermi_w`**, **`tabel_w`**, **`fig`**)





<h4 id="tabulate_data" class="doc_header"><code>tabulate_data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L505" class="source_link" style="float:right">[source]</a></h4>

> <code>tabulate_data</code>(**`data_dict`**)





<h4 id="save_data" class="doc_header"><code>save_data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L543" class="source_link" style="float:right">[source]</a></h4>

> <code>save_data</code>(**`out_w1`**, **`data_dict`**)





<h4 id="color_toggle" class="doc_header"><code>color_toggle</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L549" class="source_link" style="float:right">[source]</a></h4>

> <code>color_toggle</code>(**`tog_w`**, **`fig`**, **`rd_btn`**)





<h4 id="clear_cache" class="doc_header"><code>clear_cache</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L574" class="source_link" style="float:right">[source]</a></h4>

> <code>clear_cache</code>(**`out_w1`**, **`cache_w`**, **`tabel_w`**)





<h4 id="matplotlib_code" class="doc_header"><code>matplotlib_code</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L596" class="source_link" style="float:right">[source]</a></h4>

> <code>matplotlib_code</code>(**`rd_btn`**, **`out_w1`**, **`dict_html`**)





<h4 id="generate_summary" class="doc_header"><code>generate_summary</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L611" class="source_link" style="float:right">[source]</a></h4>

> <code>generate_summary</code>(**`paths_list`**=*`None`*)





<h4 id="show_app" class="doc_header"><code>show_app</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L659" class="source_link" style="float:right">[source]</a></h4>

> <code>show_app</code>(**`height`**=*`600`*)

Displays a GUI for visulaizing and manipulating output of vasp calculations. It only has one argument `height` which sets `min-height` of the app.



```python

```
