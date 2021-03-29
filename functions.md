# Functions Reference


```python
from nbdev import show_doc
import pivotpy as pp
all_ = pp.__all__
show_doc(pp.Dict2Data)
show_doc(pp.Dict2Data.to_dict)
show_doc(pp.Dict2Data.to_json)
show_doc(pp.Dict2Data.to_pickle)
_ = [show_doc(eval('pp.{}'.format(f))) for f in all_ if f not in ['Dict2Data','_savefig','_show']]
```


<h2 id="Dict2Data" class="doc_header"><code>class</code> <code>Dict2Data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L38" class="source_link" style="float:right">[source]</a></h2>

> <code>Dict2Data</code>(**`d`**)

- Returns a Data object with dictionary keys as attributes of Data accessible by dot notation. Once an attribute is created, it can not be changed from outside.
- **Parmeters**
    - dict : Python dictionary (nested as well) containing any python data types.
- **Methods**
    - to_dict  : Converts a Data object to dictionary if it could be made a dictionary, otherwise throws relevant error.
    - to_json  : Converts to json str or save to file if `outfil` given. Accepts `indent` as parameter.
    - to_pickle: Converts to bytes str or save to file if `outfile` given.
    - to_tuple : Converts to a named tuple.
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



<h4 id="Dict2Data.to_dict" class="doc_header"><code>Dict2Data.to_dict</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L71" class="source_link" style="float:right">[source]</a></h4>

> <code>Dict2Data.to_dict</code>()

Converts a `Dict2Data` object (root or nested level) to a dictionary.
        



<h4 id="Dict2Data.to_json" class="doc_header"><code>Dict2Data.to_json</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L81" class="source_link" style="float:right">[source]</a></h4>

> <code>Dict2Data.to_json</code>(**`outfile`**=*`None`*, **`indent`**=*`1`*)

Dumps a `Dict2Data` object (root or nested level) to json.
- **Parameters**
    - outfile : Default is None and returns string. If given, writes to file.
    - indent  : Json indent. Default is 1.



<h4 id="Dict2Data.to_pickle" class="doc_header"><code>Dict2Data.to_pickle</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L89" class="source_link" style="float:right">[source]</a></h4>

> <code>Dict2Data.to_pickle</code>(**`outfile`**=*`None`*)

Dumps a `Dict2Data` object (root or nested level) to pickle.
- **Parameters**
    - outfile : Default is None and returns string. If given, writes to file.



<h4 id="dict2tuple" class="doc_header"><code>dict2tuple</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L27" class="source_link" style="float:right">[source]</a></h4>

> <code>dict2tuple</code>(**`name`**, **`d`**)

Converts a dictionary (nested as well) to namedtuple, accessible via index and dot notation as well as by unpacking.
- **Parameters**
    - name: Name of the tuple.
    - d   : Dictionary, nested works as well.



<h4 id="read_asxml" class="doc_header"><code>read_asxml</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L125" class="source_link" style="float:right">[source]</a></h4>

> <code>read_asxml</code>(**`path`**=*`None`*)

- Reads a big vasprun.xml file into memory once and then apply commands. If current folder contains `vasprun.xml` file, it automatically picks it.

- **Parameters**
    - path : Path/To/vasprun.xml

- **Returns**
    - xml_data : Xml object to use in other functions



<h4 id="xml2dict" class="doc_header"><code>xml2dict</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L162" class="source_link" style="float:right">[source]</a></h4>

> <code>xml2dict</code>(**`xmlnode_or_filepath`**)

Convert xml node or xml file content to dictionary. All output text is in string format, so further processing is required to convert into data types/split etc.
- The only paramenter `xmlnode_or_filepath` is either a path to an xml file or an `xml.etree.ElementTree.Element` object.
- Each node has `tag,text,attr,nodes` attributes. Every text element can be accessed via
`xml2dict()['nodes'][index]['nodes'][index]...` tree which makes it simple.



<h4 id="exclude_kpts" class="doc_header"><code>exclude_kpts</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L178" class="source_link" style="float:right">[source]</a></h4>

> <code>exclude_kpts</code>(**`xml_data`**=*`None`*)

- Returns number of kpoints to exclude used from IBZKPT.
- **Parameters**
    - xml_data : From `read_asxml` function
- **Returns**
    - int      : Number of kpoints to exclude.



<h4 id="get_ispin" class="doc_header"><code>get_ispin</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L199" class="source_link" style="float:right">[source]</a></h4>

> <code>get_ispin</code>(**`xml_data`**=*`None`*)

- Returns value of ISPIN.
- **Parameters**
    - xml_data : From `read_asxml` function
- **Returns**
    - int      : Value of ISPIN.



<h4 id="get_summary" class="doc_header"><code>get_summary</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L216" class="source_link" style="float:right">[source]</a></h4>

> <code>get_summary</code>(**`xml_data`**=*`None`*)

- Returns overview of system parameters.
- **Parameters**
    - xml_data : From `read_asxml` function
- **Returns**
    - Data     : pivotpy.Dict2Data with attibutes accessible via dot notation.



<h4 id="join_ksegments" class="doc_header"><code>join_ksegments</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L256" class="source_link" style="float:right">[source]</a></h4>

> <code>join_ksegments</code>(**`kpath`**, **`kseg_inds`**=*`[]`*)

Joins a broken kpath's next segment to previous. `kseg_inds` should be list of first index of next segment



<h4 id="get_kpts" class="doc_header"><code>get_kpts</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L264" class="source_link" style="float:right">[source]</a></h4>

> <code>get_kpts</code>(**`xml_data`**=*`None`*, **`skipk`**=*`0`*, **`kseg_inds`**=*`[]`*)

Returns kpoints and calculated kpath.

Parameters:

xml_data
    From `read_asxml` function.

skipk : int
    Number of initil kpoints to skip.

kseg_inds : list
    List of indices of kpoints where path is broken.

Returns:

Data : pivotpy.Dict2Data
    with attibutes `kpath` and `kpoints`.



<h4 id="get_tdos" class="doc_header"><code>get_tdos</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L299" class="source_link" style="float:right">[source]</a></h4>

> <code>get_tdos</code>(**`xml_data`**=*`None`*, **`spin_set`**=*`1`*, **`elim`**=*`[]`*)

- Returns total dos for a spin_set (default 1) and energy limit. If spin-polarized calculations, gives SpinUp and SpinDown keys as well.
- **Parameters**
    - xml_data : From `read_asxml` function
    - spin_set : int, default is 1.and
    - elim     : List [min,max] of energy, default empty.
- **Returns**
    - Data     : pivotpy.Dict2Data with attibutes E_Fermi, ISPIN,tdos.



<h4 id="get_evals" class="doc_header"><code>get_evals</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L351" class="source_link" style="float:right">[source]</a></h4>

> <code>get_evals</code>(**`xml_data`**=*`None`*, **`skipk`**=*`None`*, **`elim`**=*`[]`*)

- Returns eigenvalues as numpy array. If spin-polarized calculations, gives SpinUp and SpinDown keys as well.
- **Parameters**
    - xml_data : From `read_asxml` function
    - skipk    : Number of initil kpoints to skip.
    - elim     : List [min,max] of energy, default empty.
- **Returns**
    - Data     : pivotpy.Dict2Data with attibutes evals and related parameters.



<h4 id="get_bands_pro_set" class="doc_header"><code>get_bands_pro_set</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L405" class="source_link" style="float:right">[source]</a></h4>

> <code>get_bands_pro_set</code>(**`xml_data`**=*`None`*, **`spin_set`**=*`1`*, **`skipk`**=*`0`*, **`bands_range`**=*`None`*, **`set_path`**=*`None`*)

- Returns bands projection of a spin_set(default 1). If spin-polarized calculations, gives SpinUp and SpinDown keys as well.
- **Parameters**
    - xml_data    : From `read_asxml` function
    - skipk       : Number of initil kpoints to skip (Default 0).
    - spin_set    : Spin set to get, default is 1.
    - bands_range : If elim used in `get_evals`,that will return bands_range to use here. Note that range(0,2) will give 2 bands 0,1 but tuple (0,2) will give 3 bands 0,1,2.
    - set_path     : path/to/_set[1,2,3,4].txt, works if `split_vasprun` is used before.
- **Returns**
    - Data     : pivotpy.Dict2Data with attibutes of bands projections and related parameters.



<h4 id="get_dos_pro_set" class="doc_header"><code>get_dos_pro_set</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L500" class="source_link" style="float:right">[source]</a></h4>

> <code>get_dos_pro_set</code>(**`xml_data`**=*`None`*, **`spin_set`**=*`1`*, **`dos_range`**=*`None`*)

- Returns dos projection of a spin_set(default 1) as numpy array. If spin-polarized calculations, gives SpinUp and SpinDown keys as well.
- **Parameters**
    - xml_data    : From `read_asxml` function
    - spin_set    : Spin set to get, default 1.
    - dos_range   : If elim used in `get_tdos`,that will return dos_range to use here..
- **Returns**
    - Data     : pivotpy.Dict2Data with attibutes of dos projections and related parameters.



<h4 id="get_structure" class="doc_header"><code>get_structure</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L541" class="source_link" style="float:right">[source]</a></h4>

> <code>get_structure</code>(**`xml_data`**=*`None`*)

- Returns structure's volume,basis,positions and rec-basis.
- **Parameters**
    - xml_data : From `read_asxml` function.
- **Returns**
    - Data     : pivotpy.Dict2Data with attibutes volume,basis,positions rec_basis and labels.



<h4 id="export_vasprun" class="doc_header"><code>export_vasprun</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L582" class="source_link" style="float:right">[source]</a></h4>

> <code>export_vasprun</code>(**`path`**=*`None`*, **`skipk`**=*`None`*, **`elim`**=*`[]`*, **`kseg_inds`**=*`[]`*, **`shift_kpath`**=*`0`*, **`try_pwsh`**=*`True`*)

- Returns a full dictionary of all objects from `vasprun.xml` file. It first try to load the data exported by powershell's `Export-VR(Vasprun)`, which is very fast for large files. It is recommended to export large files in powershell first.
- **Parameters**
    - path       : Path to `vasprun.xml` file. Default is `'./vasprun.xml'`.
    - skipk      : Default is None. Automatically detects kpoints to skip.
    - elim       : List [min,max] of energy interval. Default is [], covers all bands.
    - kseg_inds : List of indices of kpoints where path is broken.
    - shift_kpath: Default 0. Can be used to merge multiple calculations on single axes side by side.
    - try_pwsh   : Default is True and tries to load data exported by `Vasp2Visual` in Powershell.
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



<h4 id="load_export" class="doc_header"><code>load_export</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L674" class="source_link" style="float:right">[source]</a></h4>

> <code>load_export</code>(**`path`**=*`'./vasprun.xml'`*, **`kseg_inds`**=*`[]`*, **`shift_kpath`**=*`0`*, **`path_to_ps`**=*`'pwsh'`*, **`skipk`**=*`None`*, **`max_filled`**=*`10`*, **`max_empty`**=*`10`*, **`keep_files`**=*`True`*)

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



<h4 id="dump_dict" class="doc_header"><code>dump_dict</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L827" class="source_link" style="float:right">[source]</a></h4>

> <code>dump_dict</code>(**`dict_data`**=*`None`*, **`dump_to`**=*`'pickle'`*, **`outfile`**=*`None`*, **`indent`**=*`1`*)

- Dump an `export_vasprun` or `load_export`'s `Data` object or any dictionary to json or pickle string/file. It convert `Dict2Data` to dictionary before serializing to json/pickle, so json/pickle.loads() of converted Data would be a simple dictionary, pass that to `Dict2Data` to again make accessible via dot notation.
- **Parameters**
    - dict_data : Any dictionary/Dict2Data object containg numpy arrays, including `export_vasprun` or `load_export` output.
    - dump_to  : Defualt is `pickle` or `json`.
    - outfile  : Defualt is None and return string. File name does not require extension.
    - indent   : Defualt is 1. Only works for json.



<h4 id="load_from_dump" class="doc_header"><code>load_from_dump</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L857" class="source_link" style="float:right">[source]</a></h4>

> <code>load_from_dump</code>(**`file_or_str`**, **`keep_as_dict`**=*`False`*)

- Loads a json/pickle dumped file or string by auto detecting it.
- **Parameters**
    - file_or_str : Filename of pickl/json or their string.
    - keep_as_dict: Defualt is False and return `Data` object. If True, returns dictionary.



<h4 id="islice2array" class="doc_header"><code>islice2array</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L887" class="source_link" style="float:right">[source]</a></h4>

> <code>islice2array</code>(**`path_or_islice`**, **`dtype`**=*`float`*, **`delimiter`**=*`'\\s+'`*, **`include`**=*`None`*, **`exclude`**=*`'#'`*, **`raw`**=*`False`*, **`fix_format`**=*`True`*, **`start`**=*`0`*, **`nlines`**=*`None`*, **`count`**=*`-1`*, **`cols`**=*`None`*, **`new_shape`**=*`None`*)

- Reads a sliced array from txt,csv type files and return to array. Also manages if columns lengths are not equal and return 1D array. It is faster than loading  whole file into memory. This single function could be used to parse EIGENVAL, PROCAR, DOCAR and similar files with just a combination of `exclude, include,start,stop,step` arguments.
- **Parameters**
    - path_or_islice: Path/to/file or `itertools.islice(file_object)`. islice is interesting when you want to read different slices of an opened file and do not want to open it again and again. For reference on how to use it just execute `pivotpy.export_potential??` in a notebook cell or ipython terminal to see how islice is used extensively.
    - dtype: float by default. Data type of output array, it is must have argument.
    - start,nlines: The indices of lines to start reading from and number of lines after start respectively. Only work if `path_or_islice` is a file path. both could be None or int, while start could be a list to read slices from file provided that nlines is int. The spacing between adjacent indices in start should be equal to or greater than nlines as pointer in file do not go back on its own.  These parameters are in output of `slice_data`
    > Note: `start` should count comments if `exclude` is None. You can use `slice_data` function to get a dictionary of `start,nlines, count, cols, new_shape` and unpack in argument instead of thinking too much.
    - count: `np.size(output_array) = nrows x ncols`, if it is known before execution, performance is increased. This parameter is in output of `slice_data`.
    - delimiter:  Default is `\s+`. Could be any kind of delimiter valid in numpy and in the file.
    - cols: List of indices of columns to pick. Useful when reading a file like PROCAR which e.g. has text and numbers inline. This parameter is in output of `slice_data`.
    - include: Default is None and includes everything. String of patterns separated by | to keep, could be a regular expression.
    - exclude: Default is '#' to remove comments. String of patterns separated by | to drop,could be a regular expression.
    - raw    : Default is False, if True, returns list of raw strings. Useful to select `cols`.
    - fix_format: Default is True, it sepearates numbers with poor formatting like 1.000-2.000 to 1.000 2.000 which is useful in PROCAR. Keep it False if want to read string literally.
    - new_shape : Tuple of shape Default is None. Will try to reshape in this shape, if fails fallbacks to 2D or 1D. This parameter is in output of `slice_data`.
- **Examples**
    > `islice2array('path/to/PROCAR',start=3,include='k-point',cols=[3,4,5])[:2]`
    > array([[ 0.125,  0.125,  0.125],
    >        [ 0.375,  0.125,  0.125]])
    > `islice2array('path/to/EIGENVAL',start=7,exclude='E',cols=[1,2])[:2]`
    > array([[-11.476913,   1.      ],
    >        [  0.283532,   1.      ]])
> Note: Slicing a dimension to 100% of its data is faster than let say 80% for inner dimensions, so if you have to slice more than 50% of an inner dimension, then just load full data and slice after it.



<h4 id="slice_data" class="doc_header"><code>slice_data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L973" class="source_link" style="float:right">[source]</a></h4>

> <code>slice_data</code>(**`dim_inds`**, **`old_shape`**)

- Returns a dictionary that can be unpacked in arguments of isclice2array function. This function works only for regular txt/csv/tsv data files which have rectangular data written.
- **Parameters**
    - dim_inds : List of indices array or range to pick from each dimension. Inner dimensions are more towards right. Last itmes in dim_inds is considered to be columns. If you want to include all values in a dimension, you can put -1 in that dimension. Note that negative indexing does not work in file readig, -1 is s special case to fetch all items.
    - old_shape: Shape of data set including the columns length in right most place.
- **Example**
    - You have data as 3D arry where third dimension is along column.
    > 0 0
    > 0 2
    > 1 0
    > 1 2
    - To pick [[0,2], [1,2]], you need to give
    > slice_data(dim_inds = [[0,1],[1],-1], old_shape=(2,2,2))
    > {'start': array([1, 3]), 'nlines': 1, 'count': 2}
    - Unpack above dictionary in `islice2array` and you will get output array.
- Note that dimensions are packed from right to left, like 0,2 is repeating in 2nd column.



<h4 id="split_vasprun" class="doc_header"><code>split_vasprun</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L1024" class="source_link" style="float:right">[source]</a></h4>

> <code>split_vasprun</code>(**`path`**=*`None`*)

- Splits a given vasprun.xml file into a smaller _vasprun.xml file plus _set[1,2,3,4].txt files which contain projected data for each spin set.
- **Parameters**
    - path: path/to/vasprun.xml file.
- **Output**
    - _vasprun.xml file with projected data.
    - _set1.txt for projected data of colinear calculation.
    - _set1.txt for spin up data and _set2.txt for spin-polarized case.
    - _set[1,2,3,4].txt for each spin set of non-colinear calculations.



<h4 id="get_file_size" class="doc_header"><code>get_file_size</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L31" class="source_link" style="float:right">[source]</a></h4>

> <code>get_file_size</code>(**`path`**)





<h4 id="interpolate_data" class="doc_header"><code>interpolate_data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L43" class="source_link" style="float:right">[source]</a></h4>

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



<h4 id="ps2py" class="doc_header"><code>ps2py</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L72" class="source_link" style="float:right">[source]</a></h4>

> <code>ps2py</code>(**`ps_command`**=*`'Get-ChildItem'`*, **`exec_type`**=*`'-Command'`*, **`path_to_ps`**=*`'powershell.exe'`*)

- Captures powershell output in python.
- **Parameters**
    - ps_command: enclose ps_command in ' ' or " ".
    - exec_type : type of execution, default '-Command', could be '-File'.
    - path_to_ps: path to powerhell.exe if not added to PATH variables.



<h4 id="ps2std" class="doc_header"><code>ps2std</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L105" class="source_link" style="float:right">[source]</a></h4>

> <code>ps2std</code>(**`ps_command`**=*`'Get-ChildItem'`*, **`exec_type`**=*`'-Command'`*, **`path_to_ps`**=*`'powershell.exe'`*)

- Prints powershell output in python std.
- **Parameters**
    - ps_command: enclose ps_command in ' ' or " ".
    - exec_type: type of execution, default '-Command', could be '-File'.
    - path_to_ps: path to powerhell.exe if not added to PATH variables.



<h4 id="get_child_items" class="doc_header"><code>get_child_items</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L119" class="source_link" style="float:right">[source]</a></h4>

> <code>get_child_items</code>(**`path`**=*`'C:\\Users\\mass_\\AppData\\Local\\Temp'`*, **`depth`**=*`None`*, **`recursive`**=*`True`*, **`include`**=*`None`*, **`exclude`**=*`None`*, **`filesOnly`**=*`False`*, **`dirsOnly`**=*`False`*)

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



<h2 id="color" class="doc_header"><code>class</code> <code>color</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L159" class="source_link" style="float:right">[source]</a></h2>

> <code>color</code>()





<h2 id="EncodeFromNumpy" class="doc_header"><code>class</code> <code>EncodeFromNumpy</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L181" class="source_link" style="float:right">[source]</a></h2>

> <code>EncodeFromNumpy</code>(**`skipkeys`**=*`False`*, **`ensure_ascii`**=*`True`*, **`check_circular`**=*`True`*, **`allow_nan`**=*`True`*, **`sort_keys`**=*`False`*, **`indent`**=*`None`*, **`separators`**=*`None`*, **`default`**=*`None`*) :: `JSONEncoder`

- Serializes python/Numpy objects via customizing json encoder.
- **Usage**
    - `json.dumps(python_dict, cls=EncodeFromNumpy)` to get json string.
    - `json.dump(*args, cls=EncodeFromNumpy)` to create a file.json.



<h2 id="DecodeToNumpy" class="doc_header"><code>class</code> <code>DecodeToNumpy</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L208" class="source_link" style="float:right">[source]</a></h2>

> <code>DecodeToNumpy</code>(**\*`args`**, **\*\*`kwargs`**) :: `JSONDecoder`

- Deserilizes JSON object to Python/Numpy's objects.
- **Usage**
    - `json.loads(json_string,cls=DecodeToNumpy)` from string, use `json.load()` for file.



<h2 id="Vasprun" class="doc_header"><code>class</code> <code>Vasprun</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L238" class="source_link" style="float:right">[source]</a></h2>

> <code>Vasprun</code>(**`path`**=*`None`*, **`skipk`**=*`None`*, **`elim`**=*`[]`*, **`kseg_inds`**=*`[]`*, **`shift_kpath`**=*`0`*)

- All plotting functions that depend on `export_vasprun` are joined under this class and renamed.
- **Parameters**
    - path       : str: path/to/vasprun.xml. Auto picks in CWD.
    - skipk      : int: Skip initial kpoints
    - elim       : list: Energy range e.g. [-5,5]
    - kseg_inds : list: Join broken path at given indices. Could be obtained from `SEG-INDS` if used `trace_kpath`.
    - shift_kpath: float: Shift in kpath values for side by side plotting.
- **Attributes**
    - data : Return of `export_vasprun` which is auto-picked in plotting methods under this class.
- **Methods**
    - sbands    : Shortcut for `quick_bplot`.
    - sdos      : Shortcut for `quick_dos_lines`.
    - srgb      : Shortcut for `quick_rgb_lines`.
    - scolor    : Shortcut for `quick_color_lines`.
    - idos      : Shortcut for `plotly_dos_lines`.
    - irgb      : Shortcut for `plotly_rgb_lines`.
    - Each of above mathods have an attribute `kwargs` which can be accessed, modified and put back as argumnets.
- **Example**
    > vasp   = Vasprun(path='./vasprun.xml')
    > kwargs = vasp.sbands.kwargs
    > Modify kwargs dictionary as you want for input parameters and unpack back in function.
    > vasp.sbands(**kwargs)



<h4 id="nav_links" class="doc_header"><code>nav_links</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L293" class="source_link" style="float:right">[source]</a></h4>

> <code>nav_links</code>(**`current_index`**=*`0`*, **`doc_url`**=*`'https://massgh.github.io/pivotpy/'`*, **`items`**=*`['Index', 'XmlElementTree', 'StaticPlots', 'InteractivePlots', 'Utilities', 'StructureIO', 'Widgets']`*, **`horizontal`**=*`False`*, **`out_string`**=*`False`*)





<h4 id="export_outcar" class="doc_header"><code>export_outcar</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L321" class="source_link" style="float:right">[source]</a></h4>

> <code>export_outcar</code>(**`path`**=*`None`*)

- Read potential at ionic sites from OUTCAR.



<h4 id="export_potential" class="doc_header"><code>export_potential</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L371" class="source_link" style="float:right">[source]</a></h4>

> <code>export_potential</code>(**`locpot`**=*`None`*, **`e`**=*`True`*, **`m`**=*`False`*)

- Returns Data from LOCPOT and similar structure files like CHG. Loads only single set out of 2/4 magnetization data to avoid performance/memory cost while can load electrostatic and one set of magnetization together.
- **Parameters**
    - locpot: path/to/LOCPOT or similar stuructured file like CHG. LOCPOT is auto picked in CWD.
    - e     : Electric potential/charge density. Default is True.
    - m     : Magnetization density m. Default is False. If True, picks `m` for spin polarized case, and `m_x` for non-colinear case. Additionally it can take 'x','y' and 'z' in case of non-colinear calculations.
- **Exceptions**
    - Would raise index error if magnetization density set is not present in LOCPOT/CHG in case `m` is not False.



<h2 id="LOCPOT_CHG" class="doc_header"><code>class</code> <code>LOCPOT_CHG</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L451" class="source_link" style="float:right">[source]</a></h2>

> <code>LOCPOT_CHG</code>(**`path`**=*`None`*, **`e`**=*`True`*, **`m`**=*`False`*)

- Returns Data from LOCPOT and similar structure files like CHG. Loads only single set out of 2/4 magnetization data to avoid performance/memory cost while can load electrostatic and one set of magnetization together.
- **Parameters**
    - path: path/to/LOCPOT or similar stuructured file like CHG. LOCPOT is auto picked in CWD.
    - e   : Electric potential/charge density. Default is True.
    - m   : Magnetization density m. Default is False. If True, picks `m` for spin polarized case, and `m_x` for non-colinear case. Additionally it can take 'x','y' and 'z' in case of non-colinear calculations.
- **Exceptions**
    - Would raise index error if magnetization density set is not present in LOCPOT/CHG in case `m` is not False.



<h4 id="transform_color" class="doc_header"><code>transform_color</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L558" class="source_link" style="float:right">[source]</a></h4>

> <code>transform_color</code>(**`arr`**, **`s`**=*`1`*, **`c`**=*`1`*, **`b`**=*`0`*, **`mixing_matrix`**=*`None`*)

- Color transformation such as brightness, contrast, saturation and mixing of an input color array. `c = -1` would invert color,keeping everything else same.
- **Parameters**
    - arr: input array, a single RGB/RGBA color or an array with inner most dimension equal to 3 or 4. e.g. [[[0,1,0,1],[0,0,1,1]]].
    - c  : contrast, default is 1. Can be a float in [-1,1].
    - s  : saturation, default is 1. Can be a float in [-1,1]. If s = 0, you get a gray scale image.
    - b  : brightness, default is 0. Can be a float in [-1,1] or list of three brightnesses for RGB components.
    - mixing_matrix: A 3x3 matrix to mix RGB values, such as `pp.color_matrix`.

[Recoloring](https://docs.microsoft.com/en-us/windows/win32/gdiplus/-gdiplus-recoloring-use?redirectedfrom=MSDN)
[Rainmeter](https://docs.rainmeter.net/tips/colormatrix-guide/)



<h4 id="modify_axes" class="doc_header"><code>modify_axes</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L46" class="source_link" style="float:right">[source]</a></h4>

> <code>modify_axes</code>(**`ax`**=*`None`*, **`xticks`**=*`[]`*, **`ktick_vals`**=*`[]`*, **`xlim`**=*`[]`*, **`yticks`**=*`[]`*, **`yt_labels`**=*`[]`*, **`ylim`**=*`[]`*, **`xlabel`**=*`None`*, **`ylabel`**=*`None`*, **`vlines`**=*`True`*, **`zeroline`**=*`True`*)

- Returns None, applies given settings on axes. Prefered to use before other plotting.
- **Parameters**
    - ax  : Matplotlib axes object.
    - (x,y)ticks : List of positions on (x,y axes).
    - (xt,yt)_labels : List of labels on (x,y) ticks points.
    - (x,y)lim : [min, max] of (x,y) axes.
    - (x,y)label : axes labels.
    - vlines : If true, drawn when `ylim` is not empty.
    - zeroline : If True, drawn when `xlim` is not empty.



<h4 id="init_figure" class="doc_header"><code>init_figure</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L85" class="source_link" style="float:right">[source]</a></h4>

> <code>init_figure</code>(**`figsize`**=*`(3.4, 2.6)`*, **`nrows`**=*`1`*, **`ncols`**=*`1`*, **`widths`**=*`[]`*, **`heights`**=*`[]`*, **`axes_off`**=*`[]`*, **`axes_3d`**=*`[]`*, **`sharex`**=*`False`*, **`sharey`**=*`False`*, **`azim`**=*`45`*, **`elev`**=*`15`*, **`ortho3d`**=*`True`*, **\*\*`subplots_adjust_kwargs`**)

- Returns flatten axes of initialized figure, based on plt.subplots(). If you want to access parent figure, use ax.get_figure() or current figure as plt.gcf().
- **Parameters**
    - figsize   : Tuple (width, height). Default is (3.4,2.6).
    - nrows     : Default 1.
    - ncols     : Default 1.
    - widths    : List with len(widths)==nrows, to set width ratios of subplots.
    - heights   : List with len(heights)==ncols, to set height ratios of subplots.
    - share(x,y): Share axes between plots, this removes shared ticks automatically.
    - axes_off  : Turn off axes visibility, If `nrows = ncols = 1, set True/False`, If anyone of `nrows or ncols > 1`, provide list of axes indices to turn off. If both `nrows and ncols > 1`, provide list of tuples (x_index,y_index) of axes.
    - axes_3d   : Change axes to 3D. If `nrows = ncols = 1, set True/False`, If anyone of `nrows or ncols > 1`, provide list of axes indices to turn off. If both `nrows and ncols > 1`, provide list of tuples (x_index,y_index) of axes.
    azim,elev   : Matplotlib's 3D angles, defualt are 45,15.
    ortho3d     : Only works for 3D axes. If True, x,y,z are orthogonal, otherwise perspective.
    - **subplots_adjust_kwargs : These are same as `plt.subplots_adjust()`'s arguements.



<h4 id="plot_bands" class="doc_header"><code>plot_bands</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L159" class="source_link" style="float:right">[source]</a></h4>

> <code>plot_bands</code>(**`ax`**=*`None`*, **`kpath`**=*`None`*, **`bands`**=*`None`*, **`showlegend`**=*`False`*, **`E_Fermi`**=*`None`*, **`color1`**=*`(0, 0, 0.8)`*, **`style1`**=*`'solid'`*, **`lw1`**=*`0.7`*, **`color2`**=*`(0.8, 0, 0)`*, **`style2`**=*`'dashed'`*, **`lw2`**=*`0.7`*)

- Returns axes object and plot on which all matplotlib allowed actions could be performed.
- **Parameters**
    - ax         : Matplotlib axes object, if not given, one is created.
    - kpath      : 1D array from `get_kpts`().kpath or `export_vasprun`().kpath.
    - bands      : Dictionary Object from `get_evals` or `export_vasprun`().bands.
    - showlegend : Boolean, default is False, if true, gives legend for spin-polarized calculations.
    - E_Fermi    : If not given, automatically picked from bands object.
    - **kwargs   : lines color,width and style to distinguish spin Up and Down.
- **Returns**
    - ax : matplotlib axes object with plotted bands.



<h4 id="quick_bplot" class="doc_header"><code>quick_bplot</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L209" class="source_link" style="float:right">[source]</a></h4>

> <code>quick_bplot</code>(**`path_evr`**=*`None`*, **`ax`**=*`None`*, **`skipk`**=*`None`*, **`kseg_inds`**=*`[]`*, **`elim`**=*`[]`*, **`ktick_inds`**=*`[]`*, **`ktick_vals`**=*`[]`*, **`E_Fermi`**=*`None`*, **`figsize`**=*`(3.4, 2.6)`*, **`txt`**=*`None`*, **`xytxt`**=*`[0.2, 0.9]`*, **`ctxt`**=*`'black'`*)

- Returns axes object and plot on which all matplotlib allowed actions could be performed.
- **Parameters**
    - path_evr   : path/to/vasprun.xml or output of `export_vasprun`. Auto picks in CWD.
    - ax         : Matplotlib axes object, if not given, one is created.
    - skipk      : Number of kpoints to skip, default will be from IBZKPT.
    - kseg_inds : Points where kpath is broken.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from `export_vasprun`.
    - ktick_inds : High symmetry kpoints indices.abs
    - ktick_vals  : High Symmetry kpoints labels.
    - **kwargs   : figsize=(3.4,2.6). Text,its position and color.
- **Returns**
    - ax : matplotlib axes object with plotted bands.



<h4 id="add_text" class="doc_header"><code>add_text</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L263" class="source_link" style="float:right">[source]</a></h4>

> <code>add_text</code>(**`ax`**=*`None`*, **`xs`**=*`0.25`*, **`ys`**=*`0.9`*, **`txts`**=*`'[List]'`*, **`colors`**=*`'r'`*, **`transform`**=*`True`*, **\*\*`kwargs`**)

- Adds text entries on axes, given single string or list.
- **Parameters**
    - xs    : List of x coordinates relative to axes or single coordinate.
    - ys    : List of y coordinates relative to axes or single coordinate.
    - txts  : List of strings or one string.
    - colors: List of x colors of txts or one color.
    - transform: Dafault is True and positions are relative to axes, If False, positions are in data coordinates.
    - kwargs: plt.text key words arguments except bbox,ha,va, transform.



<h4 id="add_legend" class="doc_header"><code>add_legend</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L292" class="source_link" style="float:right">[source]</a></h4>

> <code>add_legend</code>(**`ax`**=*`None`*, **`colors`**=*`[]`*, **`labels`**=*`[]`*, **`styles`**=*`'solid'`*, **`widths`**=*`0.7`*, **`anchor`**=*`(0, 1)`*, **`ncol`**=*`3`*, **`loc`**=*`'lower left'`*, **`fontsize`**=*`'small'`*, **`frameon`**=*`False`*, **\*\*`legend_kwargs`**)

- Adds custom legeneds on a given axes,returns None.
- **Parameters**
    - ax       : Matplotlib axes.
    - colors   : List of colors.
    - labels   : List of labels.
    - styles   : str or list of line styles.
    - widths   : str or list of line widths.
    - **kwargs : Matplotlib's legend arguments.



<h4 id="add_colorbar" class="doc_header"><code>add_colorbar</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L325" class="source_link" style="float:right">[source]</a></h4>

> <code>add_colorbar</code>(**`ax`**=*`None`*, **`cmap_or_clist`**=*`None`*, **`N`**=*`256`*, **`ticks`**=*`[0.16666666666666666, 0.5, 0.8333333333333334]`*, **`ticklabels`**=*`['r', 'g', 'b']`*, **`vertical`**=*`False`*, **`fontsize`**=*`8`*)

- Plots colorbar on a given axes. This axes should be only for colorbar. Returns None or throws ValueError for given colors.
- **Parameters**
    - ax         : Matplotlib axes object.
    - cmap_or_clist: List/array of colors in or colormap's name. If None(default), first tries to get latest `quick_rgb_lines` colormap and if no success, then `RGB_f` colorbar is added. If nothing works, matplotlib's default colormap is plotted.
    - N          : int, number of color points Default 256.
    - ticks      : List of tick values to show on colorbar in interval [0,1].
    - ticklabels : List of labels for ticks.
    - vertical   : Boolean, default is Fasle.
    - fontsize   : Default 8. Adjustable according to plot space.
- **Note**: Use `'RGB_f'` to map colors (after) plotted in `quick_rgb_lines` and use `'RGB_m'` to plot DOS in same colors.



<h4 id="color_wheel" class="doc_header"><code>color_wheel</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L388" class="source_link" style="float:right">[source]</a></h4>

> <code>color_wheel</code>(**`ax`**=*`None`*, **`xy`**=*`(1, 1)`*, **`scale`**=*`0.12`*, **`rlim`**=*`(0.2, 1)`*, **`N`**=*`256`*, **`colormap`**=*`None`*, **`ticks`**=*`[0.16666666666666666, 0.5, 0.8333333333333334]`*, **`labels`**=*`['s', 'p', 'd']`*, **`showlegend`**=*`True`*)

- Returns cax i.e. color wheel axes.
- **Parameters**
    - ax        : Axes on which color wheel will be drawn. Auto created if not given.
    - xy        : (x,y) of center of wheel.
    - scale     : Scale of the cax internally generated by color wheel.
    - rlim      : Values in [0,1] interval, to make donut like shape.
    - N         : Number of segments in color wheel.
    - colormap : Matplotlib's color map name. Auto picks `RGB_f` and if fails, fallbacks to `viridis`.
    - ticks     : Ticks in fractions in interval [0,1].
    - labels    : Ticks labels.
    - showlegend: True or False.



<h4 id="get_pros_data" class="doc_header"><code>get_pros_data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L441" class="source_link" style="float:right">[source]</a></h4>

> <code>get_pros_data</code>(**`kpath`**=*`None`*, **`evals_set`**=*`None`*, **`pros_set`**=*`None`*, **`elements`**=*`[[0]]`*, **`orbs`**=*`[[0]]`*, **`interpolate`**=*`False`*, **`scale_data`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*)

- Returns selected elements/orbitals data.
- **Parameters**
    - kapath   : `export_vasprun`().kpath or `get_kpts`().kpath.
    - evals_set: `export_vasprun`().bands.evals or `get_evals`().evals. If calculations are spin-polarized, it will be `...evals.SpinUp/SpinDown` for both.
    - pros_set : `export_vasprun().pro_bands.pros` or `get_bands_pro_set`().pros. If calculations are spin-polarized, it will be `...pros.SpinUp/SpinDown` for both.
    - elements : Lists of list of ions to project on, could be `range(start,stop,step)` as well, remember that `stop` is not included in python. so `range(0,2)` will generate 0 and 1 indices.
    - orbs     : List of lists of orbitals indices.
    - scale_data : If True, normalizes projection data to 1.
    - interpolate: Deafult is false, if True, it will add n points between nearest kpoints.
    - n        : int, default is 5. Adds n points between nearest kpoints.
    - k        : int, order of interpolation, defualt is 3. `n > k` should be hold.
- **Returns**
    - A dictionary with keys 'kpath', 'evals', 'colors' that can be unpacked in next function's arguemnet.
        - 'kpath' : Given or interpolated kpath of shape (NKPTS,).
        - 'evals' : Given or interpolated evals_set of shape (NKPTS,NBANDS).
        - 'pros': An array of shape (NKPTS,NBANDS,n) where n is length of input elements list. If scale_data = True, normalizes this array to 1.



<h4 id="make_line_collection" class="doc_header"><code>make_line_collection</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L499" class="source_link" style="float:right">[source]</a></h4>

> <code>make_line_collection</code>(**`max_width`**=*`None`*, **`colors_list`**=*`None`*, **`rgb`**=*`False`*, **`uni_width`**=*`False`*, **`scale_color`**=*`False`*, **\*\*`pros_data`**)

- Returns a tuple of line collections. If rgb = True (at len(orbs) = 3 in `get_pros_data`), returns a tuple of two entries, multicolored line collection and RGB maximum values, otherwise return tuple of single colored multiple lines.
- **Parametrs**
    - **pros_data: Output dictionary from `get_pros_data` containing kpath, evals and colors arrays.
    - max_width  : Default is None and max linwidth = 2.5. Max inewidth is scaled to max_width if an int of float is given.
    - colors_list: List of colors for multiple lines, length equal to 3rd axis length of colors.
    - rgb        : Default is False. If True and np.shape(colors)[-1] == 3, RGB line collection is returned in a tuple of length 1. Tuple is just to support iteration.
    - uni_width  : Default is False, If True, makes linewidth uniform at width = max_width/2.
    - scale_color: If True, normalizes each point's color value, as (0,0,0.5) --> (0,0,1). If False, clips colors in range [0,1] but does not effect linewidth.



<h4 id="plot_collection" class="doc_header"><code>plot_collection</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L576" class="source_link" style="float:right">[source]</a></h4>

> <code>plot_collection</code>(**`gpd_args`**, **`mlc_args`**, **`axes`**=*`None`*)

- Plots line collection from the output of get_pros_data and make_line_collection on axes.
- **Parameters**
    - gpd_args: Dictionary of arguments from function `get_pros_data`. Do not unpack it.
    - mlc_args: Dictionary of arguments from function `make_line_collection`. Do not unpack it.
    - axes    : A single or list of matplotlib's Axes. len(list) should be equal to len(orbs) given in `get_pros_data`. If axes = None, auto generated.
- **Returns**
    - axes:  axes to return are spacially useful when axes = None, you can perform other actions on those axes. It will be a list of axes and all items could be same, depending on whether one are many axes were given/generated.



<h4 id="quick_rgb_lines" class="doc_header"><code>quick_rgb_lines</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L598" class="source_link" style="float:right">[source]</a></h4>

> <code>quick_rgb_lines</code>(**`path_evr`**=*`None`*, **`ax`**=*`None`*, **`skipk`**=*`None`*, **`kseg_inds`**=*`[]`*, **`elim`**=*`[]`*, **`elements`**=*`[[0], [], []]`*, **`orbs`**=*`[[0], [], []]`*, **`labels`**=*`['Elem0-s', '', '']`*, **`max_width`**=*`None`*, **`ktick_inds`**=*`[0, -1]`*, **`ktick_vals`**=*`['$\\Gamma$', 'M']`*, **`E_Fermi`**=*`None`*, **`figsize`**=*`(3.4, 2.6)`*, **`txt`**=*`None`*, **`xytxt`**=*`[0.2, 0.9]`*, **`ctxt`**=*`'black'`*, **`uni_width`**=*`False`*, **`interpolate`**=*`False`*, **`spin`**=*`'both'`*, **`n`**=*`5`*, **`k`**=*`3`*, **`scale_color`**=*`True`*, **`scale_data`**=*`True`*, **`colorbar`**=*`True`*, **`color_matrix`**=*`None`*)

- Returns axes object and plot on which all matplotlib allowed actions could be performed. In this function,orbs,labels,elements all have list of length 3. Inside list, sublists or strings could be any length but should be there even if empty.
- **Parameters**
    - path_evr   : path/to/vasprun.xml or output of `export_vasprun`. Auto picks in CWD.
    - ax         : Matplotlib axes object, if not given, one is created.
    - skipk      : Number of kpoints to skip, default will be from IBZKPT.
    - kseg_inds : Points where kpath is broken.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from `export_vasprun`.
    - ktick_inds : High symmetry kpoints indices.abs
    - ktick_vals  : High Symmetry kpoints labels.
    - elements   : List [[0],[],[]] by default and plots s orbital of first ion..
    - orbs       : List [[r],[g],[b]] of indices of orbitals, could be empty, but shape should be same.
    - labels     : List [str,str,str] of projection labels. empty string should exist to maintain shape. Auto adds `↑`,`↓` for ISPIN=2. If a label is empty i.e. '', it will not show up in colorbar ticks or legend.
    - max_width  : Width to scale whole projections. if `uni_width=True, width=max_width/2`. Default is None and linewidth at any point = 2.5*sum(ions+orbitals projection of all three input at that point). Linewidth is scaled to max_width if an int or float is given.
    - figsize    : Tuple (width,height) in inches. Default (3.4.2.6) is article column's width.
    - txt        : Text on figure, if None, SYSTEM's name is printed.
    - xytxt      : [x_coord,y_coord] of text relative to axes.
    - ctxt       : color of text.
    - uni_width  : If True, width of bands kept uniform.
    - spin       : Plot spin-polarized for spin {'up','down','both'}. Default is both.
    - interpolate: Default is False, if True, bands are interpolated.
    - n          : int, number of points, default is 5.
    - k          : int, order of interpolation 0,1,2,3. Defualt 3. `n > k` should be hold.
    - scale_color: Boolean. Default True, colors are scaled to 1 at each point. If False, clips colors in range [0,1] but does not effect linewidth.
    - scale_data : Default is True and normalizes projection data to 1. Has no visual effect if scale_color = True too.
    - colorbar   : Default is True. Displays a vertical RGB colorbar.
    - color_matrix: Only works if `scale_color==True`. 3x3 or 3x4 numpy array or list to transform from RGB to another space,provided that sum(color_matrix[i,:3]) <= 1. 4th column, if given can be used to control the saturation,contrast and brightness as s,c,b = color_matrix[:,3] For simply changing the color intensity use np.diag([r,g,b]) with r,g,b interval in [0,1]. Try `pivotpy.color_matrix` as suggested color matrix and modify, which at s=0 returns gray scale.!
- **Returns**
    - ax : matplotlib axes object with plotted projected bands.
    - Registers as colormap `RGB_m` to use in DOS to plot in same colors and `RGB_f` to display bands colorbar on another axes.
> Note: Two figures made by this function could be comapred quantitatively only if `scale_data=False, max_width=None, scale_color=False` as these parameters act internally on data.



<h4 id="quick_color_lines" class="doc_header"><code>quick_color_lines</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L811" class="source_link" style="float:right">[source]</a></h4>

> <code>quick_color_lines</code>(**`path_evr`**=*`None`*, **`axes`**=*`None`*, **`skipk`**=*`None`*, **`kseg_inds`**=*`[]`*, **`elim`**=*`[]`*, **`elements`**=*`[[0]]`*, **`orbs`**=*`[[0]]`*, **`labels`**=*`['s']`*, **`colormap`**=*`'gist_rainbow'`*, **`scale_data`**=*`False`*, **`max_width`**=*`None`*, **`spin`**=*`'both'`*, **`ktick_inds`**=*`[0, -1]`*, **`ktick_vals`**=*`['$\\Gamma$', 'M']`*, **`E_Fermi`**=*`None`*, **`showlegend`**=*`True`*, **`figsize`**=*`(3.4, 2.6)`*, **`txt`**=*`None`*, **`xytxt`**=*`[0.2, 0.85]`*, **`ctxt`**=*`'black'`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*, **`legend_kwargs`**=*`{'ncol': 4, 'anchor': (0, 1.05), 'handletextpad': 0.5, 'handlelength': 1, 'fontsize': 'small', 'frameon': False}`*, **\*\*`subplots_adjust_kwargs`**)

- Returns axes object and plot on which all matplotlib allowed actions could be performed. If given, elements, orbs, and labels must have same length. If not given, zeroth ion is plotted with s-orbital.
- **Parameters**
    - path_evr   : Path/to/vasprun.xml or output of `export_vasprun`. Auto picks in CWD.
    - axes       : Matplotlib axes object with one or many axes, if not given, auto created.
    - skipk      : Number of kpoints to skip, default will be from IBZKPT.
    - kseg_inds : Points where kpath is broken.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from `export_vasprun`.
    - ktick_inds : High symmetry kpoints indices.abs
    - ktick_vals  : High Symmetry kpoints labels.
    - elements   : List [[0],], by defualt and plot first ion's projections.
    - orbs       : List [[0],] lists of indices of orbitals, could be empty.
    - labels     : List [str,] of orbitals labels. len(labels)==len(orbs) must hold.  Auto adds `↑`,`↓` for ISPIN=2. If a label is empty i.e. '', it will not show up in legend.
    - colormap  : Matplotlib's standard color maps. Default is 'gist_ranibow'.
    - showlegend : True by defualt and displays legend relative to axes[0]. If False, it writes text on individual ax.
    - scale_data : Default is False, If True, normalize projection data to 1.
    - max_width  : Width to scale whole projections. Default is None and linewidth at any point on a line = 2.5*sum(ions+orbitals projection of the input for that line at that point). Linewidth is scaled to max_width if an int or float is given.
    - figsize    : Tuple (width,height) in inches. Default (3.4.2.6) is article column's width.
    - txt        : Text on figure, if None, SYSTEM's name is printed.
    - xytxt      : [x_coord,y_coord] of text relative to axes.
    - ctxt       : color of text.
    - spin       : Plot spin-polarized for spin {'up','down','both'}. Default is both.
    - interpolate: Default is False, if True, bands are interpolated.
    - n          : int, number of points, default is 5.
    - k          : int, order of interpolation 0,1,2,3. Defualt 3. `n > k` should be hold.
    - legend_kwargs: Dictionary containing legend arguments.
    - **subplots_adjust_kwargs : plt.subplots_adjust parameters.
- **Returns**
    - ax : matplotlib axes object with plotted projected bands.
> Note: Two figures made by this function could be comapred quantitatively only if `scale_data=False, max_width=None` as these parameters act internally on data.



<h4 id="select_pdos" class="doc_header"><code>select_pdos</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L977" class="source_link" style="float:right">[source]</a></h4>

> <code>select_pdos</code>(**`tdos`**=*`None`*, **`pdos_set`**=*`None`*, **`ions`**=*`[0]`*, **`orbs`**=*`[0]`*, **`E_Fermi`**=*`0`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*)

- Returns (interpolated/orginal) enrgy(N,), tdos(N,), and pdos(N,) of selected ions/orbitals.
- **Parameters**
    - tdos     : `export_vasprun`().tdos or `get_tdos`().tdos. If calculations are spin-polarized, it will be `..tdos.SpinUp/SpinDown` for both. You need to apply this function twice for SpinUp and SpinDown separately.
    - pdos_set : `export_vasprun().pro_dos.pros` or `get_dos_pro_set`().pros. If calculations are spin-polarized, it will be `...pros.SpinUp/SpinDown` for both.
    - ions     : List of ions to project on, could be `range(start,stop,step)` as well, remember that `stop` is not included in python. so `range(0,2)` will generate 0 and 1 indices.
    - orbs     : List of orbitals indices to pick.
    - E_Fermi  : Here it is zero. Needs to be input.
    - interpolate: Deafult is false, if True, it will add n points between nearest points.
    - n        : int, default is 5. Adds n points between nearest kpoints.
    - k        : int, order of interpolation, defualt is 3. `n > k` should be hold.



<h4 id="collect_dos" class="doc_header"><code>collect_dos</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L1015" class="source_link" style="float:right">[source]</a></h4>

> <code>collect_dos</code>(**`path_evr`**=*`None`*, **`elim`**=*`[]`*, **`elements`**=*`[[0]]`*, **`orbs`**=*`[[0]]`*, **`labels`**=*`['s']`*, **`E_Fermi`**=*`None`*, **`spin`**=*`'both'`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*)

- Returns lists of energy,tdos, pdos and labels. If given,elements,orbs and labels must have same length. If not given, zeroth ions is collected with s-orbital.
- **Parameters**)
    - path_evr   : Path/to/vasprun.xml or output of `export_vasprun`. Auto picks in CWD.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from `export_vasprun`.
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



<h4 id="quick_dos_lines" class="doc_header"><code>quick_dos_lines</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L1126" class="source_link" style="float:right">[source]</a></h4>

> <code>quick_dos_lines</code>(**`path_evr`**=*`None`*, **`ax`**=*`None`*, **`elim`**=*`[]`*, **`include_dos`**=*`'both'`*, **`elements`**=*`[[0]]`*, **`orbs`**=*`[[0]]`*, **`labels`**=*`['s']`*, **`colormap`**=*`'gist_rainbow'`*, **`tdos_color`**=*`(0.8, 0.95, 0.8)`*, **`linewidth`**=*`0.5`*, **`fill_area`**=*`True`*, **`vertical`**=*`False`*, **`E_Fermi`**=*`None`*, **`figsize`**=*`(3.4, 2.6)`*, **`txt`**=*`None`*, **`xytxt`**=*`[0.2, 0.85]`*, **`ctxt`**=*`'black'`*, **`spin`**=*`'both'`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*, **`showlegend`**=*`True`*, **`legend_kwargs`**=*`{'ncol': 4, 'anchor': (0, 1), 'handletextpad': 0.5, 'handlelength': 1, 'fontsize': 'small', 'frameon': False}`*)

- Returns ax object (if ax!=False) and plot on which all matplotlib allowed actions could be performed, returns lists of energy,tdos and pdos and labels. If given,elements,orbs colors, and labels must have same length. If not given, zeroth ions is plotted with s-orbital.
- **Parameters**)
    - path_evr   : Path/to/vasprun.xml or output of `export_vasprun`. Auto picks in CWD.
    - ax         : Matplotlib axes object, if None, one is created. If False, data lists are returned.
    - include_dos: One of {'both','tdos','pdos'}.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from `export_vasprun`.
    - elements   : List [[0],], by defualt and plot first ion's projections.
    - orbs       : List [[0],] lists of indices of orbitals, could be empty.
    - labels     : List [str,] of orbitals labels. len(labels)==len(orbs) must hold.  Auto adds `↑`,`↓` for ISPIN=2.
    - colormap  : Matplotlib's standard color maps. Default is 'gist_ranibow'. Use 'RGB' if want to compare with `quick_rgb_lines` with 3 projection inputs (len(orbs)==3).
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



<h4 id="plt2html" class="doc_header"><code>plt2html</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L1278" class="source_link" style="float:right">[source]</a></h4>

> <code>plt2html</code>(**`plt_fig`**=*`None`*, **`transparent`**=*`True`*, **`dash_html`**=*`None`*)

- Returns base64 encoded Image to display in notebook or HTML <svg> or plotly's dash_html_components.Img object.
- **Parameters**
    - plt_fig    : Matplotlib's figure instance, auto picks as well.
    - transparent: True of False for fig background.
    - dash_html  : Default is None which results in an image display in jupyter notebook.
        - If True, returns html.Img object for plotly's dash.
        - If False, returns <svg> object to embed in HTML DOM.



<h4 id="plt2text" class="doc_header"><code>plt2text</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L1311" class="source_link" style="float:right">[source]</a></h4>

> <code>plt2text</code>(**`plt_fig`**=*`None`*, **`width`**=*`144`*, **`vscale`**=*`0.96`*, **`colorful`**=*`True`*, **`invert`**=*`False`*, **`crop`**=*`False`*, **`outfile`**=*`None`*)

Displays matplotlib figure in terminal as text. You should use a monospcae font like `Cascadia Code PL` to display image correctly. Use before plt.show().
- **Parameters**
    - plt_fig: Matplotlib's figure instance. Auto picks if not given.
    - width  : Character width in terminal, default is 144. Decrease font size when width increased.
    - vscale : Useful to tweek aspect ratio. Default is 0.96 and prints actual aspect in `Cascadia Code PL`. It is approximately `2*width/height` when you select a single space in terminal.
    - colorful: Default is False, prints colored picture if terminal supports it, e.g Windows Terminal.
    - invert  : Defult is False, could be useful for grayscale image.
    - crop    : Default is False. Crops extra background, can change image color if top left pixel is not in background, in that case set this to False.
    - outfile: If None, prints to screen. Writes on a file.



<h4 id="plot_potential" class="doc_header"><code>plot_potential</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L1369" class="source_link" style="float:right">[source]</a></h4>

> <code>plot_potential</code>(**`basis`**=*`None`*, **`e_or_m`**=*`None`*, **`operation`**=*`'mean_z'`*, **`ax`**=*`None`*, **`period`**=*`None`*, **`lr_pos`**=*`(0.25, 0.75)`*, **`lr_widths`**=*`[0.5, 0.5]`*, **`labels`**=*`('$V(z)$', '$\\langle V \\rangle _{roll}(z)$', '$\\langle V \\rangle $')`*, **`colors`**=*`((0, 0.2, 0.7), 'b', 'r')`*, **`annotate`**=*`True`*)

- Returns tuple(ax,Data) where Data contains resultatnt parameters of averaged potential of LOCPOT.
- **Parameters**
    - basis  : `export_potential().basis`.
    - e_or_m : `epxort_potential().[e,m,m_x,m_y,m_z]` is 3D grid data. As `epxort_potential` is slow, so compute it once and then plot the output data.
    - operation: Default is 'mean_z'. What to do with provided volumetric potential data. Anyone of these 'mean_x','min_x','max_x','mean_y','min_y','max_y','mean_z','min_z','max_z'.
    - ax: Matplotlib axes, if not given auto picks.
    - period: Periodicity of potential in fraction between 0 and 1. For example if a slab is made of 4 super cells in z-direction, period=0.25.
    - lr_pos: Locations around which averages are taken.Default (0.25,0.75). Provide in fraction between 0 and 1. Center of period is located at these given fractions. Work only if period is given.
    - lr_widths: Default is [0.5,0.5], you may have slabs which have different lengths on left and right side. Provide a pair proportional to widths e.g (1,1), (1,1.1) etc. and it is auto normalized to 1. Works only if period is given.
    - labels: List of three labels for legend. Use plt.legend() or pp.add_legend() for labels to appear. First entry is data plot, second is its convolution and third is complete average.
    - colors: List of three colors for lines.
    - annotate: True by default, writes difference of right and left averages on plot.



<h4 id="get_rgb_data" class="doc_header"><code>get_rgb_data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L23" class="source_link" style="float:right">[source]</a></h4>

> <code>get_rgb_data</code>(**`kpath`**=*`None`*, **`evals_set`**=*`None`*, **`pros_set`**=*`None`*, **`elements`**=*`[[0], [], []]`*, **`orbs`**=*`[[0], [], []]`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*, **`scale_color`**=*`False`*)

- Returns a formatted RGB colored data to pass into `rgb2plotly` function. Two arguments, `elements` and `orbs` should be in one-to-one correspondence. Returned item has transpose data shape, so that main iteration is over bands.
- **Parameters**
    - kapath   : `export_vasprun`().kpath or `get_kpts`().kpath.
    - evals_set: `export_vasprun`().bands.evals or `get_evals`().evals. If calculations are spin-polarized, it will be `...evals.SpinUp/SpinDown` for both. You need to apply twice for SpinUp and SpinDown separately.
    - pros_set : `export_vasprun().pro_bands.pros` or `get_bands_pro_set`().pros. If calculations are spin-polarized, it will be `...pros.SpinUp/SpinDown` for both. You need to create collections twice for SpinUp and SpinDown separately.
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
    - widths   : An (NBAND,NKPTS) numpy arry, its actually colors summed along z-axis.



<h4 id="flip_even_patches" class="doc_header"><code>flip_even_patches</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L88" class="source_link" style="float:right">[source]</a></h4>

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



<h4 id="rgb2plotly" class="doc_header"><code>rgb2plotly</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L111" class="source_link" style="float:right">[source]</a></h4>

> <code>rgb2plotly</code>(**`rgb_data`**=*`None`*, **`mode`**=*`'markers'`*, **`max_width`**=*`None`*, **`showlegend`**=*`False`*, **`name`**=*`''`*, **`labels`**=*`['s', 'p', 'd']`*, **`symbol`**=*`0`*, **`start`**=*`0`*)

- Returns data object of plotly's figure using `get_rgb_data`. Returned data could be fed to a plolty's figure.
- ** Parameters**
    - rgb_data    : output of `get_rgb_data`.
    - mode        : Three plotting modes are available:
        - 'markers' : Plot whole data as a single scatter object. Its too fast.
        - 'bands'   : Plot data such that each band is accessible via legend.
        - 'lines'   : A replica of `matplotlib LineCollection` object. It plots at each point separately, slower than other two modes.
    - max_width  : Line/Scatter thickness is scaled to `max_width`. None by default and represent actual data.
    - name       : Name to be shown on hover text or legend.
    - labels     : Optional, show red green blue colors corresponding orbitals.
    - showlegend : Optional, only suitbale if spin up/down or 'bands' mode is ON.
    - symbol     : Plotly's marker symbol. 0 for circle, 5/6 for Up/Down.
    - start      : Start index of bands, defult is zero. Should be export_vasprun().bands.indices[0] or get_evals().indices[0].



<h4 id="plotly2html" class="doc_header"><code>plotly2html</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L178" class="source_link" style="float:right">[source]</a></h4>

> <code>plotly2html</code>(**`fig`**, **`filename`**=*`None`*, **`out_string`**=*`False`*, **`modebar`**=*`True`*)

- Writes plotly's figure as HTML file or display in IPython which is accessible when online. It is different than plotly's `fig.to_html` as it is minimal in memory. If you need to have offline working file, just use `fig.write_html('file.html')` which will be larger in size.
- **Parameters**
    - fig      : A plotly's figure object.
    - filename : Name of file to save fig. Defualt is None and show plot in Colab/Online or return hrml string.
    - out_string: If True, returns HTML string, if False displays graph if possible.



<h4 id="plotly_rgb_lines" class="doc_header"><code>plotly_rgb_lines</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L229" class="source_link" style="float:right">[source]</a></h4>

> <code>plotly_rgb_lines</code>(**`path_evr`**=*`None`*, **`elements`**=*`[[], [], []]`*, **`orbs`**=*`[[], [], []]`*, **`labels`**=*`['', '', '']`*, **`mode`**=*`'markers'`*, **`elim`**=*`[]`*, **`E_Fermi`**=*`None`*, **`skipk`**=*`None`*, **`kseg_inds`**=*`[]`*, **`max_width`**=*`6`*, **`title`**=*`None`*, **`ktick_inds`**=*`[0, -1]`*, **`ktick_vals`**=*`['Γ', 'M']`*, **`figsize`**=*`None`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*)

- Returns plotly's figure object, takes care of spin-polarized calculations automatically. `elements`,`orbs` and `labels` are required to be one-to-one lists of size 3 where each item in list could be another list or integer.
- **Parameters**
    - path_evr  : Path/to/vasprun.xml or xml output of `read_asxml`.
    - elements   : List of size 3 of list of indices of ions. If not given, picks all ions for each orbital.
    - orbs       : List of size 3 of list of orbital indices, if not gievn, s,p,d plotted.
    - labels  : List of labels for projection.
    - mode       : Three plotting modes are available:
        - 'markers' : Plot whole data as a single scatter object. Its too fast.
        - 'bands'   : Plot data such that each band is accessible via legend.
        - 'lines'   : A replica of `matplotlib LineCollection` object. It plots at each point separately, slower than other two modes.
    - **kwargs      : interpolate, ticks, figsize,elim,kseg_inds,max_width,title etc.



<h4 id="plotly_dos_lines" class="doc_header"><code>plotly_dos_lines</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L381" class="source_link" style="float:right">[source]</a></h4>

> <code>plotly_dos_lines</code>(**`path_evr`**=*`None`*, **`elim`**=*`[]`*, **`elements`**=*`[[0]]`*, **`orbs`**=*`[[0]]`*, **`labels`**=*`['s']`*, **`colormap`**=*`'gist_rainbow'`*, **`tdos_color`**=*`(0.5, 0.95, 0)`*, **`linewidth`**=*`2`*, **`fill_area`**=*`True`*, **`vertical`**=*`False`*, **`E_Fermi`**=*`None`*, **`figsize`**=*`None`*, **`spin`**=*`'both'`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*, **`title`**=*`None`*)

- Returns ax object (if ax!=False) and plot on which all matplotlib allowed actions could be performed, returns lists of energy,tdos and pdos and labels. If given,elements,orbs colors, and labels must have same length. If not given, zeroth ions is plotted with s-orbital.
- **Parameters**)
    - path_evr   : Path/to/vasprun.xml or output of `export_vasprun`. Auto picks in CWD.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from `export_vasprun`.
    - elements   : List [[0,],] of ions indices, by defualt plot first ion's projections.
    - orbs       : List [[0,],] lists of indices of orbitals, could be empty.
    - labels     : List [str,] of orbitals labels. len(labels)==len(orbs) must hold.
    - colormap  : Matplotlib's standard color maps. Default is 'gist_ranibow'. Use 'RGB' if want to compare with `plotly_rgb_lines` with 3 projection inputs (len(orbs)==3).
    - fill_area  : Default is True and plots filled area for dos. If False, plots lines only.
    - vertical   : False, If True, plots along y-axis.
    - figsize    : Tuple in pixels (width,height).
    - interpolate: Default is False, if True, bands are interpolated.
    - n          : int, number of points, default is 5.
    - k          : int, order of interpolation 0,1,2,3. Defualt 3. `n > k` should be hold.
    - legend_kwargs: Dictionary to contain legend arguments to fix.
- **Returns**
    - fig        : Plotly's figure object.



<h2 id="Arrow3D" class="doc_header"><code>class</code> <code>Arrow3D</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L36" class="source_link" style="float:right">[source]</a></h2>

> <code>Arrow3D</code>(**`x`**, **`y`**, **`z`**, **`u`**, **`v`**, **`w`**, **\*`args`**, **\*\*`kwargs`**) :: `FancyArrowPatch`

Draw 3D fancy arrow.



<h4 id="fancy_quiver3d" class="doc_header"><code>fancy_quiver3d</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L51" class="source_link" style="float:right">[source]</a></h4>

> <code>fancy_quiver3d</code>(**`X`**, **`Y`**, **`Z`**, **`U`**, **`V`**, **`W`**, **`ax`**=*`None`*, **`C`**=*`'r'`*, **`L`**=*`0.7`*, **`mutation_scale`**=*`10`*, **\*\*`kwargs`**)

Plots 3D arrows on a given ax. See [FancyArrowPatch](https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.patches.FancyArrowPatch.html).
- **Parameters**
    - X, Y, Z : 1D arrays of coordinates of arrows' tail point.
    - U, V, W : 1D arrays of dx,dy,dz of arrows.
    - ax: 3D axes, if not given, auto created.
    - C : 1D colors array mapping for arrows. Could be one color.
    - L : 1D linwidths array mapping for arrows. Could be one linewidth.
    - mutation_scale: Arrow head width/size scale. Default is 10.
    - kwargs: FancyArrowPatch's keyword arguments excluding positions,color, lw and mutation_scale, shrinkA, shrinkB which are already used. An important keyword argument is `arrowstyle` which could be '->','-|>', their inverted forms and many more. See on matplotlib.



<h4 id="save_mp_API" class="doc_header"><code>save_mp_API</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L77" class="source_link" style="float:right">[source]</a></h4>

> <code>save_mp_API</code>(**`api_key`**)

- Save materials project api key for autoload in functions.



<h4 id="load_mp_data" class="doc_header"><code>load_mp_data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L96" class="source_link" style="float:right">[source]</a></h4>

> <code>load_mp_data</code>(**`formula`**, **`api_key`**=*`None`*, **`mp_id`**=*`None`*, **`max_sites`**=*`None`*, **`min_sites`**=*`None`*)

- Returns fetched data using request api of python form materials project website.
- **Parameters**
    - formula  : Material formula such as 'NaCl'.
    - api_key  : API key for your account from material project site. Auto picks if you already used `save_mp_API` function.
    - mp_id     : Optional, you can specify material ID to filter results.
    - max_sites : Maximum number of sites. If None, sets `min_sites + 1`, if `min_sites = None`, gets all data.
    - min_sites : Minimum number of sites. If None, sets `max_sites + 1`, if `max_sites = None`, gets all data.



<h4 id="get_crystal" class="doc_header"><code>get_crystal</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L151" class="source_link" style="float:right">[source]</a></h4>

> <code>get_crystal</code>(**`formula`**, **`api_key`**=*`None`*, **`mp_id`**=*`None`*, **`max_sites`**=*`None`*, **`min_sites`**=*`None`*)

- Returns crystal information dictionary including cif data format.
- **Parameters**
    - formula  : Material formula such as 'NaCl'.
    - api_key  : API key for your account from material project site. Auto picks if you already used `save_mp_API` function.
    - mp_id    : Optional, you can specify material ID to filter results.
    - max_sites : Maximum number of sites. If None, sets `min_sites + 1`, if `min_sites = None`, gets all data.
    - min_sites : Minimum number of sites. If None, sets `max_sites + 1`, if `max_sites = None`, gets all data.



<h4 id="get_poscar" class="doc_header"><code>get_poscar</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L174" class="source_link" style="float:right">[source]</a></h4>

> <code>get_poscar</code>(**`formula`**, **`api_key`**=*`None`*, **`mp_id`**=*`None`*, **`max_sites`**=*`None`*, **`min_sites`**=*`None`*)

- Returns poscar information dictionary including cif data format.
- **Parameters**
    - formula  : Material formula such as 'NaCl'.
    - api_key  : API key for your account from material project site. Auto picks if you already used `save_mp_API` function.
    - mp_id    : Optional, you can specify material ID to filter results.
    - max_sites : Maximum number of sites. If None, sets `min_sites + 1`, if `min_sites = None`, gets all data.
    - min_sites : Minimum number of sites. If None, sets `max_sites + 1`, if `max_sites = None`, gets all data.
- **Usage**
    - `get_poscar('GaAs',api_key,**kwargs)`. Same result is returned from `Get-POSCAR` command in PowerShell terminal if Vasp2Visual module is installed.



<h4 id="get_kpath" class="doc_header"><code>get_kpath</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L245" class="source_link" style="float:right">[source]</a></h4>

> <code>get_kpath</code>(**`hsk_list`**=*`[]`*, **`labels`**=*`[]`*, **`n`**=*`5`*, **`weight`**=*`None`*, **`ibzkpt`**=*`None`*, **`outfile`**=*`None`*)

- Generate list of kpoints along high symmetry path. Options are write to file or return KPOINTS list. It generates uniformly spaced point with input `n` as just a scale factor of number of points per unit length. You can also specify custom number of kpoints in an interval by putting number of kpoints as 4th entry in left kpoint.
- **Parameters**
    - hsk_list : N x 3 list of N high symmetry points, if broken path then [[N x 3],[M x 3],...]. Optionally you can put a 4 values point where 4th entry will decide number of kpoints in current interval. Make sure that points in a connected path patch are at least two i.e. `[[x1,y1,z1],[x2,y2,z2]]` or `[[x1,y1,z1,N],[x2,y2,z2]]`.
    - n        ; int, number per unit length, this makes uniform steps based on distance between points.
    - weight : Float, if None, auto generates weights.
    - ibzkpt : Path to ibzkpt file, required for HSE calculations.
    - labels : Hight symmetry points labels. Good for keeping record of lables and points indices for later use.                - Note: If you do not want to label a point, label it as 'skip' at its index and it will be removed.
    - outfile: Path/to/file to write kpoints.
- **Attributes**
    - If `outfile = None`, a tuple is returned which consists of:
        - nkpts   : get_kmesh().nkpts.
        - kpoints : get_kmesh().kpoints.
        - weights : get_kmesh().weights.



<h4 id="export_poscar" class="doc_header"><code>export_poscar</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L313" class="source_link" style="float:right">[source]</a></h4>

> <code>export_poscar</code>(**`path`**=*`None`*)

Export POSCAR file to python objects.
- **Parameters**
    - path: Path/to/POSCAR file. Auto picks in CWD.



<h4 id="get_basis" class="doc_header"><code>get_basis</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L349" class="source_link" style="float:right">[source]</a></h4>

> <code>get_basis</code>(**`path_pos`**)

Returns given(computed) and inverted(without 2π) basis as tuple(given,inverted).
- **Parameters**
    - path_pos: path/to/POSCAR or 3 given vectors as rows of a matrix.



<h4 id="get_kmesh" class="doc_header"><code>get_kmesh</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L365" class="source_link" style="float:right">[source]</a></h4>

> <code>get_kmesh</code>(**`n_xyz`**=*`[5, 5, 5]`*, **`weight`**=*`None`*, **`ibzkpt`**=*`None`*, **`path_pos`**=*`None`*, **`outfile`**=*`None`*)

- Generates uniform mesh of kpoints. Options are write to file, or return KPOINTS list.
- **Parameters**
    - n_xyz  : List of [nx ny nz] or integer. If integer given, it represents numbers of kpoints along smallest side in reciprocal space and kmesh is autoscaled.
    - weight : Float, if None, auto generates weights.
    - ibzkpt : Path to ibzkpt file, required for HSE calculations.
    - path_pos : POSCAR file path or real space lattice vectors, if None, cubic symmetry is used and it is fast.
    - outfile: Path/to/file to write kpoints.
- **Attributes**
    - If `outfile = None`, a tuple is returned which consists of:
        - nkpts   : get_kmesh().nkpts.
        - kpoints : get_kmesh().kpoints.
        - weight  : get_kmesh().weight, its one float number, provided or calculated.



<h4 id="tan_inv" class="doc_header"><code>tan_inv</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L427" class="source_link" style="float:right">[source]</a></h4>

> <code>tan_inv</code>(**`vy`**, **`vx`**)

- Returns full angle from x-axis counter clockwise.
- **Parameters**
    - vy : Perpendicular componet of vector including sign.
    - vx : Base compoent of vector including sign.



<h4 id="order" class="doc_header"><code>order</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L457" class="source_link" style="float:right">[source]</a></h4>

> <code>order</code>(**`points`**, **`loop`**=*`True`*)

- Returns indices of counterclockwise ordered vertices of a plane in 3D.
- **Parameters**
    - points: numpy array of shape (N,3) or List[List(len=3)].
    - loop  : Default is True and appends start point at end to make a loop.
- **Example**
    > pts = np.array([[1,0,3],[0,0,0],[0,1,2]])
    > inds = order(pts)
    > pts[inds]
    ```
    array([[1, 2, 3],
           [0, 0, 0],
           [1, 0, 3]
           [0, 1, 2]])
    ```



<h4 id="out_bz_plane" class="doc_header"><code>out_bz_plane</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L502" class="source_link" style="float:right">[source]</a></h4>

> <code>out_bz_plane</code>(**`test_point`**, **`plane`**)

- Returns True if test_point is between plane and origin. Could be used to sample BZ mesh in place of ConvexHull.
- **Parameters**
    - test_points: 3D point.
    - plane      : List of at least three coplanar 3D points.



<h4 id="rad_angle" class="doc_header"><code>rad_angle</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L519" class="source_link" style="float:right">[source]</a></h4>

> <code>rad_angle</code>(**`v1`**, **`v2`**)

- Returns interier angle between two vectors.
- **Parameters**
    - v1,v2 : Two vectors/points in 3D.



<h4 id="get_bz" class="doc_header"><code>get_bz</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L534" class="source_link" style="float:right">[source]</a></h4>

> <code>get_bz</code>(**`path_pos`**=*`None`*, **`loop`**=*`True`*, **`digits`**=*`8`*, **`primitive`**=*`False`*)

- Return required information to construct first Brillouin zone in form of tuple (basis, normals, vertices, faces).
- **Parameters**
    - path_pos : POSCAR file path or list of 3 vectors in 3D aslist[list,list,list].
    - loop   : If True, joins the last vertex of a BZ plane to starting vertex in order to complete loop.
    - digits : int, rounding off decimal places, no effect on intermediate calculations, just for pretty final results.
    - primitive: Defualt is False and returns Wigner-Seitz cell, If True returns parallelipiped in rec_basis.
- **Attributes**
    - basis   : get_bz().basis, recprocal lattice basis.
    - normals : get_bz().normals, all vertors that are perpendicular BZ faces/planes.
    - vertices: get_bz().vertices, all vertices of BZ, could be input into ConvexHull for constructing 3D BZ.
    - faces   : get_bz().faces, vertices arranged into faces, could be input to Poly3DCollection of matplotlib for creating BZ from faces' patches.
    - specials : get_bz().specials, Data with attributes `coords`,`kpoints` and `near` in on-one correspondence for high symmetry KPOINTS in recirprocal coordinates space. `near` gives indices of nearest special points around a vertex. All vertices with z > 0 are included.



<h4 id="splot_bz" class="doc_header"><code>splot_bz</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L630" class="source_link" style="float:right">[source]</a></h4>

> <code>splot_bz</code>(**`path_pos_bz`**=*`None`*, **`ax`**=*`None`*, **`plane`**=*`None`*, **`color`**=*`'blue'`*, **`fill`**=*`True`*, **`vectors`**=*`True`*, **`v3`**=*`False`*, **`vname`**=*`'b'`*, **`colormap`**=*`'plasma'`*, **`light_from`**=*`(1, 1, 1)`*, **`alpha`**=*`0.4`*)

- Plots matplotlib's static figure.
- **Parameters**
    - path_pos_bz: Auto picks in CWD if POSCAR file found. This accept three kind of objects:
        - List of 3 basis vectors in real space.
        - Path/to/POSCAR.
        - Output of `get_bz` function.

    - fill       : True by defult, determines whether to fill surface of BZ or not.
    - color      : color to fill surface and stroke color.
    - vectors    : Plots basis vectors, default is True.
    - v3         : Plots 3rd vector as well. Only works in 2D and when `vectors=True`.
    - plane      : Default is None and plots 3D surface. Can take 'xy','yz','zx' to plot in 2D.
    - ax         : Auto generated by default, 2D/3D axes, auto converts in 3D on demand as well.
    - vname      : Default is `b` for reciprocal space, can set `a` for plotting cell as after `get_bz(get_bz().basis)` you get real space lattice back if `primitive=True` both times.
    - colormap  : If None, single color is applied, only works in 3D and `fill=True`. Colormap is applied along z.
    - light_from: Point from where light is thrown on BZ planes, default is (1,1,1). Only works on plane in 3D.
    - alpha    : Opacity of filling in range [0,1]. Increase for clear viewpoint.
- **Returns**
    - ax   : Matplotlib's 2D axes if `plane=None`.
    - ax3d : Matplotlib's 2D axes if `plane` is given.

> Tip: `splot_bz(rec_basis,primitive=True)` will plot cell in real space.



<h4 id="iplot_bz" class="doc_header"><code>iplot_bz</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L747" class="source_link" style="float:right">[source]</a></h4>

> <code>iplot_bz</code>(**`path_pos_bz`**=*`None`*, **`fill`**=*`True`*, **`color`**=*`'rgba(168,204,216,0.4)'`*, **`background`**=*`'rgb(255,255,255)'`*, **`vname`**=*`'b'`*, **`alpha`**=*`0.4`*, **`ortho3d`**=*`True`*, **`fig`**=*`None`*)

- Plots interactive figure showing axes,BZ surface, special points and basis, each of which could be hidden or shown.
- **Parameters**
    - path_pos_bz: Auto picks in CWD if POSCAR file found. This accept three kind of objects:
        - List of 3 basis vectors in real space.
        - Path/to/POSCAR.
        - Output of `get_bz` function.

    - fill       : True by defult, determines whether to fill surface of BZ or not.
    - color      : color to fill surface 'rgba(168,204,216,0.4)` by default.
    - background : Plot background color, default is 'rgb(255,255,255)'.
    - vname      : Default is `b` for reciprocal space, can set `a` for plotting cell as after `get_bz(get_bz().basis)` you get real space lattice back if `primitive=True` both times.
    - alpha      : Opacity of BZ planes.
    - ortho3d    : Default is True, decides whether x,y,z are orthogonal or perspective.
    - fig        : (Optional) Plotly's `go.Figure`. If you want to plot on another plotly's figure, provide that.
- **Returns**
    - fig   : plotly.graph_object's Figure instance.

> Tip: `iplot_bz(rec_basis,primitive=True)` will plot cell in real space.



<h4 id="to_R3" class="doc_header"><code>to_R3</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L854" class="source_link" style="float:right">[source]</a></h4>

> <code>to_R3</code>(**`basis`**, **`points`**)

Transforms coordinates of points (relative to non-othogonal basis) into orthogonal space.
- **Parameters**
    - basis : Non-orthogonal basis of real or reciprocal space.
    - points: 3D points relative to basis, such as KPOINTS and Lattice Points.



<h4 id="kpoints2bz" class="doc_header"><code>kpoints2bz</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L869" class="source_link" style="float:right">[source]</a></h4>

> <code>kpoints2bz</code>(**`bz`**, **`kpoints`**, **`primitive`**=*`False`*)

Brings KPOINTS inside BZ. Applies `to_R3` only if `primitive=True`.
- **Parameters**
    - bz       : Output of get_bz(), make sure use same value of `primitive` there and here.
    - kpoints  : List or array of KPOINTS to transorm into BZ or R3.
    - primitive: Default is False and brings kpoints into regular BZ. If True, returns `to_R3()`.



<h2 id="BZ" class="doc_header"><code>class</code> <code>BZ</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L911" class="source_link" style="float:right">[source]</a></h2>

> <code>BZ</code>(**`path_pos`**=*`None`*, **`loop`**=*`True`*, **`digits`**=*`8`*, **`primitive`**=*`False`*)

- Return required information to construct first Brillouin zone in form of tuple (basis, normals, vertices, faces).
- **Parameters**
    - path_pos : POSCAR file path or list of 3 vectors in 3D aslist[list,list,list].
    - loop   : If True, joins the last vertex of a BZ plane to starting vertex in order to complete loop.
    - digits : int, rounding off decimal places, no effect on intermediate calculations, just for pretty final results.
    - primitive: Defualt is False and returns Wigner-Seitz cell, If True returns parallelipiped in rec_basis.
- **Methods**
    - iplot/iplotc : Interactive plot of BZ/Cell.
    - splot/splotc : Static plot of BZ/Cell.
    - fetch/fetchc : Brings KPOINTS/Positions inside BZ/Cell.



<h4 id="fix_sites" class="doc_header"><code>fix_sites</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L955" class="source_link" style="float:right">[source]</a></h4>

> <code>fix_sites</code>(**`poscar`**, **`tol`**=*`0.01`*, **`eqv_sites`**=*`True`*)

Add equivalent sites to make a full data shape of lattice. Returns same data after fixing.
- **Parameters**
    - poscar: Output of `export_poscar` or `export_vasprun().poscar`.
    - tol   : Tolerance value. Default is 0.01.
    - eqv_sites: If True, add sites on edges and faces. If False, just fix coordinates, i.e. `pos > 1 - tol -> pos - 1`, useful for merging poscars to make slabs.



<h4 id="get_pairs" class="doc_header"><code>get_pairs</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L1004" class="source_link" style="float:right">[source]</a></h4>

> <code>get_pairs</code>(**`basis`**, **`positions`**, **`r`**, **`eps`**=*`0.01`*)

Returns a tuple of Lattice (coords,pairs), so coords[pairs] given nearest site bonds.
- **Parameters**
    - basis: Real space lattice basis.
    - positions: Array(N,3) of fractional positions of lattice sites. If coordinates positions, provide unity basis.
    - r        : Cartesian distance between the pairs in units of Angstrom e.g. 1.2 -> 1.2E-10.
    - eps      : Tolerance value. Default is 10^-2.



<h4 id="iplot_lat" class="doc_header"><code>iplot_lat</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L1028" class="source_link" style="float:right">[source]</a></h4>

> <code>iplot_lat</code>(**`poscar`**, **`sizes`**=*`10`*, **`colors`**=*`'blue'`*, **`bond_length`**=*`None`*, **`tol`**=*`0.1`*, **`eps`**=*`0.01`*, **`eqv_sites`**=*`True`*, **`line_width`**=*`4`*, **`edge_color`**=*`'black'`*, **`fill`**=*`False`*, **`alpha`**=*`0.4`*, **`ortho3d`**=*`True`*, **`fig`**=*`None`*)

Interactive plot of lattice.
- **Main Parameters**
    - poscar     : Output of export_poscar or export_vasprun().poscar.
    - sizes      : Size of sites. Either one int/float or list equal to type of ions.
    - colors     : Colors of sites. Either one colors or list equal to type of ions.
    - bond_length: Length of bond in fractional unit [0,1]. It is scaled to V^1/3 and auto calculated if not provides.
Other parameters just mean what they seem to be.



<h4 id="splot_lat" class="doc_header"><code>splot_lat</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L1104" class="source_link" style="float:right">[source]</a></h4>

> <code>splot_lat</code>(**`poscar`**, **`sizes`**=*`50`*, **`colors`**=*`[]`*, **`colormap`**=*`None`*, **`bond_length`**=*`None`*, **`tol`**=*`0.1`*, **`eps`**=*`0.01`*, **`eqv_sites`**=*`True`*, **`line_width`**=*`1`*, **`edge_color`**=*`(1, 0.5, 0, 0.4)`*, **`vectors`**=*`True`*, **`v3`**=*`False`*, **`plane`**=*`None`*, **`light_from`**=*`(1, 1, 1)`*, **`fill`**=*`False`*, **`alpha`**=*`0.4`*, **`ax`**=*`None`*)

Static plot of lattice.
- **Main Parameters**
    - poscar     : Output of export_poscar or export_vasprun().poscar.
    - sizes      : Size of sites. Either one int/float or list equal to type of ions.
    - bond_length: Length of bond in fractional unit [0,1]. It is scaled to V^1/3 and auto calculated if not provides.
    - colors: List of colos. If given, preffered over colormap, should have same length as type of ions.
Other parameters just mean what they seem to be.

> Tip: Use `plt.style.use('ggplot')` for better 3D perception.



<h4 id="join_poscars" class="doc_header"><code>join_poscars</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L1186" class="source_link" style="float:right">[source]</a></h4>

> <code>join_poscars</code>(**`poscar1`**, **`poscar2`**, **`direction`**=*`'z'`*, **`tol`**=*`0.01`*)

Joins two POSCARs in a given direction. In-plane lattice parameters are kept from `poscar1` and basis of `poscar2` parallel to `direction` is modified while volume is kept same.
- **Parameters**
    - poscar1, poscar2:  Base and secondary POSCARs respectivly. Output of `export_poscar` or similar object from other functions.
    - direction: The joining direction. It is general and can join in any direction along basis. Expect one of ['a','b','c','x','y','z'].
    - tol: Default is 0.01. It is used to bring sites near 1 to near zero in order to complete sites in plane. Vasp relaxation could move a point, say at 0.00100 to 0.99800 which is not useful while merging sites.



<h4 id="write_poscar" class="doc_header"><code>write_poscar</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L1258" class="source_link" style="float:right">[source]</a></h4>

> <code>write_poscar</code>(**`poscar`**, **`sd_list`**=*`None`*, **`outfile`**=*`None`*, **`overwrite`**=*`False`*)

Writes poscar data object to a file or returns string
- **Parameters**
    - poscar   : Output of `export_poscar`,`join_poscars` etc.
    - sd_list  : A list ['T T T','F F F',...] strings to turn on selective dynamics at required sites. len(sd_list)==len(sites) should hold.
    - outfile  : str,file path to write on.
    - overwrite: bool, if file already exists, overwrite=True changes it.



<h4 id="scale_poscar" class="doc_header"><code>scale_poscar</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L1296" class="source_link" style="float:right">[source]</a></h4>

> <code>scale_poscar</code>(**`path_poscar`**, **`scale`**=*`(1, 1, 1)`*, **`tol`**=*`0.01`*)

Create larger/smaller cell from a given POSCAR.
- **Parameters**
    - path_poscar: Path/to/POSCAR or `poscar` data object.
    - scale: Tuple of three values along (a,b,c) vectors. int or float values. If number of sites are not as expected in output, tweak `tol` instead of `scale`. You can put a minus sign with `tol` to get more sites and plus sign to reduce sites.
    - tol: It is used such that site positions are blow `1 - tol`, as 1 belongs to next cell, not previous one.
> Tip: scale = (2,2,2) enlarges a cell and next operation of (1/2,1/2,1/2) should bring original cell back.



<h4 id="css_style" class="doc_header"><code>css_style</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L35" class="source_link" style="float:right">[source]</a></h4>

> <code>css_style</code>(**`colors_dict`**)

Return style based on colors_dict available as pp.light_colors, pp.dark_colors etc



<h4 id="{'table_fg': '#ABB2BF', 'tr_odd_bg': '#282C34', 'tr_even_bg': '#21252B', 'tr_hover_bg': '#414855', 'btn_fg': '#ABB2BF', 'btn_bg': '#3D4450', 'tab_bg': '#21252B', 'tab_fg': '#61AFEF', 'tab_shadow': '#282C34', 'box_bg': '#21252B', 'box_border': '#282C34', 'text': '#ABB2BF', 'input_bg': '#282C34', 'input_fg': '#ABB2BF', 'input_border': '#282C34', 'input_hover': '#414855', 'dd_select_fg': '#ABB2BF', 'dd_select_bg': 'transparent', 'dd_hover': '#3D4450', 'dd_opt_bg': '#282C34', 'dd_opt_fg': 'whitesmoke', 'dd_focus_bg': 'green'}" class="doc_header"><code>{'table_fg': '#ABB2BF', 'tr_odd_bg': '#282C34', 'tr_even_bg': '#21252B', 'tr_hover_bg': '#414855', 'btn_fg': '#ABB2BF', 'btn_bg': '#3D4450', 'tab_bg': '#21252B', 'tab_fg': '#61AFEF', 'tab_shadow': '#282C34', 'box_bg': '#21252B', 'box_border': '#282C34', 'text': '#ABB2BF', 'input_bg': '#282C34', 'input_fg': '#ABB2BF', 'input_border': '#282C34', 'input_hover': '#414855', 'dd_select_fg': '#ABB2BF', 'dd_select_bg': 'transparent', 'dd_hover': '#3D4450', 'dd_opt_bg': '#282C34', 'dd_opt_fg': 'whitesmoke', 'dd_focus_bg': 'green'}</code><a href="" class="source_link" style="float:right">[source]</a></h4>

dict() -> new empty dictionary
dict(mapping) -> new dictionary initialized from a mapping object's
    (key, value) pairs
dict(iterable) -> new dictionary initialized as if via:
    d = {}
    for k, v in iterable:
        d[k] = v
dict(**kwargs) -> new dictionary initialized with the name=value pairs
    in the keyword argument list.  For example:  dict(one=1, two=2)



<h4 id="{'table_fg': 'black', 'tr_odd_bg': '#eaf0f0', 'tr_even_bg': 'white', 'tr_hover_bg': '#abe4ff', 'btn_fg': 'black', 'btn_bg': '#c3d4d4', 'tab_bg': '#F3F3F3', 'tab_fg': 'black', 'tab_shadow': 'gray', 'box_bg': '#F3F3F3', 'box_border': 'whitesmoke', 'text': 'black', 'input_bg': 'white', 'input_fg': 'gray', 'input_border': '#e0e8e8', 'input_hover': 'skyblue', 'dd_select_fg': 'skyblue', 'dd_select_bg': 'transparent', 'dd_hover': 'white', 'dd_opt_bg': '#eaf0f0', 'dd_opt_fg': 'gray', 'dd_focus_bg': 'red'}" class="doc_header"><code>{'table_fg': 'black', 'tr_odd_bg': '#eaf0f0', 'tr_even_bg': 'white', 'tr_hover_bg': '#abe4ff', 'btn_fg': 'black', 'btn_bg': '#c3d4d4', 'tab_bg': '#F3F3F3', 'tab_fg': 'black', 'tab_shadow': 'gray', 'box_bg': '#F3F3F3', 'box_border': 'whitesmoke', 'text': 'black', 'input_bg': 'white', 'input_fg': 'gray', 'input_border': '#e0e8e8', 'input_hover': 'skyblue', 'dd_select_fg': 'skyblue', 'dd_select_bg': 'transparent', 'dd_hover': 'white', 'dd_opt_bg': '#eaf0f0', 'dd_opt_fg': 'gray', 'dd_focus_bg': 'red'}</code><a href="" class="source_link" style="float:right">[source]</a></h4>

dict() -> new empty dictionary
dict(mapping) -> new dictionary initialized from a mapping object's
    (key, value) pairs
dict(iterable) -> new dictionary initialized as if via:
    d = {}
    for k, v in iterable:
        d[k] = v
dict(**kwargs) -> new dictionary initialized with the name=value pairs
    in the keyword argument list.  For example:  dict(one=1, two=2)



<h4 id="{'table_fg': 'black', 'tr_odd_bg': 'whitesmoke', 'tr_even_bg': 'white', 'tr_hover_bg': 'lightgray', 'btn_fg': 'black', 'btn_bg': 'whitesmoke', 'tab_bg': 'white', 'tab_fg': 'black', 'tab_shadow': 'white', 'box_bg': 'white', 'box_border': 'white', 'text': 'black', 'input_bg': 'white', 'input_fg': 'black', 'input_border': 'lightgray', 'input_hover': 'gray', 'dd_select_fg': 'black', 'dd_select_bg': 'transparent', 'dd_hover': 'whitesmoke', 'dd_opt_bg': 'white', 'dd_opt_fg': 'black', 'dd_focus_bg': 'white'}" class="doc_header"><code>{'table_fg': 'black', 'tr_odd_bg': 'whitesmoke', 'tr_even_bg': 'white', 'tr_hover_bg': 'lightgray', 'btn_fg': 'black', 'btn_bg': 'whitesmoke', 'tab_bg': 'white', 'tab_fg': 'black', 'tab_shadow': 'white', 'box_bg': 'white', 'box_border': 'white', 'text': 'black', 'input_bg': 'white', 'input_fg': 'black', 'input_border': 'lightgray', 'input_hover': 'gray', 'dd_select_fg': 'black', 'dd_select_bg': 'transparent', 'dd_hover': 'whitesmoke', 'dd_opt_bg': 'white', 'dd_opt_fg': 'black', 'dd_focus_bg': 'white'}</code><a href="" class="source_link" style="float:right">[source]</a></h4>

dict() -> new empty dictionary
dict(mapping) -> new dictionary initialized from a mapping object's
    (key, value) pairs
dict(iterable) -> new dictionary initialized as if via:
    d = {}
    for k, v in iterable:
        d[k] = v
dict(**kwargs) -> new dictionary initialized with the name=value pairs
    in the keyword argument list.  For example:  dict(one=1, two=2)



<h4 id="get_files_gui" class="doc_header"><code>get_files_gui</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L215" class="source_link" style="float:right">[source]</a></h4>

> <code>get_files_gui</code>(**`auto_fill`**=*`'vasprun.xml'`*, **`html_style`**=*`None`*, **`height`**=*`320`*)

- Creates a GUI interface for files/folders filtering.
- **Parmeters**
    - auto_fill  : Default is `vasprun.xml`, any file/folder.
    - html_style : None,Output of `css_style`.
    - height     : Height of Grid box.
- **Returns**
    - Tuple(GUI_gridbox,Files_Dropdown). Access second one by item itself.



<h2 id="InputGui" class="doc_header"><code>class</code> <code>InputGui</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L309" class="source_link" style="float:right">[source]</a></h2>

> <code>InputGui</code>(**`sys_info`**=*`None`*, **`html_style`**=*`None`*, **`height`**=*`400`*)





<h4 id="click_data" class="doc_header"><code>click_data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L422" class="source_link" style="float:right">[source]</a></h4>

> <code>click_data</code>(**`sel_en_w`**, **`fermi_w`**, **`data_dict`**, **`fig`**, **`bd_w`**)





<h4 id="tabulate_data" class="doc_header"><code>tabulate_data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L459" class="source_link" style="float:right">[source]</a></h4>

> <code>tabulate_data</code>(**`data_dict`**)





<h4 id="save_data" class="doc_header"><code>save_data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L485" class="source_link" style="float:right">[source]</a></h4>

> <code>save_data</code>(**`out_w1`**, **`data_dict`**)





<h4 id="color_toggle" class="doc_header"><code>color_toggle</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L491" class="source_link" style="float:right">[source]</a></h4>

> <code>color_toggle</code>(**`tog_w`**, **`fig`**, **`rd_btn`**)





<h4 id="generate_summary" class="doc_header"><code>generate_summary</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L516" class="source_link" style="float:right">[source]</a></h4>

> <code>generate_summary</code>(**`paths_list`**=*`None`*)





<h2 id="VasprunApp" class="doc_header"><code>class</code> <code>VasprunApp</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/widgets.py#L563" class="source_link" style="float:right">[source]</a></h2>

> <code>VasprunApp</code>(**`height`**=*`580`*)

Display a GUI for vasp output analysis. `self.theme_colors` can be used to edit custom theme.
- **Usage Example**

```python
import pivotpy as pp
va = pp.VasprunApp()
va.cache_data = False #Turn off cache globally.
va.evr_kws['elim'] = [-2,2] #Only Bands in this range will be included. Global accross project, can change anytime.
va.evr_kws['try_pwsh'] = False #Defult is True. Tries to load Powershell exported data.
va.ibands_kws['mode'] = 'bands' #Change graph mode from 'markers' to 'bands'. Setting it to 'lines' is not recommended in live graph, it could hang all UI.
va.show() #Displays App and do work!
va.theme_colors = pp.dark_colors #Set theme to dark externally and edit dictionary values to make your own theme
va.splot(**kwargs) #Get matplotlib plot of current data.
va.df #After you do some analysis and hit `Project Summary` button, get DataFrame.
va.fig #Get current fig in Notebook cell.
```



```python

```
