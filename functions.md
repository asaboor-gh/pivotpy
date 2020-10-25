# Pivotpy Functions Reference

<h2 id="Dic2Dot" class="doc_header"><code>class</code> <code>Dic2Dot</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L8" class="source_link" style="float:right">[source]</a></h2>

> <code>Dic2Dot</code>() :: `dict`

- Returns dot notation accessible if a dictionary is input.
- It is used to pack all functions in a dictionary.



<h4 id="read_asxml" class="doc_header"><code>read_asxml</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L19" class="source_link" style="float:right">[source]</a></h4>

> <code>read_asxml</code>(**`path`**=*`None`*, **`suppress_warning`**=*`False`*)

- Reads a big vasprun.xml file into memory once and then apply commands.
If current folder contains `vasprun.xml` file, it automatically picks it.

- **Parameters**
    - path             : Path/To/vasprun.xml
    - suppress_warning : False by defualt. Warns about memory usage for large files > 100 MB.
- **Returns**
    - xml_data : Xml object to use in other functions



<h4 id="exclude_kpts" class="doc_header"><code>exclude_kpts</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L69" class="source_link" style="float:right">[source]</a></h4>

> <code>exclude_kpts</code>(**`xml_data`**=*`None`*)

- Returns number of kpoints to exclude used from IBZKPT.
- **Parameters**
    - xml_data : From `read_asxml` function
- **Returns**
    - int      : Number of kpoints to exclude.



<h4 id="get_ispin" class="doc_header"><code>get_ispin</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L90" class="source_link" style="float:right">[source]</a></h4>

> <code>get_ispin</code>(**`xml_data`**=*`None`*)

- Returns value of ISPIN.
- **Parameters**
    - xml_data : From `read_asxml` function
- **Returns**
    - int      : Value of ISPIN.



<h4 id="get_summary" class="doc_header"><code>get_summary</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L107" class="source_link" style="float:right">[source]</a></h4>

> <code>get_summary</code>(**`xml_data`**=*`None`*)

- Returns overview of system parameters.
- **Parameters**
    - xml_data : From `read_asxml` function
- **Returns**
    - dict     : Dictionary that contains system information.



<h4 id="get_kpts" class="doc_header"><code>get_kpts</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L142" class="source_link" style="float:right">[source]</a></h4>

> <code>get_kpts</code>(**`xml_data`**=*`None`*, **`skipk`**=*`0`*, **`joinPathAt`**=*`[]`*)

- Returns kpoints and calculated kpath.
- **Parameters**
    - xml_data   : From `read_asxml` function
    - skipk      : Number of initil kpoints to skip
    - joinPathAt : List of indices of kpoints where path is broken
- **Returns**
    - dict     : Dictionary that contains kpoints ans kpath.



<h4 id="get_tdos" class="doc_header"><code>get_tdos</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L176" class="source_link" style="float:right">[source]</a></h4>

> <code>get_tdos</code>(**`xml_data`**=*`None`*, **`spin_set`**=*`1`*, **`elim`**=*`[]`*)

- Returns total dos for a spin_set (default 1) and energy limit. If spin-polarized calculations, gives SpinUp and SpinDown keys as well.
- **Parameters**
    - xml_data : From `read_asxml` function
    - spin_set : int, default is 1.and
    - elim     : List [min,max] of energy, default empty.
- **Returns**
    - dict     : Dictionary that contains E_Fermi, ISPIN,tdos.



<h4 id="get_evals" class="doc_header"><code>get_evals</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L229" class="source_link" style="float:right">[source]</a></h4>

> <code>get_evals</code>(**`xml_data`**=*`None`*, **`skipk`**=*`None`*, **`elim`**=*`[]`*)

- Returns eigenvalues as numpy array. If spin-polarized calculations, gives SpinUp and SpinDown keys as well.
- **Parameters**
    - xml_data : From `read_asxml` function
    - skipk    : Number of initil kpoints to skip.
    - elim     : List [min,max] of energy, default empty.
- **Returns**
    - dict     : Dictionary that contains evals and related parameters.



<h4 id="get_bands_pro_set" class="doc_header"><code>get_bands_pro_set</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L282" class="source_link" style="float:right">[source]</a></h4>

> <code>get_bands_pro_set</code>(**`xml_data`**=*`None`*, **`spin_set`**=*`1`*, **`skipk`**=*`0`*, **`bands_range`**=*`None`*)

- Returns bands projection of a spin_set(default 1) as numpy array. If spin-polarized calculations, gives SpinUp and SpinDown keys as well.
- **Parameters**
    - xml_data    : From `read_asxml` function
    - skipk       : Number of initil kpoints to skip (Default 0).
    - spin_set    : Spin set to get, default is 1.
    - bands_range : If elim used in `get_evals`,that will return bands_range to use here..
- **Returns**
    - dict     : Dictionary that contains bands projections and related parameters.



<h4 id="get_dos_pro_set" class="doc_header"><code>get_dos_pro_set</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L335" class="source_link" style="float:right">[source]</a></h4>

> <code>get_dos_pro_set</code>(**`xml_data`**=*`None`*, **`spin_set`**=*`1`*, **`dos_range`**=*`None`*)

- Returns dos projection of a spin_set(default 1) as numpy array. If spin-polarized calculations, gives SpinUp and SpinDown keys as well.
- **Parameters**
    - xml_data    : From `read_asxml` function
    - spin_set    : Spin set to get, default 1.
    - dos_range   : If elim used in `get_tdos`,that will return dos_range to use here..
- **Returns**
    - dict        : Dictionary that contains dos projections and related parameters.



<h4 id="get_structure" class="doc_header"><code>get_structure</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L378" class="source_link" style="float:right">[source]</a></h4>

> <code>get_structure</code>(**`xml_data`**=*`None`*)

- Returns structure's volume,basis,positions and rec-basis.
- **Parameters**
    - xml_data : From `read_asxml` function
- **Returns**
    - dict     : Dictionary that contains volume,basis,positions and rec_basis.



<h4 id="export_vasprun" class="doc_header"><code>export_vasprun</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L406" class="source_link" style="float:right">[source]</a></h4>

> <code>export_vasprun</code>(**`path`**=*`None`*, **`skipk`**=*`None`*, **`elim`**=*`[]`*, **`joinPathAt`**=*`[]`*, **`shift_kpath`**=*`0`*)

- Returns a full dictionary of all objects from `vasprun.xml` file.
- **Parameters**
    - path       : Path to `vasprun.xml` file. Default is `'./vasprun.xml'`.
    - skipk      : Default is None. Automatically detects kpoints to skip.
    - elim       : List [min,max] of energy interval. Default is [], covers all bands.
    - joinPathAt : List of indices of kpoints where path is broken.
    - shift_kpath: Default 0. Can be used to merge multiple calculations on single axes side by side.
- **Returns**
    - dict : Dictionary accessible via dot notation containing objects:
        - sys_info  : System Information
        - dim_info  : Contains information about dimensions of returned objects.
        - kpoints   : numpy array of kpoints with excluded IBZKPT points
        - kpath     : 1D numpy array directly accessible for plot.
        - bands     : Dictionary containing bands.
        - tdos      : Dictionary containing total dos.
        - pro_bands : Dictionary containing bands projections.
        - pro_dos   : Dictionary containing dos projections.
        - poscar    : containing basis,positions, rec_basis and volume.
        - xml       : xml root object which is iterable over nodes using xml.iter('node').



<h4 id="load_export" class="doc_header"><code>load_export</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L481" class="source_link" style="float:right">[source]</a></h4>

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
    - dict : Dictionary accessible via dot notation containing objects:
        - sys_info  : System Information
        - dim_info  : Contains information about dimensions of returned objects.
        - kpoints   : numpy array of kpoints with excluded IBZKPT points
        - kpath     : 1D numpy array directly accessible for plot.
        - bands     : Dictionary containing bands.
        - tdos      : Dictionary containing total dos.
        - pro_bands : Dictionary containing bands projections.
        - pro_dos   : Dictionary containing dos projections.
        - poscar    : containing basis,positions, rec_basis and volume.



<h4 id="make_dot_dict" class="doc_header"><code>make_dot_dict</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L631" class="source_link" style="float:right">[source]</a></h4>

> <code>make_dot_dict</code>(**`json_loaded`**)

- Returns a pivotpy.Dic2Dot object recursively and makes keys accessible via dot notation. Works only upto 4 nesting levels because it is basically created for vasprun data transfer back and forth in `pivotpy-dash` app.
- **Parameters**
    - json_load : Output of json.load/json.loads or any python dictionary.



<h4 id="dump_dict" class="doc_header"><code>dump_dict</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/vr_parser.py#L657" class="source_link" style="float:right">[source]</a></h4>

> <code>dump_dict</code>(**`dict_obj`**=*`None`*, **`dump_to`**=*`'pickle'`*, **`outfile`**=*`None`*)

- Dump an `export_vasprun` or `load_export` object to json or pickle string/file.
- **Parameters**
    - dict_obj : Any dictionary containg numpy arrays, including `export_vasprun` or `load_export` output.
    - dump_to  : Defualt is `pickle` or `json`.
    - outfile  : Defualt is None and return string. File name does not require extension.



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

> <code>select_dirs</code>(**`path`**=*`'C:\\Users\\mass_\\Desktop'`*, **`include`**=*`[]`*, **`exclude`**=*`[]`*)

- Returns selected directories recursively from a parent directory.
- **Parameters**
    - path    : path to a parent directory, default is `"."`
    - include : list of keywords to include directories, avoid wildcards.
    - exclude : list of keywords to exclude directories, avoid wildcards.
- **Returns**
    - Tuple of two elements, list of selcted directories and given path.



<h4 id="select_files" class="doc_header"><code>select_files</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L130" class="source_link" style="float:right">[source]</a></h4>

> <code>select_files</code>(**`path`**=*`'C:\\Users\\mass_\\Desktop'`*, **`include`**=*`[]`*, **`exclude`**=*`[]`*)

- Returns selected files from a given directory.
- **Parameters**
    - path    : path to a parent directory, default is `"."`
    - include : list of keywords to include files, avoid wildcards.
    - exclude : list of keywords to exclude files, avoid wildcards.
- **Returns**
    - Tuple of two elements, list of selcted files and given path.



<h4 id="get_child_items" class="doc_header"><code>get_child_items</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/g_utils.py#L154" class="source_link" style="float:right">[source]</a></h4>

> <code>get_child_items</code>(**`path`**=*`'C:\\Users\\mass_\\Desktop'`*, **`depth`**=*`None`*, **`recursive`**=*`True`*, **`include`**=*`[]`*, **`exclude`**=*`[]`*, **`filesOnly`**=*`False`*, **`dirsOnly`**=*`False`*)

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



<h4 id="plot_bands" class="doc_header"><code>plot_bands</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L8" class="source_link" style="float:right">[source]</a></h4>

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
    - path_evr   : path/to/vasprun.xml or output of `export_vasprun`. Auto picks in CWD.
    - ax         : Matplotlib axes object, if not given, one is created.
    - skipk      : Number of kpoints to skip, default will be from IBZKPT.
    - joinPathAt : Points where kpath is broken.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from `export_vasprun`.
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

> <code>add_legend</code>(**`ax`**=*`None`*, **`colors`**=*`[(1, 0, 0), (0, 1, 0), (0, 0, 1)]`*, **`labels`**=*`['s', 'p', 'd']`*, **`styles`**=*`'solid'`*, **`widths`**=*`0.7`*, **`anchor`**=*`(0, 1)`*, **`ncol`**=*`3`*, **`loc`**=*`'lower left'`*, **`fontsize`**=*`'small'`*, **`frameon`**=*`False`*, **\*\*`legend_kwargs`**)

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
    - kapath   : `export_vasprun`().kpath or `get_kpts`().kpath.
    - evals_set: `export_vasprun`().bands.evals or `get_evals`().evals. If calculations are spin-polarized, it will be `...evals.SpinUp/SpinDown` for both. You need to create collections twice for SpinUp and SpinDown separately.
    - pros_set : `export_vasprun().pro_bands.pros` or `get_bands_pro_set`().pros. If calculations are spin-polarized, it will be `...pros.SpinUp/SpinDown` for both. You need to create collections twice for SpinUp and SpinDown separately.
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
    - path_evr   : path/to/vasprun.xml or output of `export_vasprun`. Auto picks in CWD.
    - ax         : Matplotlib axes object, if not given, one is created.
    - skipk      : Number of kpoints to skip, default will be from IBZKPT.
    - joinPathAt : Points where kpath is broken.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from `export_vasprun`.
    - xt_indices : High symmetry kpoints indices.abs
    - xt_labels  : High Symmetry kpoints labels.
    - elements   : List [[0],[],[]] by default and plots s orbital of first ion..
    - orbs       : List [[r],[g],[b]] of indices of orbitals, could be empty, but shape should be same.
    - labels  : List [str,str,str] of projection labels. empty string should exist to maintain shape.
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



<h4 id="quick_color_lines" class="doc_header"><code>quick_color_lines</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L582" class="source_link" style="float:right">[source]</a></h4>

> <code>quick_color_lines</code>(**`path_evr`**=*`None`*, **`axes`**=*`None`*, **`skipk`**=*`None`*, **`joinPathAt`**=*`[]`*, **`elim`**=*`[]`*, **`elements`**=*`[[0]]`*, **`orbs`**=*`[[0]]`*, **`labels`**=*`['s']`*, **`color_map`**=*`'gist_rainbow'`*, **`max_width`**=*`2.5`*, **`xt_indices`**=*`[0, -1]`*, **`xt_labels`**=*`['$\\Gamma$', 'M']`*, **`E_Fermi`**=*`None`*, **`showlegend`**=*`True`*, **`figsize`**=*`(3.4, 2.6)`*, **`txt`**=*`None`*, **`xytxt`**=*`[0.05, 0.85]`*, **`ctxt`**=*`'black'`*, **`spin`**=*`'both'`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*, **`legend_kwargs`**=*`{'ncol': 4, 'anchor': (0, 0.85), 'handletextpad': 0.5, 'handlelength': 1, 'fontsize': 'small', 'frameon': True}`*, **\*\*`subplots_adjust_kwargs`**)

- Returns axes object and plot on which all matplotlib allowed actions could be performed. If given, axes,elements,orbs colors, and labels must have same length. If not given, zeroth ion is plotted with s-orbital.
- **Parameters**
    - path_evr   : Path/to/vasprun.xml or output of `export_vasprun`. Auto picks in CWD.
    - axes       : Matplotlib axes object with one or many axes, if not given, auto created.
    - skipk      : Number of kpoints to skip, default will be from IBZKPT.
    - joinPathAt : Points where kpath is broken.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from `export_vasprun`.
    - xt_indices : High symmetry kpoints indices.abs
    - xt_labels  : High Symmetry kpoints labels.
    - elements   : List [[0],], by defualt and plot first ion's projections.
    - orbs       : List [[0],] lists of indices of orbitals, could be empty.
    - labels     : List [str,] of orbitals labels. len(labels)==len(orbs) must hold.
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



<h4 id="init_figure" class="doc_header"><code>init_figure</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L770" class="source_link" style="float:right">[source]</a></h4>

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



<h4 id="select_pdos" class="doc_header"><code>select_pdos</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L834" class="source_link" style="float:right">[source]</a></h4>

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



<h4 id="collect_dos" class="doc_header"><code>collect_dos</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L874" class="source_link" style="float:right">[source]</a></h4>

> <code>collect_dos</code>(**`path_evr`**=*`None`*, **`elim`**=*`[]`*, **`elements`**=*`[[0]]`*, **`orbs`**=*`[[0]]`*, **`labels`**=*`['s']`*, **`E_Fermi`**=*`None`*, **`spin`**=*`'both'`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*)

- Returns lists of energy,tdos, pdos and labels. If given,elements,orbs and labels must have same length. If not given, zeroth ions is collected with s-orbital.
- **Parameters**)
    - path_evr   : Path/to/vasprun.xml or output of `export_vasprun`. Auto picks in CWD.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from `export_vasprun`.
    - elements   : List [[0],], by defualt and plot first ion's projections.
    - orbs       : List [[0],] lists of indices of orbitals, could be empty.
    - labels     : List [str,] of orbitals labels. len(labels)==len(orbs) must hold.
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



<h4 id="quick_dos_lines" class="doc_header"><code>quick_dos_lines</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L990" class="source_link" style="float:right">[source]</a></h4>

> <code>quick_dos_lines</code>(**`path_evr`**=*`None`*, **`ax`**=*`None`*, **`elim`**=*`[]`*, **`include_dos`**=*`'both'`*, **`elements`**=*`[[0]]`*, **`orbs`**=*`[[0]]`*, **`labels`**=*`['s']`*, **`color_map`**=*`'gist_rainbow'`*, **`tdos_color`**=*`(0.8, 0.95, 0.8)`*, **`linewidth`**=*`0.5`*, **`fill_area`**=*`True`*, **`vertical`**=*`False`*, **`E_Fermi`**=*`None`*, **`figsize`**=*`(3.4, 2.6)`*, **`txt`**=*`None`*, **`xytxt`**=*`[0.05, 0.85]`*, **`ctxt`**=*`'black'`*, **`spin`**=*`'both'`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*, **`showlegend`**=*`True`*, **`legend_kwargs`**=*`{'ncol': 4, 'anchor': (0, 1), 'handletextpad': 0.5, 'handlelength': 1, 'fontsize': 'small', 'frameon': True}`*)

- Returns ax object (if ax!=False) and plot on which all matplotlib allowed actions could be performed, returns lists of energy,tdos and pdos and labels. If given,elements,orbs colors, and labels must have same length. If not given, zeroth ions is plotted with s-orbital.
- **Parameters**)
    - path_evr   : Path/to/vasprun.xml or output of `export_vasprun`. Auto picks in CWD.
    - ax         : Matplotlib axes object, if None, one is created. If False, data lists are returned.
    - include_dos: One of {'both','tdos','pdos'}.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from `export_vasprun`.
    - elements   : List [[0],], by defualt and plot first ion's projections.
    - orbs       : List [[0],] lists of indices of orbitals, could be empty.
    - labels     : List [str,] of orbitals labels. len(labels)==len(orbs) must hold.
    - color_map  : Matplotlib's standard color maps. Default is 'gist_ranibow'.
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



<h4 id="plt_to_html" class="doc_header"><code>plt_to_html</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/s_plots.py#L1135" class="source_link" style="float:right">[source]</a></h4>

> <code>plt_to_html</code>(**`plt_fig`**=*`None`*, **`dpi`**=*`300`*, **`dash_html`**=*`None`*)

- Returns base64 encoded Image to display in notebook or HTML <img/> or plotly's dash_html_components.Img object.
- **Parameters**
    - plt_fig  : Matplotlib's figure instance, auto picks as well.
    - dpi      : PNG images's DPI, default is 300.
    - dash_html: Default is None which results in an image display in jupyter notebook.
        - If True, returns html.Img object for plotly's dash.
        - If False, returns html <img/> object to embed in HTML DOM.



<h4 id="get_rgb_data" class="doc_header"><code>get_rgb_data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L6" class="source_link" style="float:right">[source]</a></h4>

> <code>get_rgb_data</code>(**`kpath`**=*`None`*, **`evals_set`**=*`None`*, **`pros_set`**=*`None`*, **`elements`**=*`[[0], [], []]`*, **`orbs`**=*`[[0], [], []]`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*, **`scale_color`**=*`False`*)

- Returns a formatted RGB colored data to pass into `rgb_to_plotly` function. Two arguments, `elements` and `orbs` should be in one-to-one correspondence. Returned item has transpose data shape, so that main iteration is over bands.
- **Parameters**
    - ax       : Matplotlib axes object, if not given, linecollection is returned.
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



<h4 id="rgb_to_plotly" class="doc_header"><code>rgb_to_plotly</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L73" class="source_link" style="float:right">[source]</a></h4>

> <code>rgb_to_plotly</code>(**`rgb_data`**=*`None`*, **`mode`**=*`'markers'`*, **`max_width`**=*`5`*, **`showlegend`**=*`False`*, **`name`**=*`''`*, **`labels`**=*`['s', 'p', 'd']`*)

- Returns data object of plotly's figure using `get_rgb_data`. Returned data could be fed to a plolty's figure.
- ** Parameters**
    - rgb_data    : output of `get_rgb_data`.
    - mode        : Three plotting modes are available:
        - 'markers' : Plot whole data as a single scatter object. Its too fast.
        - 'bands'   : Plot data such that each band is accessible via legend.
        - 'lines'   : A replica of `matplotlib LineCollection` object. It plots at each point separately, slower than other two modes.
    - max_width  : Line/Scatter thickness is scaled to `max_width`.
    - name       : Name to be shown on hover text or legend.
    - labels     : Optional, show red green blue colors corresponding orbitals.
    - showlegend : Optional, only suitbale if spin up/down or 'bands' mode is ON.



<h4 id="plotly_to_html" class="doc_header"><code>plotly_to_html</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L134" class="source_link" style="float:right">[source]</a></h4>

> <code>plotly_to_html</code>(**`fig`**, **`filename`**=*`'new_plot.html'`*)

- Writes plotly's figure as HTML file which is accessible when online. If you need to have offline working file, just use `fig.write_html('file.html')` which will be larger in size.
- **Parameters**
    - fig      : A plotly's figure object.
    - filename : name of file to save fig.



<h4 id="plotly_rgb_lines" class="doc_header"><code>plotly_rgb_lines</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L163" class="source_link" style="float:right">[source]</a></h4>

> <code>plotly_rgb_lines</code>(**`path_evr`**=*`None`*, **`elements`**=*`[[], [], []]`*, **`orbs`**=*`[[], [], []]`*, **`labels`**=*`['', '', '']`*, **`mode`**=*`'markers'`*, **`elim`**=*`[]`*, **`E_Fermi`**=*`None`*, **`skipk`**=*`None`*, **`joinPathAt`**=*`[]`*, **`max_width`**=*`6`*, **`title`**=*`None`*, **`xt_indices`**=*`[0, -1]`*, **`xt_labels`**=*`['Γ', 'M']`*, **`figsize`**=*`None`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*)

- Returns plotly's figure object, takes care of spin-polarized calculations automatically. `elements`,`orbs` and `labels` are required to be one-to-one lists of size 3 where each item in list could be another list or integer.
- **Parameters**
    - path_ever  : Path/to/vasprun.xml or xml output of `read_asxml`.
    - elements   : List of size 3 of list of indices of ions. If not given, picks all ions for each orbital.
    - orbs       : List of size 3 of list of orbital indices, if not gievn, s,p,d plotted.
    - labels  : List of labels for projection.
    - mode       : Three plotting modes are available:
        - 'markers' : Plot whole data as a single scatter object. Its too fast.
        - 'bands'   : Plot data such that each band is accessible via legend.
        - 'lines'   : A replica of `matplotlib LineCollection` object. It plots at each point separately, slower than other two modes.
    - **kwargs      : interpolate, ticks, figsize,elim,joinPathAt,max_width,title etc.



<h4 id="plotly_dos_lines" class="doc_header"><code>plotly_dos_lines</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L317" class="source_link" style="float:right">[source]</a></h4>

> <code>plotly_dos_lines</code>(**`path_evr`**=*`None`*, **`elim`**=*`[]`*, **`elements`**=*`[[0]]`*, **`orbs`**=*`[[0]]`*, **`labels`**=*`['s']`*, **`color_map`**=*`'gist_rainbow'`*, **`tdos_color`**=*`(0.5, 0.95, 0)`*, **`linewidth`**=*`2`*, **`fill_area`**=*`True`*, **`vertical`**=*`False`*, **`E_Fermi`**=*`None`*, **`figsize`**=*`None`*, **`spin`**=*`'both'`*, **`interpolate`**=*`False`*, **`n`**=*`5`*, **`k`**=*`3`*, **`title`**=*`None`*)

- Returns ax object (if ax!=False) and plot on which all matplotlib allowed actions could be performed, returns lists of energy,tdos and pdos and labels. If given,elements,orbs colors, and labels must have same length. If not given, zeroth ions is plotted with s-orbital.
- **Parameters**)
    - path_evr   : Path/to/vasprun.xml or output of `export_vasprun`. Auto picks in CWD.
    - elim       : [min,max] of energy range.
    - E_Fermi    : If not given, automatically picked from `export_vasprun`.
    - elements   : List [[0,],] of ions indices, by defualt plot first ion's projections.
    - orbs       : List [[0,],] lists of indices of orbitals, could be empty.
    - labels     : List [str,] of orbitals labels. len(labels)==len(orbs) must hold.
    - color_map  : Matplotlib's standard color maps. Default is 'gist_ranibow'.
    - fill_area  : Default is True and plots filled area for dos. If False, plots lines only.
    - vertical   : False, If True, plots along y-axis.
    - figsize    : Tuple in pixels (width,height).
    - interpolate: Default is False, if True, bands are interpolated.
    - n          : int, number of points, default is 5.
    - k          : int, order of interpolation 0,1,2,3. Defualt 3. `n > k` should be hold.
    - legend_kwargs: Dictionary to contain legend arguments to fix.
- **Returns**
    - fig        : Plotly's figure object.



<h4 id="iplotfromtxt" class="doc_header"><code>iplotfromtxt</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/i_plots.py#L440" class="source_link" style="float:right">[source]</a></h4>

> <code>iplotfromtxt</code>(**`path_to_dir`**=*`'.'`*, **`ions`**=*`[0]`*, **`orbs`**=*`[[0], [], []]`*, **`labels`**=*`['s', '', '']`*, **`elim`**=*`[-5, 5]`*, **`tick_indices`**=*`[0, -1]`*, **`tick_labels`**=*`['Γ', 'M']`*, **`force_load`**=*`False`*)

- Returns plotly's figure object, Not applicable for spin-polarized calculations. It is just to keep as it was starting code. Remember where you did start!
- **Parameters**
    - path_to_dir: Path/to/directory where `Export-VaspRun`'s output from powershell is present.
    - ions       : List of indices of ions. If not given, picks all ions.
    - orbs       : List of size 3 of list of orbital indices, if not gievn, s,p,d plotted.
    - labels  : List of labels for projection.
    - elim       : Energy limit.
    - force_load : If True, it will run powershell to generate text files only if files are not already present.
    - **kwargs   : tick_indices, tick_labels, elim etc.



<h4 id="save_mp_API" class="doc_header"><code>save_mp_API</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L8" class="source_link" style="float:right">[source]</a></h4>

> <code>save_mp_API</code>(**`api_key`**)

- Save materials project api key for autoload in functions.



<h4 id="load_mp_data" class="doc_header"><code>load_mp_data</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L29" class="source_link" style="float:right">[source]</a></h4>

> <code>load_mp_data</code>(**`formula`**, **`api_key`**=*`None`*, **`mp_id`**=*`None`*, **`max_sites`**=*`None`*)

- Returns fetched data using request api of python form materials project website.
- **Parameters**
    - formula  : Material formula such as 'NaCl'.
    - api_key  : API key for your account from material project site. Auto picks if you already used `save_mp_API` function.
    - mp_id    : Optional, you can specify material ID to filter results.
    -max_sites : Option, you can set maximum number of sites to load fastly as it will not fetch all large data sets.



<h4 id="get_crystal" class="doc_header"><code>get_crystal</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L72" class="source_link" style="float:right">[source]</a></h4>

> <code>get_crystal</code>(**`formula`**, **`api_key`**=*`None`*, **`mp_id`**=*`None`*, **`max_sites`**=*`None`*)

- Returns crystal information dictionary including cif data format.
- **Parameters**
    - formula  : Material formula such as 'NaCl'.
    - api_key  : API key for your account from material project site. Auto picks if you already used `save_mp_API` function.
    - mp_id    : Optional, you can specify material ID to filter results.
    -max_sites : Option, you can set maximum number of sites to load fastly as it will not fetch all large data sets.



<h4 id="get_poscar" class="doc_header"><code>get_poscar</code><a href="https://github.com/massgh/pivotpy/tree/master/pivotpy/sio.py#L96" class="source_link" style="float:right">[source]</a></h4>

> <code>get_poscar</code>(**`formula`**, **`api_key`**=*`None`*, **`mp_id`**=*`None`*, **`max_sites`**=*`None`*)

- Returns poscar information dictionary including cif data format.
- **Parameters**
    - formula  : Material formula such as 'NaCl'.
    - api_key  : API key for your account from material project site. Auto picks if you already used `save_mp_API` function.
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



<h4 id="show" class="doc_header"><code>show</code><a href="matplotlib/pyplot.py#L251" class="source_link" style="float:right">[source]</a></h4>

> <code>show</code>(**\*`args`**, **\*\*`kw`**)

Display a figure.

When running in ipython with its pylab mode, display all
figures and return to the ipython prompt.

In non-interactive mode, display all figures and block until
the figures have been closed; in interactive mode it has no
effect unless figures were created prior to a change from
non-interactive to interactive mode (not recommended).  In
that case it displays the figures but does not block.

A single experimental keyword argument, *block*, may be
set to True or False to override the blocking behavior
described above.



<h4 id="savefig" class="doc_header"><code>savefig</code><a href="matplotlib/pyplot.py#L719" class="source_link" style="float:right">[source]</a></h4>

> <code>savefig</code>(**\*`args`**, **\*\*`kwargs`**)

Save the current figure.

Call signature::

  savefig(fname, dpi=None, facecolor='w', edgecolor='w',
          orientation='portrait', papertype=None, format=None,
          transparent=False, bbox_inches=None, pad_inches=0.1,
          frameon=None, metadata=None)

The output formats available depend on the backend being used.

Parameters
----------

fname : str or PathLike or file-like object
    A path, or a Python file-like object, or
    possibly some backend-dependent object such as
    `matplotlib.backends.backend_pdf.PdfPages`.

    If *format* is not set, then the output format is inferred from
    the extension of *fname*, if any, and from :rc:`savefig.format`
    otherwise.  If *format* is set, it determines the output format.

    Hence, if *fname* is not a path or has no extension, remember to
    specify *format* to ensure that the correct backend is used.

Other Parameters
----------------

dpi : [ *None* | scalar > 0 | 'figure' ]
    The resolution in dots per inch.  If *None*, defaults to
    :rc:`savefig.dpi`.  If 'figure', uses the figure's dpi value.

quality : [ *None* | 1 <= scalar <= 100 ]
    The image quality, on a scale from 1 (worst) to 95 (best).
    Applicable only if *format* is jpg or jpeg, ignored otherwise.
    If *None*, defaults to :rc:`savefig.jpeg_quality` (95 by default).
    Values above 95 should be avoided; 100 completely disables the
    JPEG quantization stage.

optimize : bool
    If *True*, indicates that the JPEG encoder should make an extra
    pass over the image in order to select optimal encoder settings.
    Applicable only if *format* is jpg or jpeg, ignored otherwise.
    Is *False* by default.

progressive : bool
    If *True*, indicates that this image should be stored as a
    progressive JPEG file. Applicable only if *format* is jpg or
    jpeg, ignored otherwise. Is *False* by default.

facecolor : color spec or None, optional
    The facecolor of the figure; if *None*, defaults to
    :rc:`savefig.facecolor`.

edgecolor : color spec or None, optional
    The edgecolor of the figure; if *None*, defaults to
    :rc:`savefig.edgecolor`

orientation : {'landscape', 'portrait'}
    Currently only supported by the postscript backend.

papertype : str
    One of 'letter', 'legal', 'executive', 'ledger', 'a0' through
    'a10', 'b0' through 'b10'. Only supported for postscript
    output.

format : str
    The file format, e.g. 'png', 'pdf', 'svg', ... The behavior when
    this is unset is documented under *fname*.

transparent : bool
    If *True*, the axes patches will all be transparent; the
    figure patch will also be transparent unless facecolor
    and/or edgecolor are specified via kwargs.
    This is useful, for example, for displaying
    a plot on top of a colored background on a web page.  The
    transparency of these patches will be restored to their
    original values upon exit of this function.

bbox_inches : str or `~matplotlib.transforms.Bbox`, optional
    Bbox in inches. Only the given portion of the figure is
    saved. If 'tight', try to figure out the tight bbox of
    the figure. If None, use savefig.bbox

pad_inches : scalar, optional
    Amount of padding around the figure when bbox_inches is
    'tight'. If None, use savefig.pad_inches

bbox_extra_artists : list of `~matplotlib.artist.Artist`, optional
    A list of extra artists that will be considered when the
    tight bbox is calculated.

metadata : dict, optional
    Key/value pairs to store in the image metadata. The supported keys
    and defaults depend on the image format and backend:

    - 'png' with Agg backend: See the parameter ``metadata`` of
      `~.FigureCanvasAgg.print_png`.
    - 'pdf' with pdf backend: See the parameter ``metadata`` of
      `~.backend_pdf.PdfPages`.
    - 'eps' and 'ps' with PS backend: Only 'Creator' is supported.

pil_kwargs : dict, optional
    Additional keyword arguments that are passed to `PIL.Image.save`
    when saving the figure.  Only applicable for formats that are saved
    using Pillow, i.e. JPEG, TIFF, and (if the keyword is set to a
    non-None value) PNG.



```python

```
