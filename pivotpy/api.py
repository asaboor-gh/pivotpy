# AUTOGENERATED! DO NOT EDIT! File to edit: MainAPI.ipynb (unless otherwise specified).

__all__ = ['download_structure', '__all__', 'parse_text', 'POSCAR', 'LOCPOT', 'get_axes', 'Vasprun']

# Cell
import os
import numpy as np
import plotly.graph_objects as go
try:
    from pivotpy import vr_parser as vp
    from pivotpy import splots as sp
    from pivotpy import iplots as ip
    from pivotpy import sio as sio
    from pivotpy import widgets as wdg
    from pivotpy import utils as gu
    from pivotpy import serializer
    from pivotpy import surfaces as srf
except:
    import pivotpy.vr_parser as vp
    import pivotpy.splots as sp
    import pivotpy.iplots as ip
    import pivotpy.sio as sio
    import pivotpy.widgets as wdg
    import pivotpy.utils as gu
    import pivotpy.serializer as serializer
    import pivotpy.surfaces as srf

def _sub_doc(from_func,skip_param=None,replace={}):
    """Assing __doc__ from other function. Replace words in docs where need."""
    def wrapper(func):
        docs = '\n'.join(line for line in from_func.__doc__.splitlines() if skip_param not in line)
        for k,v in replace.items():
            docs = docs.replace(k,v)
        func.__doc__ = docs
        return func
    return wrapper

# Cell
def download_structure(formula, mp_id=None, max_sites=None,min_sites=None, api_key=None,save_key = False):
    """Download structure data from Materials project website.
    - **Parameters**
        - formula: chemical formula of the material.
        - mp_id: Materials project id of material.
        - max_sites: maximum number of sites in structure to download.
        - min_sites: minimum number of sites in structure to download.
    > max_sites and min_sites are used to filter the number of sites in structure, or use mp_id to download a specific structure.
    - **One Time API Key**
        - api_key: API key from Materials project websit, if you use save_key=True, never required again.
        - save_key: Save API key to file. You can save any time of key or device changed.
    - **Return**
        List of Structure data containing attribute/method `cif`/`export_poscar, write_cif` etc.
    """
    mp = sio.InvokeMaterialsProject(api_key= api_key)
    output = mp.request(formula=formula,mp_id=mp_id,max_sites=max_sites,min_sites=min_sites) # make a request
    if save_key and isinstance(api_key,str):
        mp.save_api_key(api_key)
    if mp.success:
        return output
    else:
        raise ConnectionError('Connection was not sccessful. Try again!')

# Cell
# Direct function exports from modules
_memebers = (
    srf.SpinDataFrame,
    sio.get_kpath,
    sio.str2kpath,
    sio.fancy_quiver3d,
    sio.rotation,
    wdg.generate_summary,
    wdg.VasprunApp,
    wdg.KPathApp,
    gu.set_dir,
    gu.get_child_items,
    gu.transform_color,
    gu.interpolate_data,
    vp.split_vasprun,
    vp.xml2dict,
    ip.iplot2html,
    sp.plt2html,
    sp.plt2text,
    sp.show,
    sp.savefig,
    sp.append_axes,
    sp.join_axes
)

# Subset of functions from modules in __all__ to make exportable as *
__all__ = [*[_m.__name__ for _m in _memebers],*[a for a in __all__ if a != '__all__']]
for _m in _memebers:
    locals()[_m.__name__] = _m # Assign to local namespace that can be exported, classes only have __name__, not name

# Cell
def parse_text(path,
          dtype=float, delimiter='\s+',
          include=None, exclude='#',
          raw=False, fix_format=True,
          start=0, nlines=None, count=-1,
          new_shape=None,cols=None,
          old_shape=None, slice_rows=None #Works if both given
          ):
    """
    - Reads a sliced array from txt,csv type files and return to array.
        Also manages if columns lengths are not equal and return 1D array.
        It is faster than loading  whole file into memory. This single function could be used
        to parse EIGENVAL, PROCAR, DOCAR and similar files with just a
        combination of `exclude, include,start,stop,step` arguments.
    - **Parameters**
        - path: Path/to/file to be parsed.
        - dtype: float by default. Data type of output array, it is must have argument.
        - start,nlines: The indices of lines to start reading from and number of lines after start respectively.
            Both could be None or int, while start could be a list to read slices from file provided that nlines is int.
            The spacing between adjacent indices in start should be equal to or greater than nlines as pointer in file
            do not go back on its own.  These parameters are not required if `slice_rows` and `old_shape` are given.
            > Note: `start` should count comments if `exclude` is None.
        - count: `np.size(output_array) = nrows x ncols`, if it is known before execution, performance is increased.
        - delimiter:  Default is `\s+`. Could be any kind of delimiter valid in numpy and in the file.
        - cols: List of indices of columns to pick. Useful when reading a file like PROCAR which e.g. has text and numbers inline. This parameter is in output of `slice_data`.
        - include: Default is None and includes everything. String of patterns separated by | to keep, could be a regular expression.
        - exclude: Default is '#' to remove comments. String of patterns separated by | to drop,could be a regular expression.
        - raw    : Default is False, if True, returns list of raw strings. Useful to select `cols`.
        - fix_format: Default is True, it sepearates numbers with poor formatting like 1.000-2.000 to 1.000 2.000 which is useful in PROCAR. Keep it False if want to read string literally.
        - new_shape : Tuple of shape Default is None. Will try to reshape in this shape, if fails fallbacks to 2D or 1D. Not required if `slice_rows` and `old_shape` are given.
        - old_shape: It is required when you need to pick blocks of data from rows.
                columns should be last entry, like (2,2,3) means 3 columns and two different indexed blocks are there.
                Only works if `slice_rows` is given too.
        - slice_rows: It is required when you need to pick blocks of data from rows.
                [(0,1),(0,1)] will pick lines at index 0,1,3,4 if first dimension has size 3. It is like N1*i+j for N1=3.
                General formula to pick rows is `inner_block_index + inner_block_size*second_block_index + inner_most_size*second_block_size*third_block_index + ...`
                `i_1 + N_1*i_2 + N_1*N_2*i_3 + ...` where i_1 is inner most index.
                Only works if `old_shape` is given too.
    - **Examples**
        > `parse_text('path/to/PROCAR',start=3,include='k-point',cols=[3,4,5])[:2]`
        > array([[ 0.125,  0.125,  0.125],
        >        [ 0.375,  0.125,  0.125]])
        > `parse_text('path/to/EIGENVAL',start=7,exclude='E',cols=[1,2])[:2]`
        > array([[-11.476913,   1.      ],
        >        [  0.283532,   1.      ]])
    > Note: Slicing a dimension to 100% of its data is faster than let say 80% for inner dimensions, so if you have to slice more than 50% of an inner dimension, then just load full data and slice after it.

    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"File {path!r} does not exists")

    extra_kws = dict(dtype=dtype,delimiter=delimiter, # Data related
                     include=include, exclude= exclude, # selection related
                     raw=raw, fix_format= fix_format,new_shape=new_shape, # Output related
                     start=start,nlines=nlines,count=count,cols=cols # slicing related
                     )
    if old_shape and slice_rows:
        _cols = -1 if cols == None else cols
        extra_kws.update(vp.slice_data(dim_inds=[*slice_rows,_cols],old_shape=old_shape))
        print(extra_kws)

    return vp.islice2array(path_or_islice=path,**extra_kws)

# Cell
from contextlib import suppress
class POSCAR:
    "POSACR class to contain data and related methods, data is PoscarData, json/tuple file/string."
    def __init__(self,path = None,text_plain = None,data = None):
        """Do not use `data` yourself, it's for operations on poscar.
        Prefrence order: data, text_plain, path"""
        self.path = path
        self.text_plain = text_plain
        self._kpts_info = None # Get defualt kpts_info
        if data:
            self._data = serializer.PoscarData.validated(data)
        else:
            self._data = sio.export_poscar(path=path,text_plain=text_plain)
            with suppress(BaseException): # Only reuqired here,not in vasprun_data or spin_data
                base_dir = os.path.split(os.path.abspath(path or './POSCAR'))[0]
                self._kpts_info = vp.get_kpoints_info(os.path.join(base_dir,'KPOINTS'))
        # These after data to work with data
        self.primitive = False
        self._bz = self.get_bz(primitive = False) # Get defualt regular BZ from sio
        self._cell = self.get_cell() # Get defualt cell
        self._plane = None # Get defualt plane, changed with splot_bz
        self._ax = None # Get defualt axis, changed with splot_bz

    def get_kpoints_info(self, other_path):
        "Return kpoints info from other_path to be POSCAR path, required for kpoints to bring in bz."
        self._kpts_info = vp.get_kpoints_info(other_path)
        return self._kpts_info

    @property
    def data(self):
        "Data object in POSCAR."
        return self._data

    @property
    def bz(self):
        return self._bz

    @property
    def cell(self):
        return self._cell

    @_sub_doc(sio.get_bz,'- path_pos')
    def get_bz(self, loop=True, digits=8, primitive=False):
        self._bz = sio.get_bz(path_pos = self._data.basis, loop=loop, digits=digits, primitive=primitive)
        self.primitive = primitive
        return self._bz

    def set_bz(self,primitive=False,loop=True,digits=8):
        """Set BZ in primitive or regular shape. returns None, just set self.bz"""
        self.get_bz(primitive=primitive,loop=loop,digits=digits)

    def get_cell(self, loop=True, digits=8):
        "See docs of `get_bz`, same except space is inverted."
        self._cell = sio.get_bz(path_pos=self._data.rec_basis,loop=loop, digits=digits, primitive=True) # cell must be primitive
        return self._cell

    @_sub_doc(sio.splot_bz,'- path_pos_bz')
    def splot_bz(self, ax=None, plane=None, color='blue', fill=True, vectors=True, v3=False, vname='b', colormap='plasma', light_from=(1, 1, 1), alpha=0.4):
        self._plane = plane # Set plane for splot_kpath
        new_ax = sio.splot_bz(path_pos_bz = self._bz, ax=ax, plane=plane, color=color, fill=fill, vectors=vectors, v3=v3, vname=vname, colormap=colormap, light_from=light_from, alpha=alpha)
        self._ax = new_ax # Set ax for splot_kpath
        return new_ax

    def splot_kpath(self,vertex = 0, knn_inds = None, labels = None, color='k', line_width = 0.8,marker_size = 10,marker_style = '.',**labels_kwargs):
        """Plot k-path over existing BZ.
        - **Parameters**
            - vertex: vertex index nearby which to plot path. There are as many vertices as there are in BZ's shape.
            - knn_inds: list of indices of k nearest points e.g. [2,3,1] will trace path linking as 2-3-1.
                0 is Gamma point and 1 is the selected vertex itself. Points are taken internally from BZ, you can see from `self.bz.specials`.
            - labels: list of labels for each k-point in same order as `knn_inds`.
            - color, line_width, marker_size, marker_style are passed to `plt.plot`.

        labels_kwargs are passed to `plt.text`.

        > Tip: You can use this function multiple times to plot multiple/broken paths over same BZ.
        """
        if not self._bz or not self._ax:
            raise ValueError("BZ not found, use `splot_bz` first")

        _specials = self._bz.specials
        nearest = knn_inds

        ijk = [0,1,2]
        _mapping = {'xy':[0,1],'xz':[0,2],'yz':[1,2],'zx':[2,0],'zy':[2,1],'yx':[1,0]}
        if isinstance(self._plane, str) and self._plane in _mapping:
            ijk = _mapping[self._plane]

        inds = _specials.near[vertex]
        if nearest:
            inds = [inds[n] for n in nearest]

        if not labels:
            labels = ["[{0:6.3f}, {1:6.3f}, {2:6.3f}]".format(*_specials.kpoints[i]) for i in inds]
            if nearest:
                labels = [f"{n}: {_lab}" for n, _lab in zip(nearest, labels)]

        coords = _specials.coords[inds][:,ijk]
        self._ax.plot(*coords.T,color = color,linewidth=line_width,marker=marker_style,markersize=marker_size)

        for c,text in zip(coords, labels):
            self._ax.text(*c,text,**labels_kwargs)
        return self._ax


    def splot_cell(self, ax=None, plane=None, color='blue', fill=True, vectors=True, v3=False, vname='a', colormap='plasma', light_from=(1, 1, 1), alpha=0.4):
        "See docs of `splot_bz`, everything is same except space is inverted."
        return sio.splot_bz(path_pos_bz = self._cell, ax=ax, plane=plane, color=color, fill=fill, vectors=vectors, v3=v3, vname=vname, colormap=colormap, light_from=light_from, alpha=alpha)

    @_sub_doc(sio.iplot_bz,'- path_pos_bz')
    def iplot_bz(self, fill=True, color='rgba(168,204,216,0.4)', background='rgb(255,255,255)', vname='b', alpha=0.4, ortho3d=True, fig=None):
        return sio.iplot_bz(path_pos_bz = self._bz, fill=fill, color=color, background=background, vname=vname, alpha=alpha, ortho3d=ortho3d, fig=fig)

    def iplot_cell(self, fill=True, color='rgba(168,204,216,0.4)', background='rgb(255,255,255)', vname='a', alpha=0.4, ortho3d=True, fig=None):
        "See docs of `iplot_bz`, everything is same except space is iverted."
        return sio.iplot_bz(path_pos_bz = self._cell, fill=fill, color=color, background=background, vname=vname, alpha=alpha, ortho3d=ortho3d, fig=fig)

    @_sub_doc(sio.splot_lat,'- poscar')
    def splot_lat(self, sizes=50, colors=[], colormap=None, bond_length=None, tol=0.1, eps=0.01, eqv_sites=True, translate=None, line_width=1, edge_color=(1, 0.5, 0, 0.4), vectors=True, v3=False, plane=None, light_from=(1, 1, 1), fill=False, alpha=0.4, ax=None):
        return sio.splot_lat(poscar=self._data, sizes=sizes, colors=colors, colormap=colormap, bond_length=bond_length, tol=tol, eps=eps, eqv_sites=eqv_sites, translate=translate, line_width=line_width, edge_color=edge_color, vectors=vectors, v3=v3, plane=plane, light_from=light_from, fill=fill, alpha=alpha, ax=ax)

    @_sub_doc(sio.iplot_lat,'- poscar')
    def iplot_lat(self, sizes=10, colors='blue', bond_length=None, tol=0.1, eps=0.01, eqv_sites=True, translate=None, line_width=4, edge_color='black', fill=False, alpha=0.4, ortho3d=True, fig=None):
        return sio.iplot_lat(poscar=self._data, sizes=sizes, colors=colors, bond_length=bond_length, tol=tol, eps=eps, eqv_sites=eqv_sites, translate=translate, line_width=line_width, edge_color=edge_color, fill=fill, alpha=alpha, ortho3d=ortho3d, fig=fig)

    @_sub_doc(sio.write_poscar,'- poscar')
    def write(self, sd_list=None, outfile=None, overwrite=False):
        return sio.write_poscar(poscar = self._data, sd_list=sd_list, outfile=outfile, overwrite=overwrite)

    @_sub_doc(sio.join_poscars,'- poscar1',replace={'poscar2':'other'})
    def join(self,other, direction='z', tol=0.01):
        return self.__class__(data = sio.join_poscars(poscar1=self._data, poscar2=other.data, direction=direction, tol=tol))

    @_sub_doc(sio.scale_poscar,'- path_poscar')
    def scale(self, scale=(1, 1, 1), tol=0.01):
        return self.__class__(data = sio.scale_poscar(path_poscar=self._data, scale=scale, tol=tol))

    @_sub_doc(sio.rotate_poscar,'- path_poscar')
    def rotate(self,angle_deg,axis_vec):
        return self.__class__(data = sio.rotate_poscar(path_poscar=self._data, angle_deg = angle_deg, axis_vec=axis_vec))

    @_sub_doc(sio.fix_sites,'- poscar')
    def fix_sites(self, tol=0.01, eqv_sites=True, translate=None):
        return self.__class__(data = sio.fix_sites(poscar=self._data, tol=tol, eqv_sites=eqv_sites, translate=translate))

    def get_kmesh(self, *args, shift = 0, weight=None, cartesian = False, ibzkpt=None, outfile=None):
        """Generates uniform mesh of kpoints. Options are write to file, or return KPOINTS list. Use self.write(...) as well if you use this function to have KPOINTS and POSCAR similar way.
        - Positional arguments are 1 or 3 integers which decide shape of mesh. If 1, mesh points are equally spaced using information from POSCAR.
        - **Parameters**
            - shift  : Defualt is 0 and grid is created in interval [0,1], if given, grid is shifted to [shift,1+shift] in all directions.
            - weight : Float, if None, auto generates weights.
            - cartesian: If True, generates cartesian mesh and also reguires scale to be given.
            - ibzkpt : Path to ibzkpt file, required for HSE calculations.
            - outfile: Path/to/file to write kpoints.

        If `outfile = None`, KPOINTS file content is printed."""
        scale = 1 # All POSCARS are scaled to 1.0 if written by this class
        abc_norms = np.linalg.norm(self._data.basis, axis=1).round(4)
        return sio.get_kmesh(*args, shift = shift, weight = weight, abc_norms = abc_norms, cartesian = cartesian, scale = scale,ibzkpt= ibzkpt, outfile=outfile)

    def bring_in_cell(self,points, scale = None):
        """Brings atoms's positions inside Cell and returns their R3 coordinates.
        If points are cartesian, they are just scaled to fit inside the cell.
        The scale factor is usaully `a` which is on second line of POSCAR.
        """
        if self.data.cartesian:
            if isinstance(scale,(int,float)):
                return np.array(points) * scale
            else:
                raise RuntimeError('Found Cartesian POSCAR, tweak `scale = a`  where a is on second line of POSCAR.')
        return sio.to_R3(self._data.basis, points= points)

    @_sub_doc(sio.kpoints2bz,'- bz')
    def bring_in_bz(self,kpoints, scale = None):
        """Brings kpoints inside already set BZ, (primitive or regular).
        If kpoints are cartesian, they are just scaled to fit inside the BZ.
        The scale factor is usaully `2π/a` where a is on second line of POSCAR.
        """
        if not self._kpts_info:
            raise RuntimeError('Run `POSCAR.get_kpoints_info(other_path)` first. Only required once!')

        if self._kpts_info.cartesian:
            if isinstance(scale,(int,float)):
                return np.array(kpoints) * scale
            else:
                raise RuntimeError('Found Cartesian KPOINTS, tweak `scale ~ 2π/a`, where a is on second line of POSCAR.')

        if not self._bz:
            raise RuntimeError('No BZ found. Please run `get_bz()` first.')
        return sio.kpoints2bz(self._bz, kpoints= kpoints,primitive = self.primitive)

# Cell
class LOCPOT:
    def __init__(self,path=None,e = True,m = False):
        """
        - Returns Data from LOCPOT and similar structure files like CHG. Loads only single set out of 2/4 magnetization data to avoid performance/memory cost while can load electrostatic and one set of magnetization together.
        - **Parameters**
            - path: path/to/LOCPOT or similar stuructured file like CHG. LOCPOT is auto picked in CWD.
            - e   : Electric potential/charge density. Default is True.
            - m   : Magnetization density m. Default is False. If True, picks `m` for spin polarized case, and `m_x` for non-colinear case. Additionally it can take 'x','y' and 'z' in case of non-colinear calculations.
        - **Exceptions**
            - Would raise index error if magnetization density set is not present in LOCPOT/CHG in case `m` is not False.
        """
        self.path = path # Must be
        self.m = m # Required to put in plots.
        self._data = gu.export_potential(locpot=path, e=e,m=m)

    @property
    def data(self):
        return self._data

    @_sub_doc(sp.plot_potential,'- e_or_m')
    @_sub_doc(sp.plot_potential,'- basis')
    def splot_e(self,operation='mean_z',ax=None,period=None,
                 lr_pos=(0.25,0.75),lr_widths = [0.5,0.5],
                 labels=(r'$V(z)$',r'$\langle V \rangle _{roll}(z)$',r'$\langle V \rangle $'),
                 colors = ((0,0.2,0.7),'b','r'),annotate=True):
        return sp.plot_potential(basis=self._data.basis,e_or_m=self._data.e,operation=operation,
                                    ax=ax,period=period,lr_pos=lr_pos,lr_widths=lr_widths,
                                    labels=labels,colors=colors,annotate=annotate)

    @_sub_doc(sp.plot_potential,'- e_or_m')
    @_sub_doc(sp.plot_potential,'- basis')
    def splot_m(self,operation='mean_z',ax=None,period=None,
                lr_pos = (0.25,0.75),lr_widths = [0.5,0.5],
                labels = (r'$M(z)$',r'$\langle M \rangle _{roll}(z)$',r'$\langle M \rangle $'),
                colors = ((0,0.2,0.7),'b','r'),annotate=True):
        if self.m:
            try:
                e_or_m = self._data.m
            except:
                e_or_m = self._data.to_dict()[f'm_{self.m}']
        else:
            raise ValueError("Magnetization data set does not exist in {}".format(self.path))
        return sp.plot_potential(basis=self._data.basis,e_or_m=e_or_m,operation=operation,
                                    ax=ax,period=period,lr_pos=lr_pos,lr_widths=lr_widths,
                                    labels=labels,colors=colors,annotate=annotate)

    def view_period(self,period_guess=0.25,operation='mean_z',nslice=10,e_or_m=None,):
        """
        - Periodicity check by plotly's interactive plot.
        - **Parameters**
            - period_guess: Initial guess of period. Default is 0.25. Should be in [0,1].
            - operation   : Any of ['mean_x','min_x','max_x','mean_y','min_y','max_y','mean_z','min_z','max_z'].
            - nslice      : Default is 10. Number of periods around and including period_guess. e.g. If you give 0.25 as period_guess and nslice is 10, you will get 10 lines of rolling average over given data from where you can choose best fit or try another guess and so on.
            - e_or_m      : None by default. Not required in most cases as `view_period()` will try to get data itself from top class in order of `self._data.[e,m,m_x,m_y,m_z]` and if `self._data.e` exists it never goes to others, so you can overwrite this by setting `e_or_m = self._data.[your choice]`.
        """
        pos = period_guess
        check = ['mean_x','min_x','max_x','mean_y','min_y','max_y','mean_z','min_z','max_z']
        if operation not in check:
            raise ValueError("operation expects any of {!r}, got {}".format(check,operation))
        if e_or_m is None:
            try:
                data = self._data.e
            except:
                try:
                    data = self._data.m
                except:
                   data = self._data.to_dict()[f'm_{self.m}']
                else:
                    raise ValueError("Magnetization data set does not exist in {}".format(self.path))
        else:
            data = e_or_m

        _opr,_dir = operation.split('_')
        x_ind = 'xyz'.index(_dir)
        other_inds = tuple([i for i in [0,1,2] if i != x_ind])
        _func_ = np.min if _opr == 'min' else np.max if _opr == 'max' else np.mean

        fig = go.Figure()
        _arr = _func_(data,axis = other_inds)
        N = np.rint(pos*len(_arr)).astype(int)
        _range = range(int(N-nslice/2),int(N+nslice/2+1)) # +1 for range.
        for div in _range:
            if div > 0 and div < len(_arr):
                y = np.convolve(_arr+div,np.ones((div,))/div,mode='valid')
                x = np.linspace(0,1,len(y))
                h_text = ["{}: {:>5.3f}</br>v: {:>5.3f}".format(_dir,_h,_v-div) for _h,_v in zip(x,y)]
                fig.add_trace(go.Scatter(x=x,y=y,name="Roll_av({:>5.3f})".format(div/len(_arr)),hovertext=h_text))
        fig.update_layout(title = self._data.SYSTEM,font=dict(family="stix serif",size=14),
                          yaxis = go.layout.YAxis(title_text='No. of Points in Rolling Average'),
                          xaxis = go.layout.XAxis(title_text="{}({}<sub>max</sub>)".format(_dir,_dir)))
        return fig

# Cell
@_sub_doc(sp.get_axes,'- self',replace={'get_axes':'get_axes'})
def get_axes(figsize=(3.4, 2.6), nrows=1, ncols=1, widths=[], heights=[], axes_off=[], axes_3d=[], sharex=False, sharey=False, azim=45, elev=15, ortho3d=True, **subplots_adjust_kwargs):
    axes = sp.get_axes(figsize=figsize, nrows=nrows, ncols=ncols, widths=widths, heights=heights, axes_off=axes_off, axes_3d=axes_3d, sharex=sharex, sharey=sharey, azim=azim, elev=elev, ortho3d=ortho3d, **subplots_adjust_kwargs)
    for ax in np.array([axes]).flatten():
        for f in [sp.add_text,sp.add_legend,sp.add_colorbar,sp.color_wheel,sp.break_spines,sp.modify_axes,sp.append_axes, sp.join_axes]:
            if ax.name != '3d':
                setattr(ax,f.__name__,f.__get__(ax,type(ax)))
    return axes
get_axes.__doc__ = get_axes.__doc__ + '''
**There are extra methods added to each axes (only 2D) object.**
- add_text
- add_legend
- add_colorbar
- color_wheel
- break_spines
- modify_axes
- append_axes
- join_axes
'''

# Cell
class Vasprun:
    """
    - All plotting functions that depend on `export_vasprun` are joined under this class and renamed.

    - **Main Parameter**
        - path: str: path/to/vasprun.xml. Auto picks in CWD.

    - **Optional Parameters** (only useful if path is `vasprun.xml` file)
        - skipk      : int: Skip initial kpoints.
        - elim       : list: Energy range e.g. [-5,5].
        - shift_kpath: float: Shift in kpath values for side by side plotting.
        - try_pwsh   : bool: True by default, tries to load data exported using Powershell's `Vasp2Visual.Export-Vasprun` command.
        - data   : json/pickle file/str or VasprunData or a valid dictionary. Takes precedence over path parameter.

    - **Attributes and Methods**
        - data        : Exported data from given file. This has it's own attributes as well to save as json/pickle etc.
        - to_json     : Saves data in `.json` file. Useful for transport to other langauges.
        - to_pickle   : Saves data in `.pickle` file. Useful for fast reload in python.
        - splot_[...] : Plots data using `sp.splot_[...]` functions.
        - iplot_[...] : Plots data using `sp.iplot_[...]` functions.

    > Tip: If KPOINTS file is generated by this module, ticks on kpath are auto-picked.
    """
    def __init__(self,path = None,skipk = None,elim=[],shift_kpath=0,try_pwsh=True,data=None):
        if data: #json/pickle data strings
            self._data = serializer.VasprunData.validated(data)
        else:
            self._data = vp.export_vasprun(path=path,skipk=skipk,elim=elim,shift_kpath=shift_kpath,try_pwsh=try_pwsh)

        self.elim = elim
        self._kpath = self._data.kpath  # For info only, get updated with plot commands
        self._efermi = self._data.bands.E_Fermi   # For info only, get updated with plot commands

        if path == None:
            self.kticks = sio.read_ticks('KPOINTS')
        elif os.path.isfile(path):
            self.kticks = sio.read_ticks(os.path.join(os.path.dirname(path),'KPOINTS'))
        else:
            self.kticks = {} # no kticks available when loading from json/pickle data_str

    @property
    def poscar(self):
        """Returns POSCAR object that can be used for plotting BZ/Lattice etc.

        New in 1.1.5
        """
        return POSCAR(data = self._data.poscar)  #POSCAR class

    def __handle_kwargs(self,kwargs,dos=False):
        kwargs = {'elim': self.elim, **kwargs}
        if dos:
            return kwargs
        ticks = {k:self.kticks[k] for k in ['ktick_inds','ktick_vals','kseg_inds']}
        kwargs = {**ticks,**kwargs} #Prefer provided ones

        # Set for info only in case of bandstructure
        if  'kseg_inds' in kwargs and kwargs['kseg_inds']:
            self._kpath =  vp.join_ksegments(self._data.kpath,kwargs['kseg_inds'])
        if 'E_Fermi' in kwargs and kwargs['E_Fermi'] != None: #None is important to pick 0 as fermi as well
            self._efermi = kwargs['E_Fermi']
        return kwargs

    @_sub_doc(serializer.Dict2Data.to_json,'')
    def to_json(self,outfile=None,indent=1):
        return self._data.to_json(outfile=outfile,indent=indent)

    @_sub_doc(serializer.Dict2Data.to_pickle,'')
    def to_pickle(self,outfile=None):
        return self._data.to_pickle(outfile=outfile)

    @property
    def data(self):
        "Get exported data."
        return self._data

    def select(self,kpoints_inds = None, bands_inds = None, kseg_inds = None):
        """Seletc data based on kpoints and bands indices.
        This is useful to select only a subset of data and even reorder kpoints after calculations.
        Both   `kpoints_inds` and `bands_inds` are based on current data and should be based on zero indexing.
        `kseg_inds` is index of disconnected kpoints in `kpoints_inds`, e.g. in `kpoints_inds = [0,5,6,7]`, if 0 and 5 are disconnected, `kseg_inds = [1]`.

        **Returns** `Vasprun` object with selected data that can be plotted using `splot_[...]` or `iplot_[...]` functions.

        New in version 1.1.4
        """
        if kpoints_inds is None and bands_inds is None:
            return self

        assert isinstance(kpoints_inds,(list,tuple,range)) if kpoints_inds is not None else True
        assert isinstance(bands_inds,(list,tuple,range)) if bands_inds is not None else True

        d = self.data.to_dict()
        kpoints_inds = range(len(d['kpoints'])) if kpoints_inds is None else kpoints_inds
        bands_inds = range(len(d['bands']['indices'])) if bands_inds is None else bands_inds

        d['kpoints'] = d['kpoints'][kpoints_inds]
        d['kpath'] = [0, *np.linalg.norm(d['kpoints'][1:] - d['kpoints'][:-1],axis=1).cumsum().round(6)]
        d['kpath'] = vp.join_ksegments(d['kpath'],kseg_inds) # If broken kpath is provided, join it

        d['bands']['indices'] = tuple([d['bands']['indices'].start + b for b in bands_inds]) # It is range in original data

        if self.data.sys_info.ISPIN == 1:
            d['bands']['evals'] = d['bands']['evals'][kpoints_inds][:,bands_inds]
            d['pro_bands']['pros'] = d['pro_bands']['pros'][:,kpoints_inds][:,:,bands_inds,...]
        else:
            d['bands']['evals']['SpinUp'] = d['bands']['evals']['SpinUp'][kpoints_inds][:,bands_inds]
            d['bands']['evals']['SpinDown'] = d['bands']['evals']['SpinDown'][kpoints_inds][:,bands_inds]
            d['pro_bands']['pros']['SpinUp'] = d['pro_bands']['pros']['SpinUp'][:,kpoints_inds][:,:,bands_inds,...]
            d['pro_bands']['pros']['SpinDown'] = d['pro_bands']['pros']['SpinDown'][:,kpoints_inds][:,:,bands_inds,...]

        return self.__class__(data_str = serializer.VasprunData(d).to_json())


    @_sub_doc(sp.splot_bands,'- path_evr')
    def splot_bands(self,ax = None,**kwargs):
        kwargs = self.__handle_kwargs(kwargs)
        return sp.splot_bands(self._data,ax = ax, **kwargs)

    @_sub_doc(sp.splot_dos_lines,'- path_evr')
    def splot_dos_lines(self,elements = [[0],], orbs = [[0],], labels = ['s',], ax = None, query_data= {}, **kwargs):
        kwargs = self.__handle_kwargs(kwargs,dos=True)
        return sp.splot_dos_lines(self._data,elements = elements, orbs = orbs, labels = labels, ax = ax, query_data = query_data,**kwargs)

    @_sub_doc(sp.splot_rgb_lines,'- path_evr')
    def splot_rgb_lines(self,elements = [[],[],[]], orbs = [[],[],[]], labels = ['','',''], ax = None, query_data= {}, **kwargs):
        kwargs = self.__handle_kwargs(kwargs)
        return sp.splot_rgb_lines(self._data,elements = elements, orbs = orbs, labels = labels, ax = ax, query_data = query_data,**kwargs)

    @_sub_doc(sp.splot_color_lines,'- path_evr')
    def splot_color_lines(self,elements = [[0],], orbs = [[0],], labels = ['s',],axes = None, query_data= {}, **kwargs):
        kwargs = self.__handle_kwargs(kwargs)
        return sp.splot_color_lines(self._data,elements = elements, orbs = orbs, labels = labels, axes = axes, query_data = query_data,**kwargs)

    @_sub_doc(ip.iplot_dos_lines,'- path_evr')
    def iplot_dos_lines(self,elements = [[0],], orbs = [[0],], labels = ['s',], query_data= {}, **kwargs):
        kwargs = self.__handle_kwargs(kwargs, dos=True)
        return ip.iplot_dos_lines(self._data,elements = elements, orbs = orbs, labels = labels, query_data = query_data,**kwargs)

    @_sub_doc(ip.iplot_rgb_lines,'- path_evr')
    def iplot_rgb_lines(self,elements = [[],[],[]], orbs = [[],[],[]], labels = ['','',''], query_data= {}, **kwargs):
        kwargs = self.__handle_kwargs(kwargs)
        return ip.iplot_rgb_lines(self._data,elements = elements, orbs = orbs, labels = labels, query_data = query_data,**kwargs)

    def get_band_info(self,b_i,k_i=None):
        """Get band information for given band index `b_i`. If `k_i` is given, returns info at that point
        Fermi energy is subtracted from all energies. When a plot commnad is called, the Fermi energy is updated if provided.
        """
        def at_minmax(_bands,_pros,func,k_i=None):
            _bands_ = _bands.flatten() - self._efermi # subtract fermi energy
            if isinstance(k_i,int):
                extrema = _bands_[k_i]
                k = float(self._kpath[k_i])
                kp = self._data.kpoints[k_i]
                pros = _pros[:,k_i,:].sum(axis=0).flatten()
            else:
                extrema = func(_bands_)
                where, = np.where(_bands_ == extrema) # unpack singelton
                k, kp = [float(self._kpath[w]) for w in where], self._data.kpoints[where]
                pros = _pros[:,where[0],:].sum(axis=0).flatten()
            return serializer.Dict2Data({'e':float(extrema),'k':k,'kp':kp.tolist(),
                    'pros':{l.replace('-',''):float(p) for p,l in zip(pros,self._data.pro_bands.labels)}})

        if self._data.bands.ISPIN == 1:
            b = self._data.bands.evals[:,b_i]
            p = self._data.pro_bands.pros[:,:,b_i,:]

            if isinstance(k_i,int): # single kpoint
                return at_minmax(b,p,np.min,k_i=k_i)

            return serializer.Dict2Data({'min':at_minmax(b,p,np.min,k_i=k_i),'max':at_minmax(b,p,np.max,k_i=k_i)})

        else: # spin-polarized
            bu = self._data.bands.evals.SpinUp[:,b_i]
            pu = self._data.pro_bands.pros.SpinUp[:,:,b_i,:]

            _minu = at_minmax(bu,pu,np.min,k_i=k_i)
            _maxu = at_minmax(bu,pu,np.max,k_i=k_i)

            bd = self._data.bands.evals.SpinDown[:,b_i]
            pd = self._data.pro_bands.pros.SpinDown[:,:,b_i,:]

            _mind = at_minmax(bd,pd,np.min,k_i=k_i)
            _maxd = at_minmax(bd,pd,np.max,k_i=k_i)

            if isinstance(k_i,int): # single kpoint
                return serializer.Dict2Data({'SpinUp':_minu,'SpinDown':_mind})

            return serializer.Dict2Data({'SpinUp':{'min':_minu,'max':_maxu},'SpinDown':{'min':_mind,'max':_maxd}})

    def get_en_diff(self,b1_i,b2_i,k1_i=None,k2_i=None):
        """Get energy difference between two bands at given two kpoints indices. Index 2 is considered at higher energy.
        - b1_i, b2_i : band indices of the two bands, minimum energy difference is calculated.
        - k1_i, k2_i : k-point indices of the two bands.

        > If k1_i and k2_i are not provided, `min(b2_i) - max(b1_i)` is calculated which is equivalent to band gap.

        Returns: Data with follwoing attributes which can be used to annotate the difference on plot.
            de     : energy difference
            coords : np.array([[k1,e1],[k2,e2]]) #E_Fermi is subtracted either from system or when user provides in a plot command.
            eqv_coords: list(coords) at equivalent k-points if exit. Do not appear if k1_i and k2_i are provided.

        For spin-polarized case, 4 blocks of above data are returned which are accessible by
        `u1u2, u1d2, d1u2, d1d2` and they collects energy difference between 2 given bands at 2 different spin.

        """
        if k1_i and k2_i == None:
            raise ValueError('When you provide `k1_i`, `k2_i` cannot be None. They both can be None at same time.')
        if k1_i == None and k2_i:
            raise ValueError('When you provide `k2_i`, `k1_i` cannot be None. They both can be None at same time.')

        def format_coords(b1_max,b2_min):
            "maximum of b1 and min of b2 is taken for energy difference of two bands when kpoint not given."
            combs = []
            for k1 in b1_max.k:
                for k2 in b2_min.k:
                    combs.append(np.array([[k1,b1_max.e],[k2,b2_min.e]]))

            _out = {'coords':combs[0]}
            if combs[1:]:
                _out['eqv_coords'] = combs[1:]
            return _out

        if self._data.bands.ISPIN == 1:
            if isinstance(k1_i,int):
                b1 = self.get_band_info(b1_i,k_i=k1_i)
                b2 = self.get_band_info(b2_i,k_i=k2_i)
                return serializer.Dict2Data({'de':b2.e - b1.e, 'coords': np.array([[b1.k,b1.e],[b2.k, b2.e]])})
            else:
                b1 = self.get_band_info(b1_i,k_i=None).max
                b2 = self.get_band_info(b2_i,k_i=None).min

                return serializer.Dict2Data({'de':b2.e - b1.e, **format_coords(b1,b2)})
        else:
            if isinstance(k1_i,int):
                b1u = self.get_band_info(b1_i,k_i=k1_i).SpinUp
                b1d = self.get_band_info(b1_i,k_i=k1_i).SpinDown
                b2u = self.get_band_info(b2_i,k_i=k2_i).SpinUp
                b2d = self.get_band_info(b2_i,k_i=k2_i).SpinDown

                return serializer.Dict2Data({
                    'u1u2':{'de':b2u.e - b1u.e, 'coords':np.array([[b1u.k, b1u.e], [b2u.k, b2u.e]])},
                    'd1d2':{'de':b2d.e - b1d.e, 'coords':np.array([[b1d.k, b1d.e], [b2d.k, b2d.e]])},
                    'd1u2':{'de':b2u.e - b1d.e, 'coords':np.array([[b1d.k, b1d.e], [b2u.k, b2u.e]])},
                    'u1d2':{'de':b2d.e - b1u.e, 'coords':np.array([[b1u.k, b1u.e], [b2d.k, b2d.e]])}
                })
            else:
                b1u = self.get_band_info(b1_i,k_i=None).SpinUp.max # max in lower band
                b1d = self.get_band_info(b1_i,k_i=None).SpinDown.max
                b2u = self.get_band_info(b2_i,k_i=None).SpinUp.min # min in upper band
                b2d = self.get_band_info(b2_i,k_i=None).SpinDown.min

                return serializer.Dict2Data({
                    'u1u2':{'de':b2u.e - b1u.e, **format_coords(b1u,b2u)},
                    'd1d2':{'de':b2d.e - b1d.e, **format_coords(b1d,b2d)},
                    'd1u2':{'de':b2u.e - b1d.e, **format_coords(b1d,b2u)},
                    'u1d2':{'de':b2d.e - b1u.e, **format_coords(b1u,b2d)}
                })

    def splot_en_diff(self, coords, ax, **kwargs):
        """Plot energy difference at given ax. Provide `coords` from output of `get_en_diff().coords` or `get_en_diff().eqv_coords[i]` if exist.
        Provide `ax` on which bandstructure is plotted.
        E_Fermi is already subtracted in `coords` from system or by user input when bandstructure plot commands are run.
        kwargs are passed to `ax.step`.
        Returns ax.
        """
        kwargs = {'marker':'.',**kwargs}
        kwargs['where'] = 'mid' # override this
        coords = np.array(coords) # make sure it is np.array
        ax.step(*coords.T,**kwargs)
        return ax