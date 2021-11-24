# AUTOGENERATED! DO NOT EDIT! File to edit: MainAPI.ipynb (unless otherwise specified).

__all__ = ['__all__', 'parse_text', 'POSCAR', 'LOCPOT', 'get_axes', 'Vasprun']

# Cell
import os
import numpy as np
import plotly.graph_objects as go
try:
    from pivotpy import vr_parser as vp
    from pivotpy import s_plots as sp
    from pivotpy import i_plots as ip
    from pivotpy import sio as sio
    from pivotpy import widgets as wdg
    from pivotpy import g_utils as gu
except:
    import pivotpy.vr_parser as vp
    import pivotpy.s_plots as sp
    import pivotpy.i_plots as ip
    import pivotpy.sio as sio
    import pivotpy.widgets as wdg
    import pivotpy.g_utils as gu

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
# Direct function exports from modules
_memebers = (
    sio.InvokeMaterialsProject,
    sio.get_kpath,
    sio.str2kpath,
    sio.fancy_quiver3d,
    sio.rotation,
    wdg.generate_summary,
    wdg.VasprunApp,
    wdg.KPathApp,
    gu.set_dir,
    gu.get_child_items,
    gu.color,
    gu.transform_color,
    gu.interpolate_data,
    vp.split_vasprun,
    vp.xml2dict,
    ip.iplot2html,
    sp.plt2html,
    sp.plt2text,
    sp.show,
    sp.savefig
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
class POSCAR:
    "POSACR class to contain data and related methods"
    def __init__(self,path=None,content=None,_other_data=None):
        """Do not use `_other_data` yourself, it's for operations on poscar.
        Prefrence order: _other_data, content, path"""
        self.path = path
        self.content = content
        if _other_data:
            self._data = _other_data
        else:
            self._data = sio.export_poscar(path=path,content=content)
        # These after data to work with data
        self.primitive = False
        self._bz = self.get_bz(primitive = False) # Get defualt regular BZ
        self._cell = self.get_cell() # Get defualt cell

    def __repr__(self):
        return repr(self._data)

    @property
    def data(self):
        "Data object in POSCAR."
        return self._data
    @property
    def poscar(self):
        "poscar data."
        return self._data

    @_sub_doc(sio.get_bz,'- path_pos')
    def get_bz(self, loop=True, digits=8, primitive=False):
        self._bz = sio.get_bz(path_pos=self._data.basis, loop=loop, digits=digits, primitive=primitive)
        self.primitive = primitive
        return self._bz

    def set_bz(self,primitive=False,loop=True,digits=8):
        """Set BZ in primitive or regular shape. returns None, just set self._bz"""
        self.get_bz(primitive=primitive,loop=loop,digits=digits)

    def get_cell(self, loop=True, digits=8):
        "See docs of `get_bz`, same except space is inverted."
        self._cell = sio.get_bz(path_pos=self._data.rec_basis,loop=loop, digits=digits, primitive=True) # cell must be primitive
        return self._cell

    @_sub_doc(sio.splot_bz,'- path_pos_bz')
    def splot_bz(self, ax=None, plane=None, color='blue', fill=True, vectors=True, v3=False, vname='b', colormap='plasma', light_from=(1, 1, 1), alpha=0.4):
        return sio.splot_bz(path_pos_bz = self._bz, ax=ax, plane=plane, color=color, fill=fill, vectors=vectors, v3=v3, vname=vname, colormap=colormap, light_from=light_from, alpha=alpha)

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
        return POSCAR(_other_data = sio.join_poscars(poscar1=self._data, poscar2=other.data, direction=direction, tol=tol))

    @_sub_doc(sio.scale_poscar,'- path_poscar')
    def scale(self, scale=(1, 1, 1), tol=0.01):
        return POSCAR(_other_data = sio.scale_poscar(path_poscar=self._data, scale=scale, tol=tol))

    @_sub_doc(sio.rotate_poscar,'- path_poscar')
    def rotate(self,angle_deg,axis_vec):
        return POSCAR(_other_data = sio.rotate_poscar(path_poscar=self._data, angle_deg = angle_deg, axis_vec=axis_vec))

    @_sub_doc(sio.fix_sites,'- poscar')
    def fix_sites(self, tol=0.01, eqv_sites=True, translate=None):
        return POSCAR(_other_data = sio.fix_sites(poscar=self._data, tol=tol, eqv_sites=eqv_sites, translate=translate))

    @_sub_doc(sio.get_kmesh,'- path_pos')
    def get_kmesh(self, n_xyz=[5, 5, 5], weight=None, ibzkpt=None, outfile=None):
        return sio.get_kmesh(n_xyz=n_xyz, weight=weight, ibzkpt=ibzkpt, path_pos=self._data.basis, outfile=outfile)

    def bring_in_cell(self,points):
        """Brings atoms's positions inside Cell and returns their R3 coordinates."""
        return sio.to_R3(self._data.basis, points= points)

    @_sub_doc(sio.kpoints2bz,'- bz')
    def bring_in_bz(self,kpoints):
        if not self._bz:
            return print('No BZ found. Please run get_bz() first.')
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

    def __repr__(self):
        return repr(self._data)

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
            return print("Magnetization data set does not exist in {}".format(self.path))
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
            return print("operation expects any of {!r}, got {}".format(check,operation))
        if e_or_m is None:
            try:
                data = self._data.e
            except:
                try:
                    data = self._data.m
                except:
                   data = self._data.to_dict()[f'm_{self.m}']
                else:
                    return print("Magnetization data set does not exist in {}".format(self.path))
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
        for f in [sp.add_text,sp.add_legend,sp.add_colorbar,sp.color_wheel,sp.break_spines,sp.modify_axes]:
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
'''

# Cell
class Vasprun:
    """
    - All plotting functions that depend on `export_vasprun` are joined under this class and renamed.
    - **Parameters**
        - path       : str: path/to/vasprun.xml. Auto picks in CWD.
        - skipk      : int: Skip initial kpoints
        - elim       : list: Energy range e.g. [-5,5].
        - shift_kpath: float: Shift in kpath values for side by side plotting.
    - **Attributes**
        - data : Return of `export_vasprun` which is auto-picked in plotting methods under this class.

    > Tip: If KPOINTS file is generated by this module, ticks on kpath are auto-picked.
    [See Docs](https://massgh.github.io/pivotpy/)
    """
    def __init__(self,path=None,skipk=None,elim=[],shift_kpath=0,try_pwsh=True):
        self._data = vp.export_vasprun(path=path,skipk=skipk,elim=elim,shift_kpath=shift_kpath,try_pwsh=try_pwsh)
        self.elim = elim
        if path == None:
            kfile = 'KPOINTS'
        else:
            kfile = os.path.join(os.path.dirname(path),'KPOINTS')
        self.kticks = sio.read_ticks(kfile)

    def __handle_kwargs(self,kwargs,dos=False):
        kwargs = {'elim': self.elim, **kwargs}
        if dos:
            return kwargs
        ticks = {k:self.kticks[k] for k in ['ktick_inds','ktick_vals','kseg_inds']}
        kwargs = {**ticks,**kwargs} #Prefer provided ones
        return kwargs

    def __repr__(self):
        return repr(self._data)

    @property
    def data(self):
        return self._data

    @_sub_doc(sp.splot_bands,'- path_evr')
    def splot_bands(self,ax = None,**kwargs):
        kwargs = self.__handle_kwargs(kwargs)
        return sp.splot_bands(self._data,ax = ax, **kwargs)

    @_sub_doc(sp.splot_dos_lines,'- path_evr')
    def splot_dos_lines(self,elements = [[0],], orbs = [[0],], labels = ['s',], ax = None,**kwargs):
        kwargs = self.__handle_kwargs(kwargs,dos=True)
        return sp.splot_dos_lines(self._data,elements = elements, orbs = orbs, labels = labels, ax = ax, **kwargs)

    @_sub_doc(sp.splot_rgb_lines,'- path_evr')
    def splot_rgb_lines(self,elements = [[],[],[]], orbs = [[],[],[]], labels = ['','',''], ax = None, **kwargs):
        kwargs = self.__handle_kwargs(kwargs)
        return sp.splot_rgb_lines(self._data,elements = elements, orbs = orbs, labels = labels, ax = ax, **kwargs)

    @_sub_doc(sp.splot_color_lines,'- path_evr')
    def splot_color_lines(self,elements = [[0],], orbs = [[0],], labels = ['s',],axes = None,**kwargs):
        kwargs = self.__handle_kwargs(kwargs)
        return sp.splot_color_lines(self._data,elements = elements, orbs = orbs, labels = labels, axes = axes, **kwargs)

    @_sub_doc(ip.iplot_dos_lines,'- path_evr')
    def iplot_dos_lines(self,elements = [[0],], orbs = [[0],], labels = ['s',],**kwargs):
        kwargs = self.__handle_kwargs(kwargs, dos=True)
        return ip.iplot_dos_lines(self._data,elements = elements, orbs = orbs, labels = labels, **kwargs)

    @_sub_doc(ip.iplot_rgb_lines,'- path_evr')
    def iplot_rgb_lines(self,elements = [[],[],[]], orbs = [[],[],[]], labels = ['','',''],**kwargs):
        kwargs = self.__handle_kwargs(kwargs)
        return ip.iplot_rgb_lines(self._data,elements = elements, orbs = orbs, labels = labels, **kwargs)