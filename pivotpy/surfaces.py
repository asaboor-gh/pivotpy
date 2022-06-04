# AUTOGENERATED! DO NOT EDIT! File to edit: SpinProjectedSurfaces.ipynb (unless otherwise specified).

__all__ = ['SpinDataFrame']

# Cell
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
# Inside packages import to work both with package and jupyter notebook.
try:
    from pivotpy import parser as vp
    from pivotpy import api, serializer, splots
    from .splots import _validate_input
except:
    import pivotpy.parser as vp
    import pivotpy.api as api
    import pivotpy.serializer as serializer
    import pivotpy.splots as splots
    import pivotpy.splots._validate_input as _validate_input

# Cell
def _collect_spin_data(exported_spin_data, bands = [0], elements = [[0],], orbs = [[0],], scale_data = False, E_Fermi = None):
    if not isinstance(bands,(list,tuple)):
        raise TypeError('`bands` must be list/tuple of integer.')
    elements, orbs, _ = _validate_input(elements,orbs,[str(i) for i,o in enumerate(orbs)],sys_info = exported_spin_data.sys_info, rgb = False)

    fermi = E_Fermi or exported_spin_data.evals.E_Fermi
    def per_band_data(band):
        kpoints = exported_spin_data.kpoints
        evals = {k: v[:,band] for k,v in exported_spin_data.evals.items() if k in 'eud'}
        spins = {k: v[:,:,band,:] for k,v in exported_spin_data.spins.items() if k in 'sxyzud'}

        df_dict = {f'k{k}':v for k,v in zip('xyz',kpoints.T)}
        df_dict['band'] = [(band + exported_spin_data.evals.indices[0] + 1) for i in range(len(kpoints))]

        for k,v in evals.items():
            df_dict[k if k=='e' else f'e{k}'] = v.T.flatten() - fermi

        for i, (e, o) in enumerate(zip(elements,orbs)):
            for k,v in spins.items():
                if k in 'sxyzud':
                    key = k if k == 's' else f's{k}'
                    df_dict[f'{key}_{i}'] = v.take(e,axis=0).take(o,axis=2).sum(axis=2).sum(axis=0).T.flatten()

        return df_dict

    main_dict = per_band_data(bands[0])
    for b in bands[1:]:
        data = per_band_data(b)
        for k, v in main_dict.items():
            main_dict[k] = [*main_dict[k], *data[k]]

    if elements and scale_data == True: # Only scale if projections given
        _max = []
        for k,v in main_dict.items():
            if k.startswith('s'):
                _max.append(np.abs(v).max())

        _max = max(_max)
        for k,v in main_dict.items():
            if k.startswith('s'):
                main_dict[k] = v / (_max if _max != 0 else 1)
    return main_dict

# Cell
class SpinDataFrame(pd.DataFrame):
    """Spin data from vasprun.xml is converted to a dataframe.
    - **Parameters**:
        - path: path to `vasprun.xml` or auto picks in current directory.
        - bands: list of band indices [zero based here], In output data frame you will see corresponding band number based on full data.
        - elements: list of elements to plot. inner list contains ions indices. Can leave empty to discard projection data.
        - orbs: list of orbitals to plot. inner list contains orbitals indices. Can leave empty to discard projection data
        - scale_data: if True, spin data is scaled to -1 to 1.
        - E_Fermi: if not None, auto picked as Fermi level from vasprun.xml.
        - skipk: if not None, auto skipped unnecessary k-points.
        - elim: if not None, filtered out unnecessary bands.
        - data: if not None, data is loaded from given data/pickle/json/dict and validated. Many other parameters are ignored when data is given.

    - **Returns**:
        - SpinDataFrame: dataframe with colums as k-points, eigenvalues, spin components projected over selected ions and orbtials.

    - **Methods**:
        - sliced: Slice data in a plane orthogonal to given `column` at given `value`.
        - masked: Mask data over a constant value in a given column. Useful for plotting fermi level/surface.
        - splot: plot data in a 2D plot.
        - splot3d: plot data in a 3D plot.

        All other methods are inherited from pd.DataFrame. If you apply some method, then use `wraps` to wrap the result in a SpinDataFrame.
    """
    _metadata = ['_current_attrs','scale_data','sys_info','poscar'] # These are passed after operations to new dataframe.
    def __init__(self, *args, path = None, bands = [0], elements = [[0],], orbs = [[0],], scale_data = False, E_Fermi = None, elim = None, skipk=None, data = None, **kwargs):
        if not (path or args): # It works fine without path given, but it is not recommended.
            path = './vasprun.xml'
        if path or data: # Don't updates args otherwise
            spin_data = None # To avoid access before assignment
            if data:
                spin_data = serializer.SpinData.validated(data)
            elif isinstance(path, str):
                try: # Check if 4 sets of data are there.
                    spin_data = vp.export_spin_data(path, spins = 'sxyz',skipk = skipk, elim = elim)
                except: # If not 4 sets
                    spin_data = vp.export_spin_data(path, spins = 's',skipk = skipk, elim = elim)
            else:
                raise ValueError('Invalid path or data!')

            if spin_data:
                out_dict = _collect_spin_data(spin_data, bands = bands, elements = elements, orbs = orbs, scale_data = scale_data, E_Fermi = E_Fermi)
                super().__init__(out_dict)
                self.scale_data = scale_data
                self.sys_info = spin_data.sys_info
                # Path below is used to get kpoints info
                self.poscar = api.POSCAR(path = path, data = spin_data.poscar)
                self.poscar._kpts_info = spin_data.sys_info.kpts_info
                self._current_attrs = {'cmap':'viridis'} # To store attributes of current plot for use in colorbar.

        else: # This part is only for operations on dataframe.
            if len(args) == 1: # This gives hack to load data from a file in current directory.
                if (args[0] is None) or isinstance(args[0],str): # If path is given as positional argument
                    raise ValueError('SpinDataFrame expects no positional argument!')
            super().__init__(*args, **kwargs)

    @property
    def _constructor(self):
        "That's main hero of this class. This is called when you apply some method of slice it."
        return SpinDataFrame


    def masked(self, column, value, tol = 1e-2, n = None, band = None, method = 'cubic'):
        """Mask dataframe with a given value, using a tolerance.
        If n is given, (band should also be given to avoid multivalued interpolation error) data values are interpolated to grid of size (l,m,n) where n is longest side.
        n could be arbirarily large as mask will filter out data outside the tolerance."""
        if n and not isinstance(n,int):
            raise TypeError('`n` must be an integer to be applied to short side of grid.')
        if isinstance(n,int) and not isinstance(band,int):
            raise ValueError('A single `band`(int) from dataframe must be given to mask data.'
            'Interpolation does not give correct results for multiple values over a point (x,y,z).')

        column_vals = self[column].to_numpy()
        cmin, cmax = column_vals.min(), column_vals.max()
        if (value < cmin) or (value > cmax):
            raise ValueError('value is outside of column data range!')

        df = self.copy() # To avoid changing original dataframe
        if n and band:
            bands_exist  = np.unique(self.band)
            if band not in bands_exist:
                raise ValueError('Band {} is not in dataframe! Available: {}'.format(band, bands_exist))

            _self_ = df[df['band'] == band] # Single band dataframe
            df.drop(df.index, inplace=True)  # only columns names there and metadata
            kxyz = _self_[['kx','ky','kz']].to_numpy()
            lx,*_,hx = _self_['kx'].sort_values(inplace=False)
            ly,*_,hy = _self_['ky'].sort_values(inplace=False)
            lz,*_,hz = _self_['kz'].sort_values(inplace=False)
            vs = np.array([hx-lx,hy-ly, hz-lz])
            nx, ny, nz = nxyz = (vs/vs.max()*n).astype(int)
            nijk = [i for i,n in enumerate(nxyz) if n > 0]

            if len(nijk) < 2:
                raise ValueError('At least two of kx,ky,kz must have non-coplanar points.')

            if len(nijk) == 3:
                xyz = kxyz.T
                XYZ = [a.flatten() for a in np.mgrid[lx:hx:nx*1j,ly:hy:ny*1j,lz:hz:nz*1j]]
                for name, index in zip('xyz',range(3)):
                    df[f'k{name}'] = XYZ[index]
            else:
                [l1,l2],[h1,h2],[n1,n2] = np.array([[lx,ly,lz],[hx,hy,hz],[nx,ny,nz]])[:,nijk]
                xyz = kxyz.T[nijk]
                XYZ = [a.flatten() for a in np.mgrid[l1:h1:n1*1j,l2:h2:n2*1j]]
                for name, index in zip('xyz',range(3)):
                    if index in nijk:
                        df[f'k{name}'] = XYZ[index]
                    else:
                        df[f'k{name}'] = np.zeros_like(XYZ[0])

            for c in [_c for _c in _self_.columns if _c not in 'kxkykz']:
                df[c] = griddata(tuple(xyz),_self_[c].to_numpy(),tuple(XYZ),method = method)

            del _self_ # To free memory

            df = df.round(6).dropna()

            if self.scale_data == True:
                _max = []
                for k in [c for c in df.columns if c.startswith('s')]:
                    _max.append(np.abs(df[k]).max())

                _max = max(_max)
                for k in [c for c in df.columns if c.startswith('s')]:
                    df[k] = df[k] / (_max if _max != 0 else 1)

        # Make sure to keep metadata, it doesn't work otherwise.
        self.send_metadata(df)

        return df[np.logical_and((df[column] < value + tol),(df[column] > value-tol))]

    def send_metadata(self, target_spin_dataframe):
        "Copy metadata from this to another SpinDataFrame."
        for k in self._metadata:
            setattr(target_spin_dataframe, k, getattr(self, k))

    def sliced(self,column = 'kz', value = 0):
        "Slice data in a plane orthogonal to given `column` at given `value`"
        return self[self[column] == value]

    def _collect_arrows_data(self, arrows):
        arrows_data = []
        for arr in arrows:
            if arr not in ['',*self.columns]:
                raise ValueError(f'{arr!r} is not a column in the dataframe')
            arrows_data.append(self[arr] if arr else np.zeros_like(self['kx'].to_numpy()))

        return np.array(arrows_data).T

    def _collect_kxyz(self, *xyz, basis = None):
        "Return tuple(kxyz, k_order)"
        _kxyz = ['kx','ky','kz']
        kij = [_kxyz.index(a) for a in xyz if a in _kxyz]
        kxyz = self[['kx','ky','kz']].to_numpy()
        kxyz = self.poscar.bring_in_bz(kxyz, basis = basis)

        # Handle third axis as energy as well
        if len(xyz) == 3 and xyz[2].startswith('e'):
            kxyz[:,2] = self[xyz[2]].to_numpy()
            kij = [*kij, 2] # Add energy to kij, it must be of size 2 before

        return kxyz, kij

    def _validate_columns(self, *args):
        for arg in args:
            if arg not in self.columns:
                raise ValueError(f'{arg!r} is not a column in the dataframe')

    def splot(self,*args, arrows = [], every=1, norm = 1, marker='H', ax = None, quiver_kws = {}, basis = None, **kwargs):
        """Plot energy in 2D with/without arrows.
        - **Parameters**:
            - *args: 3 or 4 names of columns, representing [X,Y,Energy,[Anything]], from given args, last one is colormapped. If kwargs has color, that takes precedence.
            - arrows: 2 or 3 names of columns, representing [U,V,[color]]. If quiver_kws has color, that takes precedence.
            - every: every nth point is plotted as arrow.
            - norm: normalization factor for size of arrows.
            - marker: marker to use for scatter, use s as another argument to change size.
            - ax: matplotlib axes to plot on (defaults to auto create one).
            - quiver_kws: these are passed to matplotlib.pyplot.quiver.
            - basis:If given, used to transform kpoints into BZ overriding default behavior.
        **kwargs are passed to matplotlib.pyplot.scatter.

        - **Returns**:
            - ax: matplotlib axes. It has additinal method `colorbar` to plot colorbar from most recent plot.

        See examples at https://massgh.github.io/pivotpy/
        """
        if arrows and len(arrows) not in [2,3]:
            raise ValueError('`arrows ` requires 2 or 3 items form spin data [s1,s2,[color]], one of s1,s2 could be "".')
        if len(args) not in [3,4]:
            raise ValueError('splot takes 3 or 4 positional arguments [X,Y,E,[Anything]], last one is colormapped if kwargs don\'t have color.')

        self._validate_columns(*args)
        kxyz, kij = self._collect_kxyz(*args[:2], basis = basis)
        ax = ax or api.get_axes()
        minmax_c = [0,1]
        cmap = kwargs.get('cmap',self._current_attrs['cmap'])

        if arrows:
            arrows_data = self._collect_arrows_data(arrows)
            cmap = quiver_kws.get('cmap',cmap)
            if 'color' in quiver_kws:
                cmap = None # No colorbar for color only
                arrows_data = arrows_data[:,:2] # color takes precedence
            ax.quiver(*kxyz[::every].T[kij],*(norm*arrows_data[::every].T), **quiver_kws)

            if len(arrows) == 3 and 'color' not in quiver_kws:
                cmap = quiver_kws.get('cmap','viridis') # Fallback to default
                minmax_c = [arrows_data[:,2].min(),arrows_data[:,2].max()]
        else:
            _C = self[args[-1]] # Most right arg is color mapped
            kwargs['marker'] = marker # Avoid double marker
            if 'color' in kwargs:
                kwargs['c'] = kwargs['color']
                del kwargs['color'] # keep one
                cmap = None # No colorbar

            kwargs['c'] = kwargs.get('c',_C)
            ax.scatter(*kxyz.T[kij],**kwargs)
            minmax_c = [min(_C),max(_C)]

        self._current_attrs = {'ax':ax,'minmax_c':minmax_c,'cmap':cmap}
        return ax

    def splot3d(self,*args, arrows = [], every=1,norm = 1, marker='H', ax = None, quiver_kws = {'arrowstyle':'-|>','size':1}, basis = None, **kwargs):
        """Plot energy in 3D with/without arrows.
        - **Parameters**:
            - *args: 3, 4 or 5 names of columns, representing [X,Y,[Z or Energy],Energy, [Anything]], out of given args, last one is color mapped. if kwargs has color, that takes precedence.
            - arrows: 3 or 4 names of columns, representing [U,V,W,[color]]. If color is not given, magnitude of arrows is color mapped. If quiver_kws has color, that takes precedence.
            - every: every nth point is plotted as arrow.
            - norm: normalization factor for size of arrows.
            - marker: marker to use for scatter, use s as another argument to change size.
            - ax: matplotlib 3d axes to plot on (defaults to auto create one).
            - quiver_kws: these are passed to pivotpy.fancy_quiver3d.
            - basis:If given, used to transform kpoints into BZ overriding default behavior.
        **kwargs are passed to matplotlib.pyplot.scatter.

        - **Returns**:
            - ax: matplotlib 3d axes. It has additinal method `colorbar` to plot colorbar from most recent plot.

        See examples at https://massgh.github.io/pivotpy/
        """
        if arrows and len(arrows) not in [3,4]:
            raise ValueError('`arrows ` requires 3 or 4 items form spin data [s1,s2, s2, [color]], one of s1,s2,s3 could be "".')
        if len(args) not in [3,4,5]:
            raise ValueError('splot3d takes 3, 4 or 5 positional arguments [X,Y,E] or [X,Y,Z,E,[Anything]], right-most is color mapped if kwargs don\'t have color.')

        if not args[2][0] in 'ek':
            raise ValueError('Z axis must be in [kx,ky,kz, energy]!')

        self._validate_columns(*args)
        kxyz, kij = self._collect_kxyz(*args[:3], basis = basis)
        ax = ax or api.get_axes(axes_3d=True)
        minmax_c = [0,1]
        cmap = kwargs.get('cmap',self._current_attrs['cmap'])

        if arrows:
            arrows_data = self._collect_arrows_data(arrows)
            cmap = quiver_kws.get('cmap',cmap)
            if 'color' in quiver_kws: # color takes precedence
                quiver_kws['C'] = quiver_kws['color']
                quiver_kws.pop('color') # It is not in FancyArrowPatch
                cmap = None # No colorbar
            elif len(arrows) == 4:
                array = arrows_data[::every,3]
                array = (array - array.min())/np.ptp(array)
                quiver_kws['C'] = plt.get_cmap(cmap)(array)
                minmax_c = [arrows_data[:,3].min(),arrows_data[:,3].max()]
            elif len(arrows) == 3:
                array = np.linalg.norm(arrows_data[::every,:3],axis=1)
                minmax_c = [array.min(),array.max()] # Fist set then normalize
                array = (array - array.min())/np.ptp(array)
                quiver_kws['C'] = plt.get_cmap(cmap)(array)

            if 'cmap' in quiver_kws:
                quiver_kws.pop('cmap') # It is not in fancy_quiver3d

            api.fancy_quiver3d(*kxyz[::every].T[kij],*(norm*arrows_data[::every].T[:3]), **quiver_kws,ax=ax)

        else:
            _C = self[args[-1]] # Most righht arg is color mapped
            kwargs['marker'] = marker # Avoid double marker
            if 'color' in kwargs:
                kwargs['c'] = kwargs['color']
                del kwargs['color'] # keep one
                cmap = None # No colorbar

            kwargs['c'] = kwargs.get('c',_C)
            ax.scatter(*kxyz.T[kij],**kwargs)
            minmax_c = [min(_C),max(_C)]

        self._current_attrs = {'ax':ax,'minmax_c':minmax_c,'cmap':cmap}
        return ax

    def colorbar(self, cax = None, nticks = 6, digits = 2, **kwargs):
        " Add colobar to most recent plot. kwargs are passed to pivotpy.splots.add_colorbar"
        if not self._current_attrs['ax']:
            raise ValueError('No plot has been made yet by using `splot, splot3d` or already consumed by `colorbar`')
        if not self._current_attrs['cmap']:
            raise ValueError('No Mappable for colorbar found!')

        ax = self._current_attrs['ax']
        cmap = self._current_attrs['cmap']
        minmax_c = self._current_attrs['minmax_c']
        self._current_attrs['ax'] = None # Reset
        self._current_attrs['cmap'] = None # Reset
        if ax.name == '3d':
            cax = cax or plt.gcf().add_axes([0.85, 0.15, 0.03, 0.7])

        return splots.add_colorbar(ax, cax, cmap ,ticks = np.linspace(*minmax_c,nticks,endpoint=True),digits = digits, **kwargs)

