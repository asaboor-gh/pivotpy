# AUTOGENERATED! DO NOT EDIT! File to edit: SpinProjectedSurfaces.ipynb (unless otherwise specified).

__all__ = ['SpinDataFrame']

# Cell
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
# Inside packages import to work both with package and jupyter notebook.
try:
    from pivotpy import vr_parser as vp
    from pivotpy import api
    from .splots import _validate_input
except:
    import pivotpy.vr_parser as vp
    import pivotpy.api as api
    import pivotpy.splots._validate_input as _validate_input

# Cell

#export
def _collect_spin_data(exported_spin_data, band = 0, elements = [[0],], orbs = [[0],], scale_data = False, E_Fermi = None):
    if not isinstance(band,int):
        raise TypeError('`band` must be an integer.')
    elements, orbs, _ = _validate_input(elements,orbs,[str(i) for i,o in enumerate(orbs)],sys_info = exported_spin_data.sys_info, rgb = False)

    fermi = E_Fermi or exported_spin_data.evals.E_Fermi

    kpoints = exported_spin_data.kpoints
    evals = {k: v[:,band] for k,v in exported_spin_data.evals.items() if k in 'eud'}
    spins = {k: v[:,:,band,:] for k,v in exported_spin_data.spins.items() if k in 'sxyzud'}

    df_dict = {f'k{k}':v for k,v in zip('xyz',kpoints.T)}

    for k,v in evals.items():
        df_dict[k if k=='e' else f'e{k}'] = v.T.flatten() - fermi

    for i, (e, o) in enumerate(zip(elements,orbs)):
        for k,v in spins.items():
            if k in 'sxyzud':
                key = k if k == 's' else f's{k}'
                df_dict[f'{key}_{i}'] = v.take(e,axis=0).take(o,axis=2).sum(axis=2).sum(axis=0).T.flatten()

    if scale_data == True:
        _max = []
        for k,v in df_dict.items():
            if k.startswith('s'):
                _max.append(np.abs(v).max())

        _max = max(_max)
        for k,v in df_dict.items():
            if k.startswith('s'):
                df_dict[k] = v / (_max if _max != 0 else 1)
    return df_dict

class SpinDataFrame(pd.DataFrame):
    """Spin data from vasprun.xml is converted to a dataframe.
    - **Parameters**:
        - path: path to `vasprun.xml` or auto picks in current directory.
        - band: band index to plot, only one band is supported.
        - elements: list of elements to plot. inner list contains ions indices.
        - orbs: list of orbitals to plot. inner list contains orbitals indices.
        - scale_data: if True, data is scaled to -1 to 1.
        - E_Fermi: if not None, auto picked as Fermi level from vasprun.xml.
        - skipk: if not None, auto skipped unnecessary k-points.

    - **Returns**:
        - SpinDataFrame: dataframe with colums as k-points, eigenvalues, spin components projected over selected ions and orbtials.

    - **Methods**:
        - wraps: wraps a pd.DataFrame in a SpinDataFrame, that enables extra methods.
        - sliced: Slice data in a plane orthogonal to given `column` at given `value`.
        - masked: Mask data over a constant value in a given column. Useful for plotting fermi level/surface.
        - splot: plot data in a 2D plot.

        All other methods are inherited from pd.DataFrame. If you apply some method, then use `wraps` to wrap the result in a SpinDataFrame.
    """
    def __init__(self,path = None, band = 0, elements = [[0],], orbs = [[0],], scale_data = False, E_Fermi = None,skipk=None):
        if not isinstance(path,(pd.DataFrame,self.__class__)):
            try: # Check if 4 sets of data are given.
                spin_data = vp.export_spin_data(path=path,spins='sxyz',skipk=skipk)
            except: # If not 4 sets
                spin_data = vp.export_spin_data(path=path,spins='s',skipk=skipk)
            out_dict = _collect_spin_data(spin_data, band = band, elements = elements, orbs = orbs, scale_data = scale_data, E_Fermi = E_Fermi)
            super().__init__(out_dict)
            self.cmap = 'viridis'
            self.scale_data = scale_data
            self.poscar = api.POSCAR(text_plain = spin_data.poscar.text_plain)
            self.sys_info = spin_data.sys_info
        else:
            super().__init__(path.to_dict()) # Works for both SpinDataFrame and pd.DataFrame

    def wraps(self,dataframe):
        "wraps a pd.DataFrame in a SpinDataFrame, that enables extra methods"
        df = self.__class__(dataframe)
        df.poscar = self.poscar
        df.sys_info = self.sys_info
        df.cmap = self.cmap
        df.scale_data = self.scale_data
        return df

    def move2bz(self, primitive = False):
        raise Exception('Use `self.poscar.bring_in_bz(koints)` instead for your custom functions, Works fine in `splot` without it.')

    def masked(self, column, value, tol = 1e-2, n = None, method = 'cubic'):
        "Mask dataframe with a given value, using a tolerance, If n is given, data values are interpolated to grid of size (l,m,n) where n is longest side."
        if n and not isinstance(n,int):
            raise TypeError('`n` must be an integer to be applied to short side of grid.')

        df =  self
        if n:
            kxyz = self[['kx','ky','kz']].to_numpy()
            lx,*_,hx = self['kx'].sort_values(inplace=False)
            ly,*_,hy = self['ky'].sort_values(inplace=False)
            lz,*_,hz = self['kz'].sort_values(inplace=False)
            vs = np.array([hx-lx,hy-ly, hz-lz])
            nx, ny, nz = nxyz = (vs/vs.max()*n).astype(int)
            nijk = [i for i,n in enumerate(nxyz) if n > 0]

            if len(nijk) < 2:
                raise ValueError('At least two of kx,ky,kz must have non-coplanar points.')
            df_dict = {}
            if len(nijk) == 3:
                xyz = kxyz.T
                XYZ = [a.flatten() for a in np.mgrid[lx:hx:nx*1j,ly:hy:ny*1j,lz:hz:nz*1j]]
                for name, index in zip('xyz',range(3)):
                    df_dict[f'k{name}'] = XYZ[index]
            else:
                [l1,l2],[h1,h2],[n1,n2] = np.array([[lx,ly,lz],[hx,hy,hz],[nx,ny,nz]])[:,nijk]
                xyz = kxyz.T[nijk]
                XYZ = [a.flatten() for a in np.mgrid[l1:h1:n1*1j,l2:h2:n2*1j]]
                for name, index in zip('xyz',range(3)):
                    if index in nijk:
                        df_dict[f'k{name}'] = XYZ[index]
                    else:
                        df_dict[f'k{name}'] = np.zeros_like(XYZ[0])

            for c in [_c for _c in self.columns if _c not in 'kxkykz']:
                df_dict[c] = griddata(tuple(xyz),self[c].to_numpy(),tuple(XYZ),method = method)

            df = pd.DataFrame(df_dict).round(6).dropna()

            if self.scale_data == True:
                _max = []
                for k in [c for c in df.columns if c.startswith('s')]:
                    _max.append(np.abs(df[k]).max())

                _max = max(_max)
                for k in [c for c in df.columns if c.startswith('s')]:
                    df[k] = df[k] / (_max if _max != 0 else 1)

        return self.wraps(df[np.logical_and((df[column] < value + tol),(df[column] > value-tol))])

    def sliced(self,column = 'kz', value = 0):
        "Slice data in a plane orthogonal to given `column` at given `value`"
        return self.wraps(self[self[column] == value])

    def splot(self,*args, arrows = [], every=4, marker='H', ax = None, **kwargs):
        """plot a 2D energy surface with/without arrows.
        - **Parameters**:
            - *args: 3 or 4 names of columns, representing X,Y,Energy, and optionally, something to colorize data.
            - arrows: 2 or 3 names of columns, representing U,V as arrows direstcion, and optionally, something to colorize data.
            - every: every nth point is plotted as arrow.
            - marker: marker to use for scatter, use s as another argument to change size.
            - ax: matplotlib axes to plot on (defaults to auto create one).
        **kwargs are passed to matplotlib.pyplot.scatter if arrows are not given, otherwise to matplotlib.pyplot.quiver.

        - **Returns**:
            - ax: matplotlib axes. It has additinal method `colorbar` to plot colorbar from most recent plot.

        See examples at https://massgh.github.io/pivotpy/
        """
        if arrows and len(arrows) not in [2,3]:
            raise ValueError('`arrows ` requires 2 or 3 items form spin data [s1,s2,[color]], one of s1,s2 could be "".')
        if len(args) not in [3,4]:
            raise ValueError('splot takes 3 or 4 positional arguments [X,Y,E,[color]]')

        for arg in args:
            if arg not in self.columns:
                raise ValueError(f'{arg!r} is not a column in the dataframe')

        arrows_data = []
        for i, arr in enumerate(arrows):
            if arr not in ['',*self.columns]:
                raise ValueError(f'{arr!r} is not a column in the dataframe')
            arrows_data.append(self[arr] if arr else np.zeros_like(self['kx'].to_numpy()))

        arrows_data = np.array(arrows_data).T

        kij = [['kx','ky','kz'].index(a) for a in args[:2]]
        kxyz = self.poscar.bring_in_bz(self[['kx','ky','kz']].to_numpy())
        ax = ax or api.get_axes()
        minmax_c = [0,1]
        self.cmap = kwargs.get('cmap',self.cmap)

        if arrows:
            ax.quiver(*kxyz[::every].T[kij],*arrows_data[::every].T, **kwargs)
            if len(arrows) == 3:
                minmax_c = [arrows_data[:,2].min(),arrows_data[:,2].max()]
        else:
            _C = self[args[3]] if len(args) == 4 else self[args[2]]
            ax.scatter(*kxyz.T[kij],c = _C, marker=marker,**kwargs)
            minmax_c = [min(_C),max(_C)]

        cmap = kwargs.get('cmap',self.cmap)
        def colorbar(cax = None, nticks = 6, **kwargs):
            ticks = np.linspace(0,1,nticks, endpoint=True)
            labels = np.linspace(*minmax_c,nticks,endpoint=True).round(2).astype(str)
            return ax.add_colorbar(cax, cmap ,ticks = list(ticks), ticklabels=labels, **kwargs)

        ax.colobar = colorbar
        return ax