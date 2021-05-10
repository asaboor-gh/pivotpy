"""Pivotpy is a processing tool for VASP DFT input/output processing.
Author: Abdul Saboor
Licence: Apache License Version 2.0, January 2004 #See file

Modules:
-------
    pivotpy.vr_parser #Parser for vasprun.xml
    pivotpy.g_utils   # general functions
    pivotpy.s_plots   # Matplotlib's plotting functions
    pivotpy.i_plots   # Plotly's interactive plotting functions
    pivotpy.sio       # Functions for BZ, POSCAR, KPath etc.
    pivotpy.widgets   # Jupyter widgets apps for output analysis/kpath selection. Very useful in big projectes. 
    
Usage:
-----
    All modules are imported under a single namespace, you could do
    >>> import pivotpy as pp 
    ... pp.__all__ #To see what is available
    ... pp.generate_summary(paths)  # Get a dataframe for whole project after you used pp.VasprunApp
    
    If you want to acess private functions/varaiables, you need to import a submodule itself, e.g.
    >>> import pivotpy.s_plots as sp 
    ... sp._plot_bands() #This is a private function, (see leading underscore)
"""

__version__ = "1.0.3"

__all__ = []

from .g_utils import __all__ as gu_all
from .vr_parser import __all__ as vp_all
from .s_plots import __all__ as sp_all
from .i_plots import __all__ as ip_all
from .sio import __all__ as si_all
from .widgets import __all__ as wg_all
__all__.extend(vp_all)
__all__.extend(gu_all)
__all__.extend(sp_all)
__all__.extend(ip_all)
__all__.extend(si_all)
__all__.extend(wg_all)



# Access all functions through root modile pivotpy
from .g_utils import *
from .vr_parser import *
from .s_plots import *
from .i_plots import *
from .sio import *
from .widgets import *
from matplotlib.pyplot import show as _show,savefig as _savefig

mpl_imported=['_show','_savefig']
__all__.extend(mpl_imported)

# Register 'RGB' colormap in current session
from matplotlib.colors import LinearSegmentedColormap as __LSC
import matplotlib.pyplot as __plt, numpy as __np
RGB = __LSC.from_list('RGB',[(0.9,0,0),(0.9,0.9,0),(0,0.9,0),(0,0.9,0.9),(0,0,0.9)])
__plt.register_cmap('RGB',RGB)

# color_marices for quick_rgb_lines
color_matrix = __np.array([[0.5, 0, 0.5, 1],[0.5, 0.5, 0., 1],[0., 0.5, 0.5,0 ]])
gray_matrix = __np.array([[1,0,0,0],[0,1,0,1],[0,0,1,0]])

#Backward Compatibility
__mapping = {
    'quick_bplot':      'splot_bands',
    'quick_rgb_lines':  'splot_rgb_lines',
    'quick_color_lines':'splot_color_lines',
    'quick_dos_lines':  'splot_dos_lines',
    'plotly_rgb_lines': 'iplot_dos_lines',
    'plotly_dos_lines': 'iplot_dos_lines'
}
for k,v in __mapping.items():
    _code = f"""def {k}(*args,**kwargs):
    "See docs and arguments of {v!r} for reference to use in this function."
    print(color.yb("Name {k!r} is deprecated, use {v!r} in future."))
    return {v}(*args,**kwargs)"""
    exec(_code)

#Aliases
sbands = splot_bands 
srgb   = splot_rgb_lines
scolor = splot_color_lines
sdos   = splot_dos_lines
idos   = iplot_dos_lines
irgb   = iplot_rgb_lines

# Edit rcParams here
import matplotlib as __mpl
from cycler import cycler as __cycler
__mpl.rcParams.update(
    {
        'figure.dpi': 144, #Better to See
        'figure.autolayout': True,
        'axes.prop_cycle': __cycler(color=['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']),
        'axes.linewidth': 0.4, #set the value globally
        'font.serif': "STIXGeneral",
        'font.family': "serif",
        'mathtext.fontset': "stix"
    }
)
