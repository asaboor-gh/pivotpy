__version__ = "0.9.3"

__all__ = []

from .s_plots import __all__ as sp_all
from .i_plots import __all__ as ip_all
from .g_utils import __all__ as gu_all
from .vr_parser import __all__ as vp_all
from .sio import __all__ as si_all
from .widgets import __all__ as wg_all
__all__.extend(vp_all)
__all__.extend(gu_all)
__all__.extend(sp_all)
__all__.extend(ip_all)
__all__.extend(si_all)
__all__.extend(wg_all)



# Access all functions through root modile pivotpy
from .s_plots import *
from .i_plots import *
from .g_utils import *
from .vr_parser import *
from .sio import *
from .widgets import *
from matplotlib.pyplot import show,savefig

mpl_imported=['show','savefig']
__all__.extend(mpl_imported)

# Register 'RGB' colormap in current session
from matplotlib.colors import LinearSegmentedColormap as LSC
import matplotlib.pyplot as plt, numpy as np
RGB = LSC.from_list('RGB',[(0.9,0,0),(0.9,0.9,0),(0,0.9,0),(0,0.9,0.9),(0,0,0.9)])
plt.register_cmap('RGB',RGB)

# color_marices for quick_rgb_lines
color_matrix = np.array([[0.5  , 0  , 0.5],[0.5  , 0.5, 0. ],[0.  , 0.5  , 0.5 ]])
gray_matrix = np.array([[0.3  , 0.59  , 0.11],[0.3  , 0.59  , 0.11],[0.3  , 0.59  , 0.11]])