__version__ = "0.8.2"

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