Pivotpy API Reference
***********************

.. toctree::
    :maxdepth: 2

    vr

.. plot::

   import pivotpy as pp, numpy as np 
   import matplotlib.pyplot as plt 
   vr1=pp.export_vasprun('E:/Research/graphene_example/ISPIN_2/bands/vasprun.xml')
   vr2=pp.export_vasprun('E:/Research/graphene_example/ISPIN_2/dos/vasprun.xml')
   axs=pp.init_figure(ncols=3,widths=[2,1,2.2],sharey=True,wspace=0.05,figsize=(8,2.6))
   elements=[0,[0],[0,1]]
   orbs=[[0],[1],[2,3]]
   labels=['s','$p_z$','$(p_x+p_y)$']
   ti_cks=dict(xt_indices=[0,30,60,-1],xt_labels=['Γ','M','K','Γ'])
   args_dict=dict(elements=elements,orbs=orbs,labels=labels,elim=[-20,15])
   pp.quick_bplot(path_evr=vr1,ax=axs[0],**ti_cks,elim=[-20,15])
   pp.quick_rgb_lines(path_evr=vr1,ax=axs[2],**args_dict,**ti_cks,colorbar=True,)
   pp.quick_dos_lines(path_evr=vr2,ax=axs[1],vertical=True,spin='both',include_dos='pdos',**args_dict,legend_kwargs={'ncol': 3},color_map='RGB_m')
   pp.color_wheel(axs[2],xy=(0.7,1.15),scale=0.2,labels=[l+'$^{⇅}$' for l in labels])
   pp._show()


.. ipython::

    In [12]: import numpy.random

    In [13]: numpy.random.rand(10,2)


Static Plots
=================

.. automodule:: pivotpy.s_plots
    :members: init_figure, quick_bplot, quick_rgb_lines, quick_color_lines, quick_dos_lines, plt_to_html, plot_potential

Interactive Plots
=================

.. automodule:: pivotpy.i_plots
    :members: plotly_to_html, plotly_rgb_lines, plotly_dos_lines

Utilities
=================

.. automodule:: pivotpy.g_utils
    :members: get_file_size, interpolate_data, ps_to_py, ps_to_std, get_child_items, EncodeFromNumpy, DecodeToNumpy, Vasprun, export_outcar, export_potential, LOCPOT_CHG, transform_color

Interactive Widgets
===================

.. automodule:: pivotpy.widgets
    :members: get_files_gui, get_input_gui, generate_summary, show_app

Structue 
=================
.. automodule:: pivotpy.sio
    :members: save_mp_API, get_poscar, get_kpath, get_kmesh, get_bz, plot_bz