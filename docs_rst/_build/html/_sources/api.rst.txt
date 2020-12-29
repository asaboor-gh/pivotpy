Pivotpy API Reference
***********************

Vasprun Parser
==============

.. autofunction:: pivotpy.vr_parser.Dict2Data
.. autofunction:: pivotpy.vr_parser.read_asxml
.. autofunction:: pivotpy.vr_parser.get_summary
.. autofunction:: pivotpy.vr_parser.get_kpts
.. autofunction:: pivotpy.vr_parser.get_tdos
.. autofunction:: pivotpy.vr_parser.get_evals
.. autofunction:: pivotpy.vr_parser.get_bands_pro_set
.. autofunction:: pivotpy.vr_parser.get_dos_pro_set
.. autofunction:: pivotpy.vr_parser.get_structure
.. autofunction:: pivotpy.vr_parser.export_vasprun
.. autofunction:: pivotpy.vr_parser.load_export
.. autofunction:: pivotpy.vr_parser.dump_dict
.. autofunction:: pivotpy.vr_parser.load_from_dump
.. autofunction:: pivotpy.vr_parser.islice2array
.. autofunction:: pivotpy.vr_parser.slice_data
.. autofunction:: pivotpy.vr_parser.split_vasprun

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