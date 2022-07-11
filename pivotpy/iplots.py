# AUTOGENERATED! DO NOT EDIT! File to edit: InteractivePlots.ipynb (unless otherwise specified).

__all__ = ['iplot2html', 'iplot_rgb_lines', 'iplot_dos_lines']

# Cell
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import plotly.graph_objects as go

# Inside packages import to work both with package and jupyter notebook.
try:
    from pivotpy import parser as vp
    from pivotpy import splots as sp
    from pivotpy import utils as gu
except:
    import pivotpy.parser as vp
    import pivotpy.splots as sp
    import pivotpy.utils as gu

# Cell
def _get_rgb_data(
    kpath       = None,
    evals_set   = None,
    pros_set    = None,
    elements    = [[0],[],[]],
    orbs        = [[0],[],[]],
    interp_nk   = {},
    occs_set    = None,
    kpoints     = None,
    ):
    """
    - Returns a formatted RGB colored data to pass into `_rgb2plotly` function. Two arguments, `elements` and `orbs` should be in one-to-one correspondence. Returned item has transpose data shape, so that main iteration is over bands.
    - **Parameters**
        - kapath   : `export_vasprun`().kpath or `get_kpts`().kpath.
        - evals_set: `export_vasprun`().bands.evals or `get_evals`().evals. If calculations are spin-polarized, it will be `...evals.SpinUp/SpinDown` for both. You need to apply twice for SpinUp and SpinDown separately.
        - pros_set : `export_vasprun().pro_bands.pros` or `get_bands_pro_set`().pros. If calculations are spin-polarized, it will be `...pros.SpinUp/SpinDown` for both. You need to create collections twice for SpinUp and SpinDown separately.
        - elements : List three lists of ions to project on, each element could be `range(start,stop,step)` as well, remember that `stop` is not included in python. so `range(0,2)` will generate 0 and 1 indices.
        - orbs     : List of three lists of orbitals indices. `[[red],[green],[blue]]`, you can create any color by this combination. For example, to get `s-orbital in yellow color`, you will use `[[0],[0],[]]`. Do not remove empty list from there, it will not effect your orbital selection.
        - interp_nk   : Dictionary with keys 'n' and 'k' for interpolation.
        - occs_set: `export_vasprun`().bands.occs or `get_evals`().occs. If calculations are spin-polarized, it will be `...occs.SpinUp/SpinDown` for both. You need to apply twice for SpinUp and SpinDown separately.
    - **Returns**
        - kpath    : List of path of NKPTS. (interpolated if given.)
        - evals    : An (NBAND,NKPTS) numpy arry.
        - occs     : An (NBAND,NKPTS) numpy arry.
        - colors   : An (NBANDS,NKPTS,3) numpy array.
        - widths   : An (NBAND,NKPTS) numpy arry, its actually colors summed along z-axis.
        - norms    : An (NBAND,NKPTS) numpy arry, normalized as 100%.
        - kpoints  : List of NKPTS. (interpolated if given.)
    """
    if pros_set == []:
        raise ValueError("Can not plot an empty eigenvalues object\n"
                         "Try with large energy range.")
    if len(orbs) < 3 :
        raise ValueError("orbs have structure [[],[],[]], do not reduce structure even if it is empty.")
    elif len(elements) <3:
        raise ValueError("elements have structure [[],[],[]], do not reduce structure even if it is empty.")
    else:
       r = np.take(pros_set,orbs[0],axis=3).sum(axis=3)
       r = np.take(r,list(elements[0]),axis=0).sum(axis=0)
       g = np.take(pros_set,orbs[1],axis=3).sum(axis=3)
       g = np.take(g,list(elements[1]),axis=0).sum(axis=0)
       b = np.take(pros_set,orbs[2],axis=3).sum(axis=3)
       b = np.take(b,list(elements[2]),axis=0).sum(axis=0)
       if interp_nk:
           from pivotpy import utils as gu
           knew,evals = gu.interpolate_data(kpath,evals_set,**interp_nk)
           _,occs = gu.interpolate_data(kpath,occs_set,**interp_nk)
           _, kpoints = gu.interpolate_data(kpath,kpoints,**interp_nk)
           r  = gu.interpolate_data(kpath,r,**interp_nk)[1].clip(min=0)
           g  = gu.interpolate_data(kpath,g,**interp_nk)[1].clip(min=0)
           b  = gu.interpolate_data(kpath,b,**interp_nk)[1].clip(min=0)
       else:
           knew,evals, occs, kpoints = kpath,evals_set, occs_set, kpoints

       evals = np.transpose(evals)
       occs = np.transpose(occs)
       max_c = max(max(map(max,r[:,:])),max(map(max,g[:,:])),max(map(max,b[:,:])))
       if max_c == 0:
           max_c=1 # Avoid divide by 0.
       _cl = np.concatenate((r,g,b)).reshape((3,-1,np.shape(r)[1])).swapaxes(0,2)
       rgb = _cl/max_c # Normalized overall data
       norms = (rgb*100).astype(int) # Normalized overall data
       widths = 0.0001+ np.sum(rgb,axis=2) #should be before scale colors

       # Now scale colors to 1 at each point.
       cl_max=np.max(_cl,axis=2)
       cl_max[cl_max==0.0] = 1 # avoid divide by zero. Contributions are 4 digits only.
       rgb = _cl/cl_max[:,:,np.newaxis] # Normalized per point

       return knew,evals,occs, rgb, widths, norms, kpoints

# Cell
def _flip_even_patches(array_1d, patch_length):
    """
    - When you reshape bands data to 1D array, you may need to draw lines which do not link ends of plot, for that, it is required to flip patches, so that next band start from where 1st end and so one.
    - **Parameters**
        - array_1d     : Numpy 1d array or list.
        - patch_length : length of xaxis patches, e.g NKPTS.
    - **Returns**
        - 1D list
    - ** Example**
        > k=[1,2,3,1,2,3]
        > _flip_even_patches(k,3)
        > [1,2,3,3,2,1]
    """
    out_put = []
    n= patch_length
    for i in range(len(array_1d)):
        if i % 2 != 0:
            out_put.extend(np.flip(array_1d[i*n:(i+1)*n]))
        else:
            out_put.extend(array_1d[i*n:(i+1)*n])
    return out_put

# Cell
def _rgb2plotly(rgb_data=None,mode='markers',max_width=None,showlegend=False,name='',labels=['s','p','d'],symbol=0,bands_indices = None):
    """
    - Returns data object of plotly's figure using `_get_rgb_data`. Returned data could be fed to a plolty's figure.
    - ** Parameters**
        - rgb_data    : output of `_get_rgb_data`.
        - mode        : Three plotting modes are available:
            - 'markers' : Plot whole data as a single scatter object. Its too fast.
            - 'bands'   : Plot data such that each band is accessible via legend.
            - 'lines'   : A replica of `matplotlib LineCollection` object. It plots at each point separately, slower than other two modes.
        - max_width  : Line/Scatter thickness is scaled to `max_width`. None by default and represent actual data.
        - name       : Name to be shown on hover text or legend.
        - labels     : Optional, show red green blue colors corresponding orbitals.
        - showlegend : Optional, only suitbale if spin up/down or 'bands' mode is ON.
        - symbol     : Plotly's marker symbol. 0 for circle, 5/6 for Up/Down.
        - bands_indices : `len(bands_indices) == NBANDS` should hold. `export_vasprun().bands.indices` or `get_evals().indices`.
    """
    if mode not in ('markers','bands','lines'):
        raise TypeError("Argument `mode` expects one of ['markers','bands','lines'], got '{}'.".format(mode))
    if rgb_data:
        k,en,occ, rgb,lws, norms, kpoints = rgb_data
        _indices = range(np.shape(en)[1]) if bands_indices is None else bands_indices

        _names = [["<sub>{}</sub>".format(int(1+i)) for j in range(len(en[0]))]
                  for i in _indices]
        clrs=(255*rgb).astype(int).clip(min=0,max=255) #clip in color range
        colors=[['rgb({},{},{})'.format(*i) for i in _c] for _c in clrs]
        txts=[['Projection: [{0}]</br>Norm (%): [{1}, {2}, {3}]'.format(', '.join(labels),*_n) for _n in _norm] for _norm in norms]

        if max_width:
            lws = max_width*lws/np.max(lws)

        ktext = [f'{x:>7.3f}{y:>7.3f}{z:>7.3f}' for (x,y,z) in kpoints]
        h_template = "</br>{} </br></br>Band: {}{}  Occ: {:>7.4f}</br>K<sub>{}</sub>: {}"
        h_text=[[h_template.format(t,name,l,o,i, ktext[i]) for i,(t,l,o) in enumerate(zip(ts,ns,os))]
                for ts,ns,os in zip(txts,_names,occ)]
        data=[]
        if mode == 'lines':
            lc = [[[k[i:i+2],ke[i:i+2]] for i in range(len(ke)-1)] for ke in en]
            lws = 0.5*(lws[:,:-1] + lws[:,1:])
            # Zip auto picks colors,lws,h_text in same number as en
            data = [
             go.Scatter(x=pt[0],y=pt[1],name=name, mode="lines",line=dict(color=c, width=w),
             showlegend=showlegend,hovertext=h)
                for pts,cs,ws,hs in zip(lc,colors,lws,h_text)
                    for pt,c,w,h in zip(pts,cs,ws,hs)]

        if mode == 'markers':
            nkpts = len(k) # should not move up or down.
            k = [k for bands in range(len(en))]
            k,en,colors,lws,h_text = [np.reshape(item,(-1)) for item in [k,en,colors,lws,h_text]]
            # Make lines reverse at even line.
            k,en,colors,lws,h_text =[_flip_even_patches(arr,nkpts) for arr in [k,en,colors,lws,h_text]]
            data.append(go.Scatter(x=k,y=en,mode='markers', name=name,
                        marker=dict(color=colors,size=lws,symbol=symbol),line=dict(width=0.001,
                        color='rgba(255,255,250,0)'),showlegend=showlegend,hovertext=h_text)
                        )
        if mode == 'bands':
            for i,e,c,w,t in zip(_indices,en,colors,lws,h_text):
                line_color = 'rgb({},{},{})'.format(*np.max(clrs[i - _indices[0]],axis=0))
                line_width = np.max(w)/8
                data.append(go.Scatter(x=k,y=e,mode='markers+lines',
                            name="{}<sub>{}</sub>".format(name,str(i+1)),
                            marker=dict(color=c,size=w,symbol=symbol),line_width= line_width,
                            line_color=line_color,showlegend=showlegend,hovertext=t)
                           )
        return data

# Cell
def iplot2html(fig,filename=None,out_string=False,modebar=True):
    """
    - Writes plotly's figure as HTML file or display in IPython which is accessible when online. It is different than plotly's `fig.to_html` as it is minimal in memory. If you need to have offline working file, just use `fig.write_html('file.html')` which will be larger in size.
    - **Parameters**
        - fig      : A plotly's figure object.
        - filename : Name of file to save fig. Defualt is None and show plot in Colab/Online or return hrml string.
        - out_string: If True, returns HTML string, if False displays graph if possible.
    """
    import uuid # Unique div-id required,otherwise jupyterlab renders at one place only and overwite it.
    div_id = "graph-{}".format(uuid.uuid1())
    fig_json = fig.to_json()
    # a simple HTML template
    if filename:
        _filename = gu.prevent_overwrite(filename)
        template = """<html>
        <head>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        </head>
        <body>
            <div id='{}'></div>
            <script>
                var fig_data = {}
                Plotly.react('{}', fig_data.data, fig_data.layout);
            </script>
        </body>
        </html>"""

        # write the JSON to the HTML template
        with open(_filename, 'w') as f:
            f.write(template.format(div_id,fig_json,div_id))

    else:
        if modebar==True: #Only for docs issue
            config = "{displayModeBar: true,scrollZoom: true}"
        else:
            config = "{displayModeBar: false,scrollZoom: true}"
        template = """<div>
        <script src='https://cdn.plot.ly/plotly-latest.min.js'></script>
            <div id='{}'><!-- Plotly chart DIV --></div>
            <script>
                var data = {};
                var config = {};
                Plotly.newPlot('{}', data.data,data.layout,config);
            </script>
        </div>""".format(div_id,fig_json,config,div_id)
        if out_string is True:
            return template
        else:
            from IPython.display import HTML
            return HTML(template)

# Cell
def iplot_rgb_lines(
    path_evr    = None,
    elements    = [[],[],[]],
    orbs        = [[],[],[]],
    labels      = ['','',''],
    mode        = 'markers',
    elim        = [],
    E_Fermi     = None,
    skipk       = None,
    kseg_inds  = [],
    max_width   = 6,
    title       = None,
    ktick_inds  = [0,-1],
    ktick_vals  = ['Γ','M'],
    figsize     = None,
    interp_nk   = {},
    query_data  = {}
    ):
    """
    - Returns plotly's figure object, takes care of spin-polarized calculations automatically. `elements`,`orbs` and `labels` are required to be one-to-one lists of size 3 where each item in list could be another list or integer.
    - **Parameters**
        - path_evr  : Path/to/vasprun.xml or xml output of `read_asxml`.
        - elements   : List of size 3 of list of indices of ions. If not given, picks all ions for each orbital.
        - orbs       : List of size 3 of list of orbital indices, if not gievn, s,p,d plotted.
        - labels  : List of labels for projection.
        - mode       : Three plotting modes are available:
            - 'markers' : Plot whole data as a single scatter object. Its too fast.
            - 'bands'   : Plot data such that each band is accessible via legend.
            - 'lines'   : A replica of `matplotlib LineCollection` object. It plots at each point separately, slower than other two modes.

        - interp_nk   : Dictionary with keys 'n' and 'k' for interpolation.
        - figsize   : Tuple(width,height) in pixels, e.g. (700,400).
        - query_data : Dictionary with keys as label and values as list of length 2. len(query_data) <=3 should hold for RGB plots. If given, used in place of elements, orbs and labels arguments.
                        Example: {'s':([0,1],[0]),'p':([0,1],[1,2,3]),'d':([0,1],[4,5,6,7,8])} will pick up s,p,d orbitals of first two ions of system.
        - **Other Parameters**
            - ktick_inds, ktick_vals,elim,kseg_inds,max_width,title etc.
    """
    if mode not in ('markers','bands','lines'):
        raise TypeError("Argument `mode` expects one of ['markers','bands','lines'], got '{}'.".format(mode))

    vr = vp._validate_evr(path_evr=path_evr,skipk=skipk,elim=elim)

    # Fix orbitals, elements and labels lengths very early.
    if query_data:
        elements, orbs, labels = sp._format_input(query_data,rgb=False) # prefer over elements, orbs and labels
    elements,orbs,labels = sp._validate_input(elements,orbs,labels,vr.sys_info,rgb=True)

    ## Main working here.
    if vr.pro_bands == None:
        raise ValueError("Can not plot an empty eigenvalues object.\n"
                         "Try with large energy range.")

    if E_Fermi == None:
        E_Fermi=vr.bands.E_Fermi
    K = vp.join_ksegments(vr.kpath,kseg_inds = kseg_inds)
    xticks=[K[i] for i in ktick_inds]
    xlim=[min(K),max(K)]
    if elim:
        ylim = [min(elim),max(elim)]
    else:
        ylim=[-10,10]
    #====Title Name======
    SYSTEM=vr.sys_info.SYSTEM
    if title == None:
        title= "{}[{}]".format(SYSTEM,','.join(labels))
    # After All Fixing
    ISPIN=vr.sys_info.ISPIN
    args_dict=dict(orbs=orbs,elements=elements,interp_nk=interp_nk) # Do not scale color there, scale here.
    data,showlegend,name=[],False,'' # Place holder

    if mode=='bands':
            showlegend=True
    kpoints = vr.kpoints
    if ISPIN == 1:
        En=vr.bands.evals-E_Fermi
        Oc = vr.bands.occs
        Pros=vr.pro_bands.pros
        new_args=dict(kpath=K, evals_set=En, pros_set=Pros, occs_set = Oc, **args_dict)
        rgb_lines=_get_rgb_data(kpoints = kpoints,**new_args)
        data=_rgb2plotly(rgb_data=rgb_lines,mode=mode,showlegend=showlegend,
                           labels=labels,name='B',max_width=max_width,bands_indices= vr.bands.indices)
    if ISPIN == 2:
        if mode == 'markers':
            showlegend=True
        En1=vr.bands.evals.SpinUp-E_Fermi
        En2=vr.bands.evals.SpinDown-E_Fermi
        Oc1 = vr.bands.occs.SpinUp
        Oc2 = vr.bands.occs.SpinDown
        Pros1=vr.pro_bands.pros.SpinUp
        Pros2=vr.pro_bands.pros.SpinDown
        new_args1=dict(kpath=K, evals_set=En1, pros_set=Pros1,occs_set = Oc1,**args_dict)
        rgb_lines1=_get_rgb_data(kpoints = kpoints,**new_args1)
        data1=_rgb2plotly(rgb_data=rgb_lines1,mode=mode,symbol=0,showlegend=showlegend,
                            labels=labels,name='B<sup>↑</sup>',max_width=max_width,bands_indices= vr.bands.indices)
        new_args2=dict(kpath=K, evals_set=En2, pros_set=Pros2,occs_set = Oc2,**args_dict)
        rgb_lines2=_get_rgb_data(kpoints = kpoints,**new_args2)
        data2=_rgb2plotly(rgb_data=rgb_lines2,mode=mode,symbol=100,showlegend=showlegend,
                            labels=labels,name='B<sup>↓</sup>',max_width=max_width,bands_indices= vr.bands.indices)
        data=[[d1,d2] for d1,d2 in zip(data1,data2)]
        data=[d for ds in data for d in ds]

    # Initiate figure
    fig=go.Figure(data=data)
    fig.update_layout(title=title,margin=go.layout.Margin(l=60,r=50,b=40,t=75,pad=0),
            yaxis=go.layout.YAxis(title_text='Energy (eV)',range=ylim),
            xaxis=go.layout.XAxis(ticktext=ktick_vals, tickvals=xticks,tickmode="array",range=xlim),
            font=dict(family="stix, serif",size=14))
    if(figsize!=None):
        fig.update_layout(width=figsize[0],height=figsize[1],autosize=False)
    #Draw lines at breakpoints
    if kseg_inds:
        kargs_dict = dict(mode='lines',line=dict(color='rgb(222,222,222)',width=1.2),showlegend=False)
        for pt in kseg_inds:
            fig.add_trace(go.Scatter(x=[K[pt],K[pt]],y=ylim,**kargs_dict))
            fig.add_trace(go.Scatter(x=[K[pt],K[pt]],y=ylim,**kargs_dict))
    update_args = dict(linewidth=0.1,linecolor='rgba(222,222,222,0.1)', mirror=True)
    fig.update_xaxes(showgrid=True,zeroline=False,showline=True,**update_args)
    fig.update_yaxes(showgrid=False,zeroline=True,showline=True,**update_args)
    return fig

# Cell
def iplot_dos_lines(
    path_evr      = None,
    elements      = [[0,],],
    orbs          = [[0],],
    labels        = ['s',],
    elim          = [],
    colormap      = 'gist_rainbow',
    tdos_color    = (0.5,0.95,0),
    linewidth     = 2,
    fill_area     = True,
    vertical      = False,
    E_Fermi       = None,
    figsize       = None,
    spin          = 'both',
    interp_nk     = {},
    title         = None,
    query_data    = {}
    ):
        """
        - Returns plotly's figure. If given,elements,orbs colors, and labels must have same length. If not given, zeroth ions is plotted with s-orbital.
        - **Parameters**)
            - path_evr   : Path/to/vasprun.xml or output of `export_vasprun`. Auto picks in CWD.
            - elements   : List [[0,],] of ions indices, by defualt plot first ion's projections.
            - orbs       : List [[0,],] lists of indices of orbitals, could be empty.
            - labels     : List [str,] of orbitals labels. len(labels) == len(orbs) must hold.
            - elim       : [min,max] of energy range.
            - E_Fermi    : If not given, automatically picked from `export_vasprun`.
            - colormap   : Matplotlib's standard color maps. Default is 'gist_ranibow'. Use 'RGB' if want to compare with `iplot_rgb_lines` with 3 projection inputs (len(orbs)==3).
            - fill_area  : Default is True and plots filled area for dos. If False, plots lines only.
            - vertical   : False, If True, plots along y-axis.
            - interp_nk  : Dictionary with keys 'n' and 'k' for interpolation.
            - figsize    : Tuple(width,height) in pixels, e.g. (700,400).
            - query_data : Dictionary with keys as label and values as list of length 2. If given, used in place of elements, orbs and labels arguments.
                        Example: {'s':([0,1],[0]),'p':([0,1],[1,2,3]),'d':([0,1],[4,5,6,7,8])} will pick up s,p,d orbitals of first two ions of system.
        - **Returns**
            - fig        : Plotly's figure object.
        """
        if query_data:
            elements,orbs,labels = sp._format_input(query_data,rgb=False) # prefer query_data over elements,orbs,labels
            # These are validated inisde _collect_dos, no need here
        en,tdos,pdos,vr=None,None,None,None # Place holders for defining
        cl_dos = sp._collect_dos(path_evr=path_evr,elim=elim,
                            elements=elements, orbs=orbs,labels=labels,
                            E_Fermi=E_Fermi, spin='both',interp_nk=interp_nk)
        try:
            en,tdos,pdos,labels,vr = cl_dos
        except:
            raise ValueError("Try with large energy range.")

        labels=[label.replace('$','').replace('^↑','<sup>↑</sup>').replace('^↓','<sup>↓</sup>') for label in labels]
        # Make additional colors for spin down. Inverted colors are better.
        if(elim):
            ylim=[min(elim),max(elim)]
        else:
            ylim=[-10,10]
        # Fix elements and colors length
        if colormap in plt.colormaps():
            from matplotlib.pyplot import cm
            if len(tdos) == 2:
                c_map   = cm.get_cmap(colormap)
                c_vals  = np.linspace(0,1,2*len(orbs))
                colors  = c_map(c_vals)
            else:
                c_map   = cm.get_cmap(colormap)
                c_vals  = np.linspace(0,1,len(orbs))
                colors  = c_map(c_vals)
            # Fix for RGB comparison
            if len(tdos) == 2 and 'both' in spin and len(orbs)==3:
                colors[[-1,-2]]= colors[[-2,-1]] #Flip last two colors only
        else:
            raise ValueError("`colormap` expects one of the follwoing:\n{}".format(plt.colormaps()))
        # Total DOS colors
        t_color=mpl.colors.to_rgb(tdos_color)
        it_color=gu.transform_color(t_color,c = -1) #inverts for c = -1
        #========Title Name========
        SYSTEM=vr.sys_info.SYSTEM
        if(title==None):
            title="{}".format(SYSTEM)

        fig = go.Figure()
        fig.update_layout(title=title,margin=go.layout.Margin(l=60,r=50,b=40,t=75,pad=0),\
                          font=dict(family="stix, serif",size=14))
        if(figsize!=None):
            fig.update_layout(width=figsize[0],height=figsize[1],autosize=False)
        if(vertical==False):
            if(fill_area==False):
                fill=None
            if(fill_area==True):
                fill='tozeroy'
            args_dic=dict(mode='lines',line_width=linewidth,fill=fill)
            fig.update_xaxes(range=ylim,title='Energy (eV)')
            if(len(tdos)==2):   # Spin polarized.
                fig.add_scatter(x=en,y=tdos[0],line_color='rgb({},{},{})'.format(*[int(255*i) for i in t_color]),\
                                 name='TDOS<sup>↑</sup>',**args_dic)
                fig.add_scatter(x=en,y=tdos[1],line_color='rgb({},{},{})'.format(*[int(255*i) for i in it_color]),\
                                 name='TDOS<sup>↓</sup>',**args_dic)
            else:   # unpolarized.
                fig.add_trace(go.Scatter(x=en,y=tdos,line_color='rgb({},{},{})'.format(*[int(255*i) for i in t_color]),\
                              name='TDOS',**args_dic))
            for p,l,c in zip(pdos,labels,colors):
                fig.add_trace(go.Scatter(x=en,y=p,line_color='rgb({},{},{})'.format(*[int(255*i) for i in c]),\
                               name=l,**args_dic))
        if(vertical==True):
            if(fill_area==False):
                fill=None
            if(fill_area==True):
                fill='tozerox'
            args_dic=dict(mode='lines',line_width=linewidth,fill=fill)
            fig.update_yaxes(range=ylim,title='Energy (eV)')
            if(len(tdos)==2):   # Spin polarized.
                fig.add_scatter(y=en,x=tdos[0],line_color='rgb({},{},{})'.format(*[int(255*i) for i in t_color]),\
                                name='TDOS<sup>↑</sup>',**args_dic)
                fig.add_scatter(y=en,x=tdos[1],line_color='rgb({},{},{})'.format(*[int(255*i) for i in it_color]),\
                                name='TDOS<sup>↓</sup>',**args_dic)
            else:   # unpolarized.
                fig.add_trace(go.Scatter(y=en,x=tdos,line_color='rgb({},{},{})'.format(*[int(255*i) for i in t_color]),\
                                name='TDOS',**args_dic))
            for p,l,c in zip(pdos,labels,colors):
                fig.add_trace(go.Scatter(y=en,x=p,line_color='rgb({},{},{})'.format(*[int(255*i) for i in c]),\
                                name=l,**args_dic))
        fig.update_xaxes(showgrid=True, zeroline=True,showline=True, linewidth=0.1, linecolor='rgba(222,222,222,0.1)', mirror=True)
        fig.update_yaxes(showgrid=True, zeroline=True,showline=True, linewidth=0.1, linecolor='rgba(222,222,222,0.1)', mirror=True)
        return fig