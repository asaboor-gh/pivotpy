# AUTOGENERATED! DO NOT EDIT! File to edit: InteractivePlots.ipynb (unless otherwise specified).

__all__ = ['get_rgb_data', 'flip_even_patches', 'rgb_to_plotly', 'plotly_to_html', 'plotly_rgb_lines',
           'plotly_dos_lines']

# Cell
def get_rgb_data(   kpath       = None,
                    evals_set   = None,
                    pros_set    = None,
                    elements    = [[0],[],[]],
                    orbs        = [[0],[],[]],
                    interpolate = False,
                    n           = 5,
                    k           = 3,
                    scale_color = False
    ):
    """
    - Returns a formatted RGB colored data to pass into `rgb_to_plotly` function. Two arguments, `elements` and `orbs` should be in one-to-one correspondence. Returned item has transpose data shape, so that main iteration is over bands.
    - **Parameters**
        - ax       : Matplotlib axes object, if not given, linecollection is returned.
        - kapath   : `export_vasprun`().kpath or `get_kpts`().kpath.
        - evals_set: `export_vasprun`().bands.evals or `get_evals`().evals. If calculations are spin-polarized, it will be `...evals.SpinUp/SpinDown` for both. You need to apply twice for SpinUp and SpinDown separately.
        - pros_set : `export_vasprun().pro_bands.pros` or `get_bands_pro_set`().pros. If calculations are spin-polarized, it will be `...pros.SpinUp/SpinDown` for both. You need to create collections twice for SpinUp and SpinDown separately.
        - elements : List three lists of ions to project on, each element could be `range(start,stop,step)` as well, remember that `stop` is not included in python. so `range(0,2)` will generate 0 and 1 indices.
        - orbs     : List of three lists of orbitals indices. `[[red],[green],[blue]]`, you can create any color by this combination. For example, to get `s-orbital in yellow color`, you will use `[[0],[0],[]]`. Do not remove empty list from there, it will not effect your orbital selection.
        - interpolate: Deafult is false, if True, it will add n points between nearest kpoints.
        - n        : int, default is 5. Adds n points between nearest kpoints.
        - k        : int, order of interpolation, defualt is 3. `n > k` should be hold.
        - scale_color: If True, colors are scaled to 1 at each points.
    - **Returns**
        - kpath    : List of NKPTS. (interpolated if given.)
        - evals    : An (NBAND,NKPTS) numpy arry.
        - colors   : An (NBANDS,NKPTS,3) numpy array.
    """
    if(pros_set==[]):
        import pivotpy.g_utils as gu
        gu.printy("Can not plot an empty eigenvalues object.")
        return gu.printg("Try with large energy range.")
    if len(orbs)<3 :
        raise ValueError("orbs have structure [[],[],[]], do not reduce structure even if it is empty.")
    elif len(elements) <3:
        raise ValueError("elements have structure [[],[],[]], do not reduce structure even if it is empty.")
    else:
       import numpy as np
       r = np.take(pros_set,orbs[0],axis=3).sum(axis=3)
       r = np.take(r,list(elements[0]),axis=0).sum(axis=0)
       g = np.take(pros_set,orbs[1],axis=3).sum(axis=3)
       g = np.take(g,list(elements[1]),axis=0).sum(axis=0)
       b = np.take(pros_set,orbs[2],axis=3).sum(axis=3)
       b = np.take(b,list(elements[2]),axis=0).sum(axis=0)
       if(interpolate==True):
           from pivotpy import g_utils as gu
           knew,evals = gu.interpolate_data(kpath,evals_set,n=n,k=k)
           r  = gu.interpolate_data(kpath,r,n=n,k=k)[1].clip(min=0)
           g  = gu.interpolate_data(kpath,g,n=n,k=k)[1].clip(min=0)
           b  = gu.interpolate_data(kpath,b,n=n,k=k)[1].clip(min=0)
       else:
           knew,evals = kpath,evals_set
       evals = np.transpose(evals)
       max_c = max(max(map(max,r[:,:])),max(map(max,g[:,:])),max(map(max,b[:,:])))
       if(max_c==0):
           max_c=1 # Avoid divide by 0.
       _cl = np.concatenate((r,g,b)).reshape((3,-1,np.shape(r)[1])).swapaxes(0,2)
       rgb = _cl/max_c # Normalized overall
       widths = 0.0001+ np.sum(rgb,axis=2) #should be before scale colors
       if(scale_color==True):
           cl_max=np.max(_cl,axis=2)
           # avoid divide by zero. Contributions are 4 digits only.
           cl_max[cl_max==0.0] = 1
           rgb = _cl/cl_max[:,:,np.newaxis] # Normalized per point
       return knew,evals,rgb, widths

# Cell
def flip_even_patches(array_1d, patch_length):
    """
    - When you reshape bands data to 1D array, you may need to draw lines which do not link ends of plot, for that, it is required to flip patches, so that next band start from where 1st end and so one.
    - **Parameters**
        - array_1d     : Numpy 1d array or list.
        - patch_length : length of xaxis patches, e.g NKPTS.
    - **Returns**
        - 1D list
    - ** Example**
        > k=[1,2,3,1,2,3]
        > flip_even_patches(k,3)
        > [1,2,3,3,2,1]
    """
    import numpy as np
    out_put = []
    n= patch_length
    for i in range(len(array_1d)):
        if i % 2 != 0:
            out_put.extend(np.flip(array_1d[i*n:(i+1)*n]))
        else:
            out_put.extend(array_1d[i*n:(i+1)*n])
    return out_put

# Cell
def rgb_to_plotly(rgb_data=None,mode='markers',max_width=5,showlegend=False,name='',labels=['s','p','d'],symbol=0):
    """
    - Returns data object of plotly's figure using `get_rgb_data`. Returned data could be fed to a plolty's figure.
    - ** Parameters**
        - rgb_data    : output of `get_rgb_data`.
        - mode        : Three plotting modes are available:
            - 'markers' : Plot whole data as a single scatter object. Its too fast.
            - 'bands'   : Plot data such that each band is accessible via legend.
            - 'lines'   : A replica of `matplotlib LineCollection` object. It plots at each point separately, slower than other two modes.
        - max_width  : Line/Scatter thickness is scaled to `max_width`.
        - name       : Name to be shown on hover text or legend.
        - labels     : Optional, show red green blue colors corresponding orbitals.
        - showlegend : Optional, only suitbale if spin up/down or 'bands' mode is ON.
        - symbol     : Plotly's marker symbol. 0 for circle, 5/6 for Up/Down.
    """
    if mode not in ('markers','bands','lines'):
        raise TypeError("Argument `mode` expects one of ['markers','bands','lines'], got '{}'.".format(mode))
        return
    if(rgb_data):
        import plotly.graph_objects as go
        import numpy as np
        k,en,rgb,lw = rgb_data

        _names = [["<sub>{}</sub>".format(int(1+i)) for j in range(len(en[0]))]
                  for i in range(len(en))]
        clrs=(255*rgb).astype(int)
        _txt=(100*rgb).astype(int)
        colors=[['rgb({},{},{})'.format(*i) for i in _c] for _c in clrs]
        txts=[['RGB({},{},{})'.format(*i) for i in _t] for _t in _txt]
        lws = max_width*lw
        h_template = "Proj: {} <br>{} </br>Band: {}{}"
        h_text=[[h_template.format([*labels],t,name,l) for t,l in zip(ts,ns)]
                for ts,ns in zip(txts,_names)]
        data=[]
        if(mode=='lines'):
            lc = [[[k[i:i+2],ke[i:i+2]] for i in range(len(ke)-1)] for ke in en]
            lws = 0.5*(lws[:,:-1] + lws[:,1:])
            # Zip auto picks colors,lws,h_text in same number as en
            data = [
             go.Scatter(x=pt[0],y=pt[1],name=name, mode="lines",line=dict(color=c, width=w),
             showlegend=showlegend,hovertext=h)
                for pts,cs,ws,hs in zip(lc,colors,lws,h_text)
                    for pt,c,w,h in zip(pts,cs,ws,hs)]

        if(mode=='markers'):
            nkpts = len(k) # should not move up or down.
            k = [k for bands in range(len(en))]
            k,en,colors,lws,h_text = [np.reshape(item,(-1)) for item in [k,en,colors,lws,h_text]]
            # Make lines reverse at even line.
            k,en,colors,lws,h_text =[flip_even_patches(arr,nkpts) for arr in [k,en,colors,lws,h_text]]
            data.append(go.Scatter(x=k,y=en,mode='markers', name=name,
                        marker=dict(color=colors,size=lws,symbol=symbol),line=dict(width=0.001,
                        color='rgba(255,255,250,0)'),showlegend=showlegend,hovertext=h_text)
                        )
        if(mode=='bands'):
            for i,values in enumerate(zip(en,colors,lws,h_text)):
                e,c,w,t = values
                data.append(go.Scatter(x=k,y=e,mode='markers+lines',
                            name="{}<sub>{}</sub>".format(name,str(i+1)),
                            marker=dict(color=c,size=w,symbol=symbol),line=dict(width=0.001,
                            color='rgba(255,255,250,0)'),showlegend=showlegend,hovertext=t)
                           )
        return data

# Cell
def plotly_to_html(fig,filename=None,out_string=False,modebar=True):
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
        with open(filename, 'w') as f:
            f.write(template.format(div_id,fig_json,div_id))
        f.close()
    else:
        if modebar==True: #Only for docs issue
            config = "{displayModeBar: true,editable: true,scrollZoom: true}"
        else:
            config = "{displayModeBar: false,editable: true,scrollZoom: true}"
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
def plotly_rgb_lines(path_evr    = None,
                    elements    = [[],[],[]],
                    orbs        = [[],[],[]],
                    labels      = ['','',''],
                    mode        = 'markers',
                    elim        = [],
                    E_Fermi     = None,
                    skipk       = None,
                    joinPathAt  = [],
                    max_width   = 6,
                    title       = None,
                    xt_indices  = [0,-1],
                    xt_labels   = ['Γ','M'],
                    figsize     = None,
                    interpolate = False,
                    n           = 5,
                    k           = 3

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
        - **kwargs      : interpolate, ticks, figsize,elim,joinPathAt,max_width,title etc.
    """
    import numpy as np
    import pivotpy.vr_parser as vp
    import pivotpy.s_plots as sp
    import pivotpy.i_plots as ip
    import plotly.graph_objects as go
    if mode not in ('markers','bands','lines'):
        raise TypeError("Argument `mode` expects one of ['markers','bands','lines'], got '{}'.".format(mode))
        return
    if(len(orbs) < 3 or len(elements) < 3):
        raise ValueError("orbs/elements have structure [[],[],[]], do not reduce structure even if it is empty.")
        return
    #checking type of given path.
    if(path_evr==None):
        vr=vp.export_vasprun(path=path_evr,skipk=skipk,elim=elim,joinPathAt=joinPathAt)
    if(path_evr!=None):
        from os import path as pt
        if(type(path_evr)==vp.Dict2Data):
            vr=path_evr
        elif(pt.isfile(path_evr)):
            vr=vp.export_vasprun(path=path_evr,skipk=skipk,elim=elim,joinPathAt=joinPathAt)
        else:
            return print("path_evr = `{}` does not exist".format(path_evr))
    # Apply a robust final check.
    try:
        vr.bands;vr.kpath
    except:
        return print("Object: \n{} \nis like a lower tree of export_vasprun(). Expects top tree.".format(vr))
    else:
        ## Main working here.
        if(vr.pro_bands==None):
            import pivotpy.g_utils as gu
            gu.printy("Can not plot an empty eigenvalues object.")
            return gu.printg("Try with large energy range.")
        #=====================================================
        orbs = [[item] if type(item)==int else item for item in orbs] #Fix if integer given.
        elem_inds = vr.sys_info.ElemIndex
        max_ind   = elem_inds[-1]-1 # Last index is used for range in ElemIndex, not python index.

        nfields=len(vr.sys_info.fields)

        # Fix int elements and orbs
        for i,e in enumerate(elements):
            if type(e)==int and e < elem_inds[-1]:
                elements[i] = range(elem_inds[e],elem_inds[e+1])
        _elements_inds = [e for es in elements for e in es]
        if _elements_inds and max(_elements_inds) > max_ind:
            return print("index {} is out of bound for {} elements.".format(max(_elements_inds),max_ind+1))
        _orb_inds = [p for orb in orbs for p in orb]
        if _orb_inds and max(_orb_inds) > nfields-1:
            return print("index {} is out of bound for {} orbitals.".format(max(_orb_inds),nfields))

        if(E_Fermi==None):
            E_Fermi=vr.bands.E_Fermi
        K=vr.kpath
        xticks=[K[i] for i in xt_indices]
        xlim=[min(K),max(K)]
        if(elim):
            ylim=[min(elim),max(elim)]
        else:
            ylim=[-10,10]
        # If elements not given, get whole system
        if _elements_inds==[]:
            elements = [range(0,max_ind+1),range(0,max_ind+1),range(0,max_ind+1)]
        # If orbs not given, get whole projections.
        if(_orb_inds==[]):
            if(nfields==3):
                orbs=[[0],[1],[2]]
            if(nfields==9 or nfields==16):
                orbs=[[0],[1,2,3],[4,5,6,7,8]]
        if _elements_inds==[] and _orb_inds == []:
            labels=['sys-s','sys-p','sys-d']
        #====Title Name======
        SYSTEM=vr.sys_info.SYSTEM
        if(title==None):
            title= "{}[{}]".format(SYSTEM,','.join(labels))

        # After All Fixing
        ISPIN=vr.sys_info.ISPIN
        args_dict=dict(orbs=orbs,elements=elements,interpolate=interpolate,n=n,k=k,scale_color=True) # Do not scale color there, scale here.
        data,showlegend,name=[],False,'' # Place holder
        if(mode=='bands'):
                showlegend=True
        if(ISPIN==1):
            En=vr.bands.evals-E_Fermi
            Pros=vr.pro_bands.pros
            new_args=dict(kpath=K, evals_set=En, pros_set=Pros,**args_dict)
            rgb_lines=get_rgb_data(**new_args)
            data=rgb_to_plotly(rgb_data=rgb_lines,mode=mode,showlegend=showlegend,labels=labels,name='B',max_width=max_width)
        if(ISPIN==2):
            if(mode=='markers'):
                showlegend=True
            En1=vr.bands.evals.SpinUp-E_Fermi
            En2=vr.bands.evals.SpinDown-E_Fermi
            Pros1=vr.pro_bands.pros.SpinUp
            Pros2=vr.pro_bands.pros.SpinDown
            new_args1=dict(kpath=K, evals_set=En1, pros_set=Pros1,**args_dict)
            rgb_lines1=get_rgb_data(**new_args1)
            data1=rgb_to_plotly(rgb_data=rgb_lines1,mode=mode,symbol=0,showlegend=showlegend,labels=labels,name='B<sup>↑</sup>',max_width=max_width)
            new_args2=dict(kpath=K, evals_set=En2, pros_set=Pros2,**args_dict)
            rgb_lines2=get_rgb_data(**new_args2)
            data2=rgb_to_plotly(rgb_data=rgb_lines2,mode=mode,symbol=100,showlegend=showlegend,labels=labels,name='B<sup>↓</sup>',max_width=max_width)
            data=[[d1,d2] for d1,d2 in zip(data1,data2)]
            data=[d for ds in data for d in ds]
        # Initiate figure
        fig=go.Figure(data=data)
        fig.update_layout(title=title,
            margin=go.layout.Margin(l=60,r=50,b=40,t=75,pad=0),#paper_bgcolor="whitesmoke",
            yaxis=go.layout.YAxis(title_text='Energy (eV)',range=ylim),
            xaxis=go.layout.XAxis(ticktext=xt_labels, tickvals=xticks,
            tickmode="array",range=xlim),font=dict(family="stix, serif",size=14))
        if(figsize!=None):
            fig.update_layout(width=figsize[0],height=figsize[1],autosize=False)
        #Draw lines at breakpoints
        if(joinPathAt):
            for pt in joinPathAt:
                fig.add_trace(go.Scatter(x=[K[pt],K[pt]],y=ylim,mode='lines',line=dict(color='rgb(0,0,0)',width=2),showlegend=False))
                fig.add_trace(go.Scatter(x=[K[pt],K[pt]],y=ylim,mode='lines',line=dict(color='rgb(222,222,222)',width=1.2),showlegend=False))
        fig.update_xaxes(showgrid=True, zeroline=False,showline=True, linewidth=0.1, linecolor='rgba(222,222,222,0.1)', mirror=True)
        fig.update_yaxes(showgrid=False, zeroline=True,showline=True, linewidth=0.1, linecolor='rgba(222,222,222,0.1)', mirror=True)
        return fig

# Cell
def plotly_dos_lines(path_evr     = None,
                    elim          = [],
                    elements      = [[0,],],
                    orbs          = [[0],],
                    labels        = ['s',],
                    color_map     = 'gist_rainbow',
                    tdos_color    = (0.5,0.95,0),
                    linewidth     = 2,
                    fill_area     = True,
                    vertical      = False,
                    E_Fermi       = None,
                    figsize       = None,
                    spin          = 'both',
                    interpolate   = False,
                    n             = 5,
                    k             = 3,
                    title         = None,
                    ):
        """
        - Returns ax object (if ax!=False) and plot on which all matplotlib allowed actions could be performed, returns lists of energy,tdos and pdos and labels. If given,elements,orbs colors, and labels must have same length. If not given, zeroth ions is plotted with s-orbital.
        - **Parameters**)
            - path_evr   : Path/to/vasprun.xml or output of `export_vasprun`. Auto picks in CWD.
            - elim       : [min,max] of energy range.
            - E_Fermi    : If not given, automatically picked from `export_vasprun`.
            - elements   : List [[0,],] of ions indices, by defualt plot first ion's projections.
            - orbs       : List [[0,],] lists of indices of orbitals, could be empty.
            - labels     : List [str,] of orbitals labels. len(labels)==len(orbs) must hold.
            - color_map  : Matplotlib's standard color maps. Default is 'gist_ranibow'. Use 'RGB' if want to compare with `plotly_rgb_lines` with 3 projection inputs (len(orbs)==3).
            - fill_area  : Default is True and plots filled area for dos. If False, plots lines only.
            - vertical   : False, If True, plots along y-axis.
            - figsize    : Tuple in pixels (width,height).
            - interpolate: Default is False, if True, bands are interpolated.
            - n          : int, number of points, default is 5.
            - k          : int, order of interpolation 0,1,2,3. Defualt 3. `n > k` should be hold.
            - legend_kwargs: Dictionary to contain legend arguments to fix.
        - **Returns**
            - fig        : Plotly's figure object.
        """
        import pivotpy.s_plots as sp
        import matplotlib.pyplot as plt
        import numpy as np
        import pivotpy.g_utils as gu
        import matplotlib as mpl
        import plotly.graph_objects as go

        en,tdos,pdos,vr=None,None,None,None # Place holders for defining
        cl_dos=sp.collect_dos(path_evr=path_evr,elim=elim, elements=elements, orbs=orbs,\
                              labels=labels, E_Fermi=E_Fermi, spin='both', interpolate=interpolate, n=n, k=k)
        try:
            en,tdos,pdos,labels,vr = cl_dos
        except TypeError:
            import pivotpy.g_utils as gu
            return gu.printg("Try with large energy range.")

        labels=[label.replace('$','').replace('^↑','<sup>↑</sup>').replace('^↓','<sup>↓</sup>') for label in labels]
        # Make additional colors for spin down. Inverted colors are better.
        if(elim):
            ylim=[min(elim),max(elim)]
        else:
            ylim=[-10,10]
        # Fix elements and colors length
        if color_map in plt.colormaps():
            from matplotlib.pyplot import cm
            if len(tdos) == 2:
                c_map   = cm.get_cmap(color_map)
                c_vals  = np.linspace(0,1,2*len(orbs))
                colors  = c_map(c_vals)
            else:
                c_map   = cm.get_cmap(color_map)
                c_vals  = np.linspace(0,1,len(orbs))
                colors  = c_map(c_vals)
            # Fix for RGB comparison
            if len(tdos) == 2 and 'both' in spin and len(orbs)==3:
                colors[[-1,-2]]= colors[[-2,-1]] #Flip last two colors only
        else:
            return print("`color_map` expects one of the follwoing:\n{}".format(plt.colormaps()))
        # Total DOS colors
        t_color=mpl.colors.to_rgb(tdos_color)
        it_color=gu.invert_color(color=t_color)
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