# PivotPy
> A Python Processing Tool for Vasp Input/Output. A CLI is available in Powershell, see <a href='https://github.com/massgh/Vasp2Visual'>Vasp2Visual</a>.


[![Run in Azure](https://notebooks.azure.com/launch.png)](https://testazurenotebooks-massaz.notebooks.azure.com/j/notebooks/test.ipynb)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/massgh/pivotpy/blob/master/test.ipynb)




<style>a{text-decoration: none !important;color:lightkblue;font-weight:bold;}
                a:focus,a:active,a:hover{color:hotpink !important;}</style>
> [&nbsp;`▶` Index●&nbsp;](https://massgh.github.io/pivotpy/)  
> [&nbsp;`▶` Example&nbsp;](https://massgh.github.io/pivotpy/Example)  
> [&nbsp;`▶` StaticPlots&nbsp;](https://massgh.github.io/pivotpy/StaticPlots)  
> [&nbsp;`▶` InteractivePlots&nbsp;](https://massgh.github.io/pivotpy/InteractivePlots)  
> [&nbsp;`▶` SpinProjectedSurfaces&nbsp;](https://massgh.github.io/pivotpy/SpinProjectedSurfaces)  
> [&nbsp;`▶` StructureIO&nbsp;](https://massgh.github.io/pivotpy/StructureIO)  
> [&nbsp;`▶` Widgets&nbsp;](https://massgh.github.io/pivotpy/Widgets)  
> [&nbsp;`▶` MainAPI&nbsp;](https://massgh.github.io/pivotpy/MainAPI)  




## Install
`pip install pivotpy`

## How to use
- Use commnad `pivotpy` in regular terminal to quickly launch documentation any time. 
- See [Full Documentation](https://massgh.github.io/pivotpy/).
- See a [Calculation Example](https://massgh.github.io/pivotpy/Example.html)
- For CLI, use [Vasp2Visual](https://github.com/massgh/Vasp2Visual).
- See [PDF Slides](https://github.com/massgh/InteractiveHTMLs/tree/master/docs/IPySlides.pdf) for detailed introduction.

## CLI commnds
- Use `pivotpy` in system terminal to launch DOCS. 
- Use `pivotpy_get_poscar` to download POSCAR.
- Use `pivotpy_get_kpath` to create fine controlled KPATH. 

## Plot in Terminal without GUI
Use `pp.plt2text(colorful=True/False)` after matplotlib's code and your figure will appear in terminal. You need to zoom out alot to get a good view like below.

Tip: Use file [matplotlib2terminal.py](https://gist.github.com/massgh/d5cc44ad32510d3ff58cfefd75c6884e) on github independent of this package to plot in terminal. 
![IMG](terminal.jpg)

# Ipywidgets-based GUI
See GIF here:
![GIF](widget.gif) 

# Live Slides in Jupyter Notebook
Navigate to  [ipyslides](https://github.com/massgh/ipyslides) or do `pip install ipyslides` to create beautiful data driven presentation in Jupyter Notebook.

```python
import os, pivotpy as pp
with pp.set_dir('E:/Research/graphene_example/ISPIN_1/bands'):
    vr = pp.Vasprun(elim=[-5,5])

print('Try Follwing Methods:')
for v in dir(vr):
    if not v.startswith('_'):
        print('vr.'+v)
```

    Try Follwing Methods:
    vr.Fermi
    vr.data
    vr.elim
    vr.fermi
    vr.get_band_info
    vr.get_en_diff
    vr.get_fermi
    vr.iplot_dos_lines
    vr.iplot_rgb_lines
    vr.kticks
    vr.poscar
    vr.select
    vr.splot_bands
    vr.splot_color_lines
    vr.splot_dos_lines
    vr.splot_en_diff
    vr.splot_rgb_lines
    vr.to_json
    vr.to_pickle
    

```python
import matplotlib.pyplot as plt 
ax1,ax2 = pp.get_axes((6,3),ncols=2)
ax1.plot(vr.data.scsteps['e_fr_energy'],lw=3, label = 'e_fr_energy',color='k')
ax1.plot(vr.data.scsteps['e_0_energy'],lw=0.7,ls='dashed',label='e_0_energy',color='skyblue')
ax1.set_ylabel('Energy (eV)')
ax1.set_xlabel('Iteration Number')
ax1.legend()

vr.poscar.splot_lat(ax=ax2,plane='xy')
X, Y, Z = vr.data.poscar.coords.T
q = ax2.quiver(X,Y,*vr.data.force[:,:2].T,scale=25,color='r')
ax2.quiverkey(q, 0.7, 1, 7, 'Force (arb. units)')
ax2.add_legend()
```


![svg](docs/images/output_8_0.svg)


## Matplotlib's static plots
- Add anything from legend,colorbar, colorwheel. In below figure, all three are shown.
- Use aliases such as sbands, sdos,srgb,irgb,scolor,idos for plotting. 

```python
#collapse_input
import pivotpy as pp, numpy as np
import matplotlib.pyplot as plt 
vr1=pp.Vasprun('E:/Research/graphene_example/ISPIN_2/bands/vasprun.xml')
vr2=pp.Vasprun('E:/Research/graphene_example/ISPIN_2/dos/vasprun.xml')
axs = pp.get_axes(ncols=3,widths=[2,1,2.2],sharey=True,wspace=0.05,figsize=(8,2.6))
elements=[0,[0],[0,1]]
orbs=[[0],[2],[1,3]]
labels=['s','$p_z$','$(p_x+p_y)$']
ti_cks=dict(ktick_inds=[0,30,60,-1],ktick_vals=['Γ','M','K','Γ'])
args_dict=dict(elements=elements,orbs=orbs,labels=labels,elim=[-20,15],colormap='viridis',)
vr1.splot_bands(ax=axs[0],**ti_cks,elim=[-20,15])
vr1.splot_rgb_lines(ax=axs[2],**args_dict,**ti_cks,colorbar=False)
vr2.splot_dos_lines(ax=axs[1],vertical=True,spin='both',include_dos='pdos',**args_dict,legend_kwargs={'ncol': 3})
axs[2].color_cube(loc=(0.7,0.25),size=0.35)
pp._show() 
```

    [0;92m Given 0 at position 1 of sequence => 'C': range(0, 2). To just pick one ion, write it as [0].[00m
    


![svg](docs/images/output_10_1.svg)


## Interactive plots using plotly

```python
args_dict['labels'] = ['s','p_z','p_x+p_y']
args_dict.pop('colormap')
fig1 = vr1.iplot_rgb_lines(**args_dict)
#pp.iplot2html(fig1) #Do inside Google Colab, fig1 inside Jupyter
from IPython.display import Markdown
Markdown("[See Interactive Plot](https://massgh.github.io/InteractiveHTMLs/iGraphene.html)")
```




[See Interactive Plot](https://massgh.github.io/InteractiveHTMLs/iGraphene.html)



## Brillouin Zone (BZ) Processing
- Look in `pivotpy.sio` module or `pivotpy.api.POSCAR` class for details on generating mesh and path of KPOINTS as well as using Materials Projects' API to get POSCAR right in the working folder. Below is a screenshot of interactive BZ plot. You can `double click` on blue points and hit `Ctrl + C` to copy the high symmetry points relative to reciprocal lattice basis vectors. 
- Same color points lie on a sphere, with radius decreasing as red to blue and  gamma point in gold color. These color help distinguishing points but the points not always be equivalent, for example in FCC, there are two points on mid of edges connecting square-hexagon and hexagon-hexagon at equal distance from center but not the same points. 
- Any colored point's hover text is in gold background.      
#### Look the output of `pivotpy.sio.splot_bz`.
![BZ](docs\images\3bz.jpg)
![Energy](docs\images\energy.jpg)
[See Interactive BZ Plot](https://massgh.github.io/InteractiveHTMLs/BZ.html)

## Interpolation 
Amost every bandstructure and DOS plot function has an argument `interp_nk` which is a dictionary with keys `n` (Number of additional points between adjacent points) and `k` (order of interpolation 0-3). `n > k` must hold.

```python
#collapse_input
import pivotpy as pp, matplotlib.pyplot as plt
plt.style.use('ggplot')
k = vr1.data.kpath
ef = vr1.data.bands.Fermi
evals = vr1.data.bands.evals.SpinUp - ef
#Let's interpolate our graph to see effect. It is useful for colored graphs.
knew,enew=pp.interpolate_data(x=k,y=evals,n=10,k=3)
plot = plt.plot(k,evals,'m',lw=5,label='real data')
plot = plt.plot(k,evals,'w',lw=1,label='interpolated',ls='dashed')
pp.splots.add_text(ax=plt.gca(),txts='Graphene')
```


![svg](docs/images/output_15_0.svg)


## LOCPOT,CHG Visualization
check out the class `pivotpy.LOCPOT` to visulize local potential/charge and magnetization in a given direction.

## Running powershell commands from python.
Some tasks are very tideious in python while just a click way in powershell. See below, and try to list processes in python yourself to see the difference!

```python
pp.utils.ps2std(ps_command='(Get-Process)[0..4]')
```

    [32;1m NPM(K)    PM(M)      WS(M)     CPU(s)      Id  SI P[0m
    [32;1m                                                   r[0m
    [32;1m                                                   o[0m
    [32;1m                                                   c[0m
    [32;1m                                                   e[0m
    [32;1m                                                   s[0m
    [32;1m                                                   s[0m
    [32;1m                                                   N[0m
    [32;1m                                                   a[0m
    [32;1m                                                   m[0m
    [32;1m                                                   e[0m
    [32;1m ------    -----      -----     ------      --  -- -[0m
    22     6.89      13.18       0.39   11856   2 A
    6     1.39       5.23       0.00    7388   0 A
    19     8.45      22.41       0.00    5096   0 A
    27    26.51      47.55       9.78   17340   2 A
    10     1.80       6.54       0.00    5080   0 a
    

## Advancaed: Poweshell Cell/Line Magic `%%ps/%ps`
- You can create a IPython cell magic to run powershell commands directly in IPython Shell/Notebook (Powershell core installation required).
- Cell magic can be assigned to a variable `foo` by `%%ps --out foo`
- Line magic can be assigned to a variable by `foo = %ps powershell_command`

Put below code in ipython profile's startup file (create one) "~/.ipython/profile_default/startup/powershell_magic.py"
```python
from IPython.core.magic import register_line_cell_magic
from IPython import get_ipython
@register_line_cell_magic
def ps(line, cell=None):
    if cell:
        return get_ipython().run_cell_magic('powershell',line,cell)
    else:
        get_ipython().run_cell_magic('powershell','--out posh_output',line)
        return posh_output.splitlines()
``` 
Additionally you need to add following lines in "~/.ipython/profile_default/ipython_config.py" file to make above magic work.
```python
from traitlets.config.application import get_config
c = get_config()
c.ScriptMagics.script_magics = ['powershell']
c.ScriptMagics.script_paths = {
    'powershell' : 'powershell.exe -noprofile -command -',
    'pwsh': 'pwsh.exe -noprofile -command -'
}
```

```python
%%ps 
Get-ChildItem 'E:\Research\graphene_example\'
```

    
    
        Directory: E:\Research\graphene_example
    
    
    Mode                 LastWriteTime         Length Na
                                                      me
    ----                 -------------         ------ --
    da----          6/9/2022  10:33 AM                IS
                                                      PI
                                                      N_
                                                      1 
    da----          5/9/2020   1:05 PM                IS
                                                      PI
                                                      N_
                                                      2 
    -a----          5/9/2020   1:01 PM          75331 OU
                                                      TC
                                                      AR
    -a----          5/9/2020   1:01 PM         240755 va
                                                      sp
                                                      ru
                                                      n.
                                                      xm
                                                      l 
    
    


```python
x = %ps (Get-ChildItem 'E:\Research\graphene_example\').Name
x
```




    ['ISPIN_1', 'ISPIN_2', 'OUTCAR', 'vasprun.xml']



[Functions Reference](functions.md)




<style>a{text-decoration: none !important;color:lightkblue;font-weight:bold;}
                a:focus,a:active,a:hover{color:hotpink !important;}</style>
> [&nbsp;`▶` Index●&nbsp;](https://massgh.github.io/pivotpy/)  
> [&nbsp;`▶` Example&nbsp;](https://massgh.github.io/pivotpy/Example)  
> [&nbsp;`▶` StaticPlots&nbsp;](https://massgh.github.io/pivotpy/StaticPlots)  
> [&nbsp;`▶` InteractivePlots&nbsp;](https://massgh.github.io/pivotpy/InteractivePlots)  
> [&nbsp;`▶` SpinProjectedSurfaces&nbsp;](https://massgh.github.io/pivotpy/SpinProjectedSurfaces)  
> [&nbsp;`▶` StructureIO&nbsp;](https://massgh.github.io/pivotpy/StructureIO)  
> [&nbsp;`▶` Widgets&nbsp;](https://massgh.github.io/pivotpy/Widgets)  
> [&nbsp;`▶` MainAPI&nbsp;](https://massgh.github.io/pivotpy/MainAPI)  



