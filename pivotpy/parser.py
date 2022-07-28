import re
import os
from io import StringIO
from itertools import islice, chain, product
from collections import namedtuple
import textwrap
import xml.etree.ElementTree as ET

import numpy as np

from . import utils, serializer
#from pivotpy import utils, serializer


def dict2tuple(name,d):
    """Converts a dictionary (nested as well) to namedtuple, accessible via index and dot notation as well as by unpacking.
    - **Parameters**
        - name: Name of the tuple.
        - d   : Dictionary, nested works as well.
    """
    return namedtuple(name,d.keys())(
           *(dict2tuple(k.upper(),v) if isinstance(v,dict) else v for k,v in d.items())
           )

def read_asxml(path = None):
    """
    - Reads a big vasprun.xml file into memory once and then apply commands. If current folder contains `vasprun.xml` file, it automatically picks it.

    - **Parameters**
        - path : Path/To/vasprun.xml

    - **Returns**
        - xml_data : Xml object to use in other functions
    """
    path = path or './vasprun.xml'
    if not os.path.isfile(path):
        raise FileNotFoundError("File: '{}'' does not exist!".format(path))

    elif 'vasprun.xml' not in path:
        raise Exception("File name should end with 'vasprun.xml'")

    fsize = utils.get_file_size(path)
    value = float(fsize.split()[0])
    print_str = """
    Memory Consumption Warning!
    ---------------------------
    File: {} is large ({}). It may consume a lot of memory (generally 3 times the file size).
        An alternative way is to use `pivotpy.split_vasprun()` to split the file into multiple files and then read resulting '_vasprun.xml' file.
    """.format(path,fsize)
    if 'MB' in fsize and value > 200:
        print(utils.color.y(textwrap.dedent(print_str)))
    elif 'GB' in fsize and value > 1:
        print(utils.color.y(textwrap.dedent(print_str)))
        
    nt = namedtuple('VasprunXML',['path','root']) # Safe and full of info
    return nt(os.path.abspath(path), ET.parse(path).getroot()) # THis is xml_data for other functions

def xml2dict(xmlnode_or_filepath):
    """Convert xml node or xml file content to dictionary. All output text is in string format, so further processing is required to convert into data types/split etc.
    - The only paramenter `xmlnode_or_filepath` is either a path to an xml file or an `xml.etree.ElementTree.Element` object.
    - Each node has `tag,text,attr,nodes` attributes. Every text element can be accessed via
    `xml2dict()['nodes'][index]['nodes'][index]...` tree which makes it simple.
    """
    if isinstance(xmlnode_or_filepath,str):
        node = read_asxml(xmlnode_or_filepath)
    else:
        node = xmlnode_or_filepath

    text = node.text.strip() if node.text else ''
    nodes = [xml2dict(child) for child in list(node)]
    return {'tag': node.tag,'text': text, 'attr':node.attrib, 'nodes': nodes}

def exclude_kpts(xml_data):
    for kpts in xml_data.root.iter('varray'):
        if(kpts.attrib=={'name': 'weights'}):
            weights=[float(arr.text.strip()) for arr in kpts.iter('v')]
    exclude=[]
    [exclude.append(item) for item in weights if item!=weights[-1]];
    skipk=len(exclude) #that much to skip
    return skipk

def get_ispin(xml_data):
    for item in xml_data.root.iter('i'):
        if(item.attrib=={'type': 'int', 'name': 'ISPIN'}):
            return int(item.text)
        
def get_force(xml_data):
    "Reads force on each ion from vasprun.xml"
    forces = []
    for node in xml_data.root.iter('varray'):
        if 'name' in node.attrib and node.attrib['name'] == 'forces':
            for subnode in node.iter('v'):
                forces = [*forces, *subnode.text.split()]
    
    forces = np.fromiter(forces,dtype=float).reshape((-1,3))
    return forces

def get_scsteps(xml_data):
    "Reads scsteps energies from vasprun.xml"
    steps = []
    for node in xml_data.root.iter('scstep'):
        for e in node.iter('energy'):
            _d = {}
            for i in e.iter('i'):
                if '_en' in i.attrib['name']:
                    _d[i.attrib['name']] = float(i.text)
            steps.append(_d)
    if steps:
        arrays = {k:[] for k in steps[0].keys()}
        for step in steps:
            for k,v in step.items():
                arrays[k].append(v)
                
        return {k:np.array(v) for k,v in arrays.items()}


def get_space_info(xml_data):
    base_dir = os.path.split(os.path.abspath(xml_data.path))[0]
    path = os.path.join(base_dir,'OUTCAR')
    if not os.path.isfile(path):
        raise FileNotFoundError(f"File {path!r} not found. OUTCAR file should be in same directory as of {xml_data.path!r}.")
    
    info_dict = {}
    pos_line = islice2array(path,include = '^positions|^\s+positions',raw=True)
    info_dict['cartesian_positions'] = True if 'cartesian' in pos_line.lower() else False
    
    k_line = islice2array(path,include = '^k-points|^\s+k-points',raw=True).splitlines()[0]
    info_dict['cartesian_kpoints'] = True if 'cartesian' in k_line.lower() else False
    
    with open(path, 'r') as f:
        text = f.read()
        two_pi, reciprocal = re.findall(r'2pi/SCALE.*position',text, flags = re.DOTALL)[0].split('k-points')
        reciprocal = reciprocal.split('position')[0].strip().splitlines()[1:]
        two_pi = two_pi.strip().splitlines()[1:]
    
    k_rec = np.array([[float(v) for v in r.strip().split()[:3]] for r in reciprocal])
    k_scale = np.array([[float(v) for v in r.strip().split()[:3]] for r in two_pi])
    
    for final in xml_data.root.iter('structure'):
        if final.attrib.get('name','') == 'finalpos':
            for arr in final.iter('varray'):
                if arr.attrib.get('name','') == 'rec_basis':
                    rec_basis = np.array([[float(a) for a in v.text.split()] for v in arr.iter('v')])
                    k_coords = k_rec.dot(rec_basis)
                    ks_norm = np.linalg.norm(k_scale,axis=1)
                    kc_norm = np.linalg.norm(k_coords,axis=1)
                    where = (kc_norm > 1e-6) # Some kpoints may yield zero norm
                    info_dict['scale'] = (ks_norm[where]/kc_norm[where]).mean().round(12).astype(float)
    
    cartesian_kpath = np.cumsum(np.linalg.norm(k_scale[1:]-k_scale[:-1],axis=1))
    info_dict['cartesian_kpath'] = np.insert(cartesian_kpath,0,0)
    return serializer.Dict2Data(info_dict)
    

def get_summary(xml_data):
    for i_car in xml_data.root.iter('incar'):
        incar={car.attrib['name']:car.text.strip() for car in i_car}
    n_ions=[int(atom.text) for atom in xml_data.root.iter('atoms')][0]
    type_ions=[int(atom_types.text) for atom_types in xml_data.root.iter('types')][0]
    elem=[info[0].text.strip() for info in xml_data.root.iter('rc')]
    elem_name=[]; #collect IONS names
    [elem_name.append(item) for item in elem[:-type_ions] if item not in elem_name]
    elem_index=[0]; #start index
    [elem_index.append((int(entry)+elem_index[-1])) for entry in elem[-type_ions:]];
    ISPIN=get_ispin(xml_data=xml_data)
    NELECT = int([i.text.strip().split('.')[0] for i in xml_data.root.iter('i') if i.attrib['name']=='NELECT'][0])
    # Fields
    try:
        for pro in xml_data.root.iter('partial'):
            dos_fields=[field.text.strip() for field in pro.iter('field')]
            dos_fields = [field for field in dos_fields if 'energy' not in field]
    except:
        dos_fields = []
    for i in xml_data.root.iter('i'): #efermi for condition required.
        if(i.attrib=={'name': 'efermi'}):
            efermi=float(i.text)
    
    #Writing information to a dictionary
    space_info = get_space_info(xml_data=xml_data)
    info_dic={'SYSTEM':incar['SYSTEM'],'NION':n_ions,'NELECT':NELECT,'TypeION':type_ions,
              'ElemName':elem_name,'ElemIndex':elem_index,'Fermi': efermi,'ISPIN':ISPIN,
              'fields':dos_fields,'incar':incar, 'space_info':space_info}
    return serializer.Dict2Data(info_dic)

def join_ksegments(kpath,kseg_inds=[]):
    """Joins a broken kpath's next segment to previous. `kseg_inds` should be list of first index of next segment"""
    path_list = np.array(kpath)
    if kseg_inds:
        for ind in kseg_inds:
            path_list[ind:] -= path_list[ind] - path_list[ind-1]
    return list(path_list)

def get_kpts(xml_data, skipk = 0):
    for kpts in xml_data.root.iter('varray'):
        if(kpts.attrib=={'name': 'kpointlist'}):
            kpoints=[[float(item) for item in arr.text.split()] for arr in kpts.iter('v')]
    kpoints=np.array(kpoints[skipk:])
    #KPath solved.
    kpath = get_space_info(xml_data).cartesian_kpath[skipk:]
    kpath = kpath - kpath[0] # Shift to start at 0
    kpath = kpath/kpath[-1] # Normaliz to 1 to see all bandstuructres on common axis in full range
    # Do Not Join KPath if it is broken, leave that to plotting functions
    return serializer.Dict2Data({'NKPTS':len(kpoints),'kpoints':kpoints,'kpath':kpath})

def get_tdos(xml_data,spin_set=1,elim=[]):
    tdos=[]; #assign for safely exit if wrong spin set entered.
    ISPIN = get_ispin(xml_data=xml_data)
    for neighbor in xml_data.root.iter('dos'):
        for item in neighbor[1].iter('set'):
            if(ISPIN==1 and spin_set==1):
                if(item.attrib=={'comment': 'spin 1'}):
                    tdos=np.array([[float(entry) for entry in arr.text.split()] for arr in item])
            if(ISPIN==2 and spin_set==1):
                if(item.attrib=={'comment': 'spin 1'}):
                    tdos_1=np.array([[float(entry) for entry in arr.text.split()] for arr in item])
                if(item.attrib=={'comment': 'spin 2'}):
                    tdos_2=np.array([[float(entry) for entry in arr.text.split()] for arr in item])
                    tdos = {'SpinUp':tdos_1,'SpinDown':tdos_2}
            if(spin_set!=1): #can get any
                if(item.attrib=={'comment': 'spin {}'.format(spin_set)}):
                    tdos=np.array([[float(entry) for entry in arr.text.split()] for arr in item])
    for i in xml_data.root.iter('i'): #efermi for condition required.
        if(i.attrib=={'name': 'efermi'}):
            efermi=float(i.text)
    dos_dic= {'Fermi':efermi,'ISPIN':ISPIN,'tdos':tdos}
    #Filtering in energy range.
    if elim: #check if elim not empty
        if(ISPIN==1 and spin_set==1):
            up_ind=np.max(np.where(tdos[:,0]-efermi<=np.max(elim)))+1
            lo_ind=np.min(np.where(tdos[:,0]-efermi>=np.min(elim)))
            tdos=tdos[lo_ind:up_ind,:]
        if(ISPIN==2 and spin_set==1):
            up_ind=np.max(np.where(tdos['SpinUp'][:,0]-efermi<=np.max(elim)))+1
            lo_ind=np.min(np.where(tdos['SpinUp'][:,0]-efermi>=np.min(elim)))
            tdos = {'SpinUp':tdos_1[lo_ind:up_ind,:],'SpinDown':tdos_2[lo_ind:up_ind,:]}
        if(spin_set!=1):
            up_ind=np.max(np.where(tdos[:,0]-efermi<=np.max(elim)))+1
            lo_ind=np.min(np.where(tdos[:,0]-efermi>=np.min(elim)))
            tdos=tdos[lo_ind:up_ind,:]
        dos_dic= {'Fermi':efermi,'ISPIN':ISPIN,'grid_range':range(lo_ind,up_ind),'tdos':tdos}
    return serializer.Dict2Data(dos_dic)

def get_evals(xml_data, skipk = None, elim = []):
    evals, occs = [], [] #assign for safely exit if wrong spin set entered.
    ISPIN = get_ispin(xml_data=xml_data)
    if skipk != None:
        skipk=skipk
    else:
        skipk = exclude_kpts(xml_data=xml_data) #that much to skip by default
    for neighbor in xml_data.root.iter('eigenvalues'):
        for item in neighbor[0].iter('set'):
            if ISPIN == 1:
                if(item.attrib=={'comment': 'spin 1'}):
                    evals = np.array([[[float(t) for t in th.text.split()] for th in thing] for thing in item])[skipk:]
                    evals, occs = evals[:,:,0], evals[:,:,1]
                    NBANDS = len(evals[0])
            if ISPIN == 2:
                if(item.attrib=={'comment': 'spin 1'}):
                    eval_1 = np.array([[[float(t) for t in th.text.split()] for th in thing] for thing in item])[skipk:]
                    eval_1, occs_1 = eval_1[:,:,0], eval_1[:,:,1]
                if(item.attrib=={'comment': 'spin 2'}):
                    eval_2=np.array([[[float(t) for t in th.text.split()] for th in thing] for thing in item])[skipk:]
                    eval_2, occs_2 = eval_2[:,:,0], eval_2[:,:,1]
                    evals = {'SpinUp':eval_1,'SpinDown':eval_2}
                    occs = {'SpinUp':occs_1,'SpinDown':occs_2}
                    NBANDS=len(eval_1[0])

    for i in xml_data.root.iter('i'): #efermi for condition required.
        if i.attrib == {'name': 'efermi'}:
            efermi=float(i.text)
    evals_dic={'Fermi':efermi,'ISPIN':ISPIN,'NBANDS':NBANDS,'evals':evals,'indices': range(NBANDS),'occs':occs}
    if elim: #check if elim not empty
        if ISPIN == 1:
            up_ind = np.max(np.where(evals[:,:]-efermi <= np.max(elim))[1])+1
            lo_ind = np.min(np.where(evals[:,:]-efermi >= np.min(elim))[1])
            evals = evals[:,lo_ind:up_ind]
            occs = occs[:,lo_ind:up_ind]
        if ISPIN == 2:
            up_ind=np.max(np.where(eval_1[:,:]-efermi <= np.max(elim))[1])+1
            lo_ind=np.min(np.where(eval_1[:,:]-efermi >= np.min(elim))[1])
            evals={'SpinUp':eval_1[:,lo_ind:up_ind],'SpinDown':eval_2[:,lo_ind:up_ind]}
            occs = {'SpinUp':occs_1[:,lo_ind:up_ind],'SpinDown':occs_2[:,lo_ind:up_ind]}
        NBANDS = int(up_ind - lo_ind) #update Bands
        evals_dic['NBANDS'] = NBANDS
        evals_dic['indices'] = range(lo_ind,up_ind)
        evals_dic['evals'] = evals
        evals_dic['occs'] = occs
        
    return serializer.Dict2Data(evals_dic)

def get_bands_pro_set(xml_data, spin_set=1, skipk=0, bands_range=None, set_path=None):
    """Returns bands projection of a spin_set(default 1). If spin-polarized calculations, gives SpinUp and SpinDown keys as well.
    - **Parameters**
        - xml_data    : From `read_asxml` function
        - skipk       : Number of initil kpoints to skip (Default 0).
        - spin_set    : Spin set to get, default is 1.
        - bands_range : If elim used in `get_evals`,that will return bands_range to use here. Note that range(0,2) will give 2 bands 0,1 but tuple (0,2) will give 3 bands 0,1,2.
        - set_path    : path/to/_set[1,2,3,4].txt, works if `split_vasprun` is used before.
    - **Returns**
        - Data     : pivotpy.Dict2Data with attibutes of bands projections and related parameters.
    """
    if bands_range != None:
        check_list = list(bands_range)
        if check_list==[]:
            raise ValueError("No bands prjections found in given energy range.")
    # Try to read _set.txt first. instance check is important.
    if isinstance(set_path,str) and os.path.isfile(set_path):
        _header = islice2array(set_path,nlines=1,raw=True,exclude=None)
        _shape = [int(v) for v in _header.split('=')[1].strip().split(',')]
        NKPTS, NBANDS, NIONS, NORBS = _shape
        if NORBS == 3:
            fields = ['s','p','d']
        elif NORBS == 9:
            fields = ['s','py','pz','px','dxy','dyz','dz2','dxz','x2-y2']
        else:
            fields = [str(i) for i in range(NORBS)] #s,p,d in indices.
        COUNT = NIONS*NBANDS*(NKPTS-skipk)*NORBS
        start = NBANDS*NIONS*skipk
        nlines = None # Read till end.
        if bands_range:
            _b_r = list(bands_range)
            # First line is comment but it is taken out by exclude in islice2array.
            start = [[NIONS*NBANDS*k + NIONS*b for b in _b_r] for k in range(skipk,NKPTS)]
            start = [s for ss in start for s in ss] #flatten
            nlines = NIONS # 1 band has nions
            NBANDS = _b_r[-1]-_b_r[0]+1 # upadte after start

        NKPTS = NKPTS-skipk # Update after start, and bands_range.
        COUNT = NIONS*NBANDS*NKPTS*NORBS
        data = islice2array(set_path,start=start,nlines=nlines,count=COUNT)
        data = data.reshape((NKPTS,NBANDS,NIONS,NORBS)).transpose([2,0,1,3])
        return serializer.Dict2Data({'labels':fields,'pros':data})

    #Collect Projection fields
    fields=[];
    for pro in xml_data.root.iter('projected'):
        for arr in pro.iter('field'):
            if('eig' not in arr.text and 'occ' not in arr.text):
                fields.append(arr.text.strip())
    NORBS = len(fields)
    #Get NIONS for reshaping data
    NIONS=[int(atom.text) for atom in xml_data.root.iter('atoms')][0]

    for spin in xml_data.root.iter('set'):
        if spin.attrib=={'comment': 'spin{}'.format(spin_set)}:
            k_sets = [kp for kp in spin.iter('set') if 'kpoint' in kp.attrib['comment']]
    k_sets = k_sets[skipk:]
    NKPTS = len(k_sets)
    band_sets = []
    for k_s in k_sets:
        b_set = [b for b in k_s.iter('set') if 'band' in b.attrib['comment']]
        if bands_range == None:
            band_sets.extend(b_set)
        else:
            b_r = list(bands_range)
            band_sets.extend(b_set[b_r[0]:b_r[-1]+1])
    NBANDS = int(len(band_sets)/len(k_sets))
    try:
        # Error prone solution but 5 times fater than list comprehension.
        bands_pro = (float(t) for band in band_sets for l in band.iter('r') for t in l.text.split())
        COUNT = NKPTS*NBANDS*NORBS*NIONS # Must be counted for performance.
        data = np.fromiter(bands_pro,dtype=float,count=COUNT)
    except:
        # Alternate slow solution
        print("Error using `np.fromiter`.\nFalling back to (slow) list comprehension...",end=' ')
        bands_pro = (l.text for band in band_sets for l in band.iter('r'))
        bands_pro = [[float(t) for t in text.split()] for text in bands_pro]
        data = np.array(bands_pro)
        del bands_pro # Release memory
        print("Done.")

    data = data.reshape((NKPTS,NBANDS,NIONS,NORBS)).transpose((2,0,1,3))
    return serializer.Dict2Data({'labels':fields,'pros':data})

def get_dos_pro_set(xml_data,spin_set=1,dos_range=None):
    """Returns dos projection of a spin_set(default 1) as numpy array. If spin-polarized calculations, gives SpinUp and SpinDown keys as well.
    - **Parameters**
        - xml_data    : From `read_asxml` function
        - spin_set    : Spin set to get, default 1.
        - dos_range   : If elim used in `get_tdos`,that will return dos_range to use here..
    - **Returns**
        - Data     : pivotpy.Dict2Data with attibutes of dos projections and related parameters.
    """
    if dos_range != None:
        check_list = list(dos_range)
        if check_list == []:
            raise ValueError("No DOS prjections found in given energy range.")

    n_ions=get_summary(xml_data=xml_data).NION
    for pro in xml_data.root.iter('partial'):
        dos_fields=[field.text.strip()for field in pro.iter('field')]
        #Collecting projections.
        dos_pro=[]; set_pro=[]; #set_pro=[] in case spin set does not exists
        for ion in range(n_ions):
            for node in pro.iter('set'):
                if(node.attrib=={'comment': 'ion {}'.format(ion+1)}):
                    for spin in node.iter('set'):
                        if(spin.attrib=={'comment': 'spin {}'.format(spin_set)}):
                            set_pro=[[float(entry) for entry in r.text.split()] for r in spin.iter('r')]
            dos_pro.append(set_pro)
    if dos_range==None: #full grid computed.
        dos_pro=np.array(dos_pro) #shape(NION,e_grid,pro_fields)
    else:
        dos_range=list(dos_range)
        min_ind=dos_range[0]
        max_ind=dos_range[-1]+1
        dos_pro=np.array(dos_pro)[:,min_ind:max_ind,:]
    final_data=np.array(dos_pro) #shape(NION,e_grid,pro_fields)
    return serializer.Dict2Data({'labels':dos_fields,'pros':final_data})


def get_structure(xml_data):
    SYSTEM = [i.text for i in xml_data.root.iter('i') if i.attrib['name'] == 'SYSTEM'][0]

    for final in xml_data.root.iter('structure'):
        if(final.attrib=={'name': 'finalpos'}):
            for i in final.iter('i'):
                volume=float(i.text)
            for arr in final.iter('varray'):
                if(arr.attrib=={'name': 'basis'}):
                    basis=[[float(a) for a in v.text.split()] for v in arr.iter('v')]
                if(arr.attrib=={'name': 'rec_basis'}):
                    rec_basis=[[float(a) for a in v.text.split()] for v in arr.iter('v')]
                if(arr.attrib=={'name': 'positions'}):
                    positions=[[float(a) for a in v.text.split()] for v in arr.iter('v')]
    # element labels
    types  = [int(_type.text) for _type in xml_data.root.iter('types')][0]
    elems  = [info[0].text.strip() for info in xml_data.root.iter('rc')]
    _inds  = np.array([int(a) for a in elems[-types:]])
    _nums  = [k + 1 for i in _inds for k in range(i)]
    labels = [f"{e} {i}" for i, e in zip(_nums,elems)]

    INDS = np.cumsum([0,*_inds]).astype(int)
    Names = list(np.unique(elems[:-types]))
    unique_d = {e:range(INDS[i],INDS[i+1]) for i,e in enumerate(Names)}
    space_info = get_space_info(xml_data)
    scale = space_info.scale
    cartesian = space_info.cartesian_positions
    st_dic={'SYSTEM':SYSTEM,'volume': volume,'basis': np.array(basis),'rec_basis': np.array(rec_basis), 'scale': 1.0,
            'extra_info': {'comment':'Exported from vasprun.xml','cartesian':cartesian,'scale':scale},
            'positions': np.array(positions),'labels':labels,'unique': unique_d}
    return serializer.PoscarData(st_dic)

def export_vasprun(path = None, skipk = None, elim = [], dos_only = False):
    """
    - Returns a full dictionary of all objects from `vasprun.xml` file. It first try to load the data exported by powershell's `Export-VR(Vasprun)`, which is very fast for large files. It is recommended to export large files in powershell first.
    - **Parameters**
        - path       : Path to `vasprun.xml` file. Default is `'./vasprun.xml'`.
        - skipk      : Default is None. Automatically detects kpoints to skip.
        - elim       : List [min,max] of energy interval. Default is [], covers all bands.
        
    **Returns**: `pivotpy.serializer.VasprunData` object.
    """
    path = path or './vasprun.xml'

    xml_data = read_asxml(path=path)

    base_dir = os.path.split(os.path.abspath(path))[0]
    set_paths = [os.path.join(base_dir,"_set{}.txt".format(i)) for i in (1,2)]
    #First exclude unnecessary kpoints. Includes only same weight points
    if skipk!=None:
        skipk=skipk
    else:
        skipk = exclude_kpts(xml_data) #that much to skip by default
    info_dic = get_summary(xml_data) #Reads important information of system.
    #KPOINTS
    kpts = get_kpts(xml_data,skipk=skipk)
    #EIGENVALS
    eigenvals = get_evals(xml_data,skipk=skipk,elim=elim)
    #TDOS
    tot_dos = get_tdos(xml_data,spin_set=1,elim=elim)
    #Bands and DOS Projection
    if elim:
        bands_range = eigenvals.indices #indices in range form.
        grid_range=tot_dos.grid_range
    else:
        bands_range = None #projection function will read itself.
        grid_range = None
        
    if dos_only:
        bands_range = range(1) # Just one band
        skipk = len(kpts.kpath) + skipk - 2 # Just Single kpoint
        
    if info_dic.ISPIN == 1:
        pro_bands = get_bands_pro_set(xml_data=xml_data,spin_set=1,skipk=skipk,bands_range=bands_range,set_path=set_paths[0])
        pro_dos = get_dos_pro_set(xml_data=xml_data,spin_set=1,dos_range=grid_range)
    if info_dic.ISPIN == 2:
        pro_1 = get_bands_pro_set(xml_data=xml_data,spin_set=1,skipk=skipk,bands_range=bands_range,set_path=set_paths[0])
        pro_2 = get_bands_pro_set(xml_data=xml_data,spin_set=2,skipk=skipk,bands_range=bands_range,set_path=set_paths[1])
        pros={'SpinUp': pro_1.pros,'SpinDown': pro_2.pros}#accessing spins in dictionary after .pro.
        pro_bands={'labels':pro_1.labels,'pros': pros}
        pdos_1 = get_dos_pro_set(xml_data=xml_data,spin_set=1,dos_range=grid_range)
        pdos_2 = get_dos_pro_set(xml_data=xml_data,spin_set=1,dos_range=grid_range)
        pdos={'SpinUp': pdos_1.pros,'SpinDown': pdos_2.pros}#accessing spins in dictionary after .pro.
        pro_dos={'labels':pdos_1.labels,'pros': pdos}
    
    # Forces and steps
    force = get_force(xml_data)
    scsteps = get_scsteps(xml_data)

    #Structure
    poscar = get_structure(xml_data = xml_data)
    poscar = {'SYSTEM':info_dic.SYSTEM,**poscar.to_dict()}
    #Dimensions dictionary.
    dim_dic={'kpoints':'(NKPTS,3)','kpath':'(NKPTS,1)','bands':'⇅(NKPTS,NBANDS)','dos':'⇅(grid_size,3)','pro_dos':'⇅(NION,grid_size,en+pro_fields)','pro_bands':'⇅(NION,NKPTS,NBANDS,pro_fields)'}
    #Writing everything to be accessible via dot notation
    full_dic={'sys_info':info_dic,'dim_info':dim_dic,'kpoints':kpts.kpoints,'kpath':kpts.kpath,'bands':eigenvals,
             'tdos':tot_dos,'pro_bands':pro_bands,'pro_dos':pro_dos,'poscar': poscar,
             'force':force,'scsteps':scsteps}
    return serializer.VasprunData(full_dic)

def _validate_evr(path_evr=None,**kwargs):
    "Validates data given for plotting functions. Returns a tuple of (Boolean,data)."
    if type(path_evr) == serializer.VasprunData:
        return path_evr

    path_evr = path_evr or './vasprun.xml' # default path.

    if isinstance(path_evr,str):
        if os.path.isfile(path_evr):
            # kwargs -> skipk=skipk,elim=elim
            return export_vasprun(path=path_evr,**kwargs)
        else:
            raise FileNotFoundError(f'File {path_evr!r} not found!')
    # Other things are not valid.
    raise ValueError('path_evr must be a path string or output of export_vasprun function.')

def islice2array(path_or_islice,dtype=float,delimiter='\s+',
                include=None,exclude='#',raw=False,fix_format = True,
                start=0,nlines=None,count=-1,cols=None,new_shape=None
                ):
    """
    - Reads a sliced array from txt,csv type files and return to array. Also manages if columns lengths are not equal and return 1D array. It is faster than loading  whole file into memory. This single function could be used to parse EIGENVAL, PROCAR, DOCAR and similar files with just a combination of `exclude, include,start,stop,step` arguments.
    - **Parameters**
        - path_or_islice: Path/to/file or `itertools.islice(file_object)`. islice is interesting when you want to read different slices of an opened file and do not want to open it again and again. For reference on how to use it just execute `pivotpy.export_locpot??` in a notebook cell or ipython terminal to see how islice is used extensively.
        - dtype: float by default. Data type of output array, it is must have argument.
        - start,nlines: The indices of lines to start reading from and number of lines after start respectively. Only work if `path_or_islice` is a file path. both could be None or int, while start could be a list to read slices from file provided that nlines is int. The spacing between adjacent indices in start should be equal to or greater than nlines as pointer in file do not go back on its own.  These parameters are in output of `slice_data`
        > Note: `start` should count comments if `exclude` is None. You can use `slice_data` function to get a dictionary of `start,nlines, count, cols, new_shape` and unpack in argument instead of thinking too much.
        - count: `np.size(output_array) = nrows x ncols`, if it is known before execution, performance is increased. This parameter is in output of `slice_data`.
        - delimiter:  Default is `\s+`. Could be any kind of delimiter valid in numpy and in the file.
        - cols: List of indices of columns to pick. Useful when reading a file like PROCAR which e.g. has text and numbers inline. This parameter is in output of `slice_data`.
        - include: Default is None and includes everything. String of patterns separated by | to keep, could be a regular expression.
        - exclude: Default is '#' to remove comments. String of patterns separated by | to drop,could be a regular expression.
        - raw    : Default is False, if True, returns list of raw strings. Useful to select `cols`.
        - fix_format: Default is True, it sepearates numbers with poor formatting like 1.000-2.000 to 1.000 2.000 which is useful in PROCAR. Keep it False if want to read string literally.
        - new_shape : Tuple of shape Default is None. Will try to reshape in this shape, if fails fallbacks to 2D or 1D. This parameter is in output of `slice_data`.
    - **Examples**
        > `islice2array('path/to/PROCAR',start=3,include='k-point',cols=[3,4,5])[:2]`
        > array([[ 0.125,  0.125,  0.125],
        >        [ 0.375,  0.125,  0.125]])
        > `islice2array('path/to/EIGENVAL',start=7,exclude='E',cols=[1,2])[:2]`
        > array([[-11.476913,   1.      ],
        >        [  0.283532,   1.      ]])
    > Note: Slicing a dimension to 100% of its data is faster than let say 80% for inner dimensions, so if you have to slice more than 50% of an inner dimension, then just load full data and slice after it.
    """
    if nlines is None and isinstance(start,(list,np.ndarray)):
        print("`nlines = None` with `start = array/list` is useless combination.")
        return np.array([]) # return empty array.

    def _fixing(_islice,include=include, exclude=exclude,fix_format=fix_format,nlines=nlines,start=start):
        if include:
            _islice = (l for l in _islice if re.search(include,l))

        if exclude:
            _islice = (l for l in _islice if not re.search(exclude,l))

        _islice = (l.strip() for l in _islice) # remove whitespace and new lines

        # Make slices here after comment excluding.
        if isinstance(nlines,int) and isinstance(start,(list,np.ndarray)):
            #As islice moves the pointer as it reads, start[1:]-nlines-1
            # This confirms spacing between two indices in start >= nlines
            start = [start[0],*[s2-s1-nlines for s1,s2 in zip(start,start[1:])]]
            _islice = chain(*(islice(_islice,s,s+nlines) for s in start))
        elif isinstance(nlines,int) and isinstance(start,int):
            _islice = islice(_islice,start,start+nlines)
        elif nlines is None and isinstance(start,int):
            _islice = islice(_islice,start,None)

        # Negative connected digits to avoid, especially in PROCAR
        if fix_format:
            _islice = (re.sub(r"(\d)-(\d)",r"\1 -\2",l) for l in _islice)
        return _islice

    def _gen(_islice,cols=cols):
        for line in _islice:
            line = line.strip().replace(delimiter,'  ').split()
            if line and cols is not None: # if is must here.
                line = [line[i] for i in cols]
            for chars in line:
                yield dtype(chars)

    #Process Now
    if isinstance(path_or_islice,str) and os.path.isfile(path_or_islice):
        with open(path_or_islice,'r') as f:
            _islice = islice(f,0,None) # Read full, Will fix later.
            _islice = _fixing(_islice)
            if raw:
                return '\n'.join(_islice)
            # Must to consume islice when file is open
            data = np.fromiter(_gen(_islice),dtype=dtype,count=count)
    else:
        _islice = _fixing(path_or_islice)
        if raw:
            return '\n'.join(_islice)
        data = np.fromiter(_gen(_islice),dtype=dtype,count=count)

    if new_shape:
        try: data = data.reshape(new_shape)
        except: pass
    elif cols: #Otherwise single array.
        try: data = data.reshape((-1,len(cols)))
        except: pass
    return data

def slice_data(dim_inds,old_shape):
    """
    - Returns a dictionary that can be unpacked in arguments of isclice2array function. This function works only for regular txt/csv/tsv data files which have rectangular data written.
    - **Parameters**
        - dim_inds : List of indices array or range to pick from each dimension. Inner dimensions are more towards right. Last itmes in dim_inds is considered to be columns. If you want to include all values in a dimension, you can put -1 in that dimension. Note that negative indexing does not work in file readig, -1 is s special case to fetch all items.
        - old_shape: Shape of data set including the columns length in right most place.
    - **Example**
        - You have data as 3D arry where third dimension is along column.
        > 0 0
        > 0 2
        > 1 0
        > 1 2
        - To pick [[0,2], [1,2]], you need to give
        > slice_data(dim_inds = [[0,1],[1],-1], old_shape=(2,2,2))
        > {'start': array([1, 3]), 'nlines': 1, 'count': 2}
        - Unpack above dictionary in `islice2array` and you will get output array.
    - Note that dimensions are packed from right to left, like 0,2 is repeating in 2nd column.
    """
    # Columns are treated diffiernetly.
    if dim_inds[-1] == -1:
        cols = None
    else:
        cols = list(dim_inds[-1])

    r_shape = old_shape[:-1]
    dim_inds = dim_inds[:-1]
    for i,ind in enumerate(dim_inds.copy()):
        if ind == -1:
            dim_inds[i] = range(r_shape[i])
    nlines = 1
    #start = [[NIONS*NBANDS*k + NIONS*b for b in _b_r] for k in range(skipk,NKPTS)] #kind of thing.
    _prod_ = product(*dim_inds)
    _mult_ = [np.product(r_shape[i+1:]) for i in range(len(r_shape))]
    _out_ = np.array([np.dot(p,_mult_) for p in _prod_]).astype(int)
    # check if innermost dimensions could be chunked.
    step = 1
    for i in range(-1,-len(dim_inds),-1):
        _inds = np.array(dim_inds[i]) #innermost
        if np.max(_inds[1:] - _inds[:-1]) == 1: # consecutive
            step = len(_inds)
            _out_ = _out_[::step] # Pick first indices
            nlines = step*nlines
            # Now check if all indices picked then make chunks in outer dimensions too.
            if step != r_shape[i]: # Can't make chunk of outer dimension if inner is not 100% picked.
                break # Stop more chunking
    new_shape = [len(inds) for inds in dim_inds] #dim_inds are only in rows.
    new_shape.append(old_shape[-1])
    return {'start':_out_,'nlines':nlines,'count': nlines*len(_out_),'cols':cols,'new_shape':tuple(new_shape)}

def split_vasprun(path=None):
    """
    - Splits a given vasprun.xml file into a smaller _vasprun.xml file plus _set[1,2,3,4].txt files which contain projected data for each spin set.
    - **Parameters**
        - path: path/to/vasprun.xml file.
    - **Output**
        - _vasprun.xml file with projected data.
        - _set1.txt for projected data of colinear calculation.
        - _set1.txt for spin up data and _set2.txt for spin-polarized case.
        - _set[1,2,3,4].txt for each spin set of non-colinear calculations.
    """
    if not path:
        path = './vasprun.xml'
    if not os.path.isfile(path):
        raise FileNotFoundError("{!r} does not exist!".format(path))
    base_dir = os.path.split(os.path.abspath(path))[0]
    out_file = os.path.join(base_dir,'_vasprun.xml')
    out_sets = [os.path.join(base_dir,'_set{}.txt'.format(i)) for i in range(1,5)]
    # process
    with open(path,'r') as f:
        lines = islice(f,None)
        indices = [i for i,l in enumerate(lines) if re.search('projected|/eigenvalues',l)]
        f.seek(0)
        print("Writing {!r} ...".format(out_file),end=' ')
        with open(out_file,'w') as outf:
            outf.write(''.join(islice(f,0,indices[1])))
            f.seek(0)
            outf.write(''.join(islice(f,indices[-1]+1,None)))
            print('Done')

        f.seek(0)
        middle = islice(f,indices[-2]+1,indices[-1]) #projected words excluded
        spin_inds = [i for i,l in enumerate(middle) if re.search('spin',l)][1:] #first useless.
        if len(spin_inds)>1:
            set_length = spin_inds[1]-spin_inds[0] # Must define
        else:
            set_length = indices[-1]-indices[-2] #It is technically more than set length, but fine for 1 set
        f.seek(0) # Must be at zero
        N_sets = len(spin_inds)
        # Let's read shape from out_file as well.
        xml_data = read_asxml(out_file)
        _summary = get_summary(xml_data)
        NIONS  = _summary.NION
        NORBS  = len(_summary.fields)
        NBANDS = get_evals(xml_data).NBANDS
        NKPTS  = get_kpts(xml_data).NKPTS
        del xml_data # free meory now.
        for i in range(N_sets): #Reads every set
            print("Writing {!r} ...".format(out_sets[i]),end=' ')
            start = (indices[-2]+1+spin_inds[0] if i==0 else 0) # pointer is there next time.
            stop_ = start + set_length # Should move up to set length only.
            with open(out_sets[i],'w') as setf:
                setf.write("  # Set: {} Shape: (NKPTS[NBANDS[NIONS]],NORBS) = {},{},{},{}\n".format(i+1,NKPTS,NBANDS,NIONS,NORBS))
                middle = islice(f,start,stop_)
                setf.write(''.join(l.lstrip().replace('/','').replace('<r>','') for l in middle if '</r>' in l))
                print('Done')

def export_spin_data(path = None, spins = 's', skipk = None, elim = None):
    """
    - Returns Data with selected spin sets. For spin polarized calculations, it returns spin up and down data.
    - **Parameters**
        - path   : Path to `vasprun.xml` file. Default is `'./vasprun.xml'`.
        - skipk  : Default is None. Automatically detects kpoints to skip.
        - elim   : List [min,max] of energy interval. Default is [], covers all bands.
        - spins  : Spin components to include from 'sxyz', e.g. 'sx' will pick <S> and <S_x> if present.
                   Only works if ISPIN == 1, otherwise it will be two sets for spin up and down.
    
    **Returns**: `pivotpy.serializer.SpinData` object.
    """
    if not isinstance(spins,str):
        raise TypeError(f"`spins` must be a string from 'sxyz', got {spins}!")

    if False in [comp in 'sxyz' for comp in spins]:
        raise ValueError(f"`spins` must be in 'sxyz', got {spins!r}!")

    xml_data = read_asxml(path = path or './vasprun.xml')

    base_dir = os.path.split(os.path.abspath(path or './vasprun.xml'))[0]
    set_paths = [os.path.join(base_dir,"_set{}.txt".format(i)) for i in (1,2,3,4)]

    skipk = skipk or exclude_kpts(xml_data=xml_data) #that much to skip by default
    full_dic = {'sys_info':get_summary(xml_data)}

    ISPIN = full_dic['sys_info'].ISPIN
    LSORBIT = getattr(full_dic['sys_info'].incar, 'LSORBIT', 'FALSE')
    if 'f' in LSORBIT.lower() and ISPIN == 1:
        for comp in spins:
            if comp in 'xyz':
                raise ValueError(f"LSORBIT = {LSORBIT} does not include spin component {comp!r}!")

    full_dic['dim_info'] = {'kpoints':'(NKPTS,3)','evals.<e,u,d>':'⇅(NKPTS,NBANDS)','spins.<u,d,s,x,y,z>':'⇅(NION,NKPTS,NBANDS,pro_fields)'}
    full_dic['kpoints']= get_kpts(xml_data, skipk = skipk).kpoints

    bands = get_evals(xml_data, skipk = skipk,elim = elim).to_dict()
    evals = bands['evals']
    bands.update({'u': evals['SpinUp'], 'd': evals['SpinDown']} if ISPIN == 2 else {'e': evals})
    
    del bands['evals'] # Do not Delete occupancies here
    full_dic['evals'] = bands

    bands_range = full_dic['bands'].indices if elim else None #indices in range form.

    spin_sets = {}
    if ISPIN == 1:
        for n, s in enumerate('sxyz', start = 1):
            if s in spins:
                spin_sets[s] = get_bands_pro_set(xml_data, spin_set = n, skipk = skipk, bands_range = bands_range, set_path = set_paths[n-1]).pros

    if ISPIN == 2:
        print(utils.color.g(f"Found ISPIN = 2, output data got attributes spins.<u,d> instead of spins.<{','.join(spins)}>"))
        pro_1 = get_bands_pro_set(xml_data, spin_set = 1, skipk = skipk, bands_range = bands_range, set_path = set_paths[0])
        pro_2 = get_bands_pro_set(xml_data, spin_set = 2, skipk = skipk, bands_range = bands_range, set_path = set_paths[1])
        spin_sets = {'u': pro_1.pros,'d': pro_2.pros}

    full_dic['spins'] = spin_sets
    full_dic['spins']['labels'] = full_dic['sys_info'].fields
    full_dic['poscar'] = {'SYSTEM':full_dic['sys_info'].SYSTEM,**(get_structure(xml_data).to_dict())}
    return serializer.SpinData(full_dic)

def export_outcar(path=None):
    "Read potential at ionic sites from OUTCAR file."
    if path is None:
        path = './OUTCAR'
    if not os.path.isfile(path):
        raise FileNotFoundError("{} does not exist!".format(path))
    # Raeding it
    with open(r'{}'.format(path),'r') as f:
        lines = f.readlines()
    # Processing
    for i,l in enumerate(lines):
        if 'NIONS' in l:
            N = int(l.split()[-1])
            nlines = np.ceil(N/5).astype(int)
        if 'electrostatic' in l:
            start_index = i+3
            stop_index = start_index+nlines
        if 'fractional' in l:
            first = i+1
        if 'vectors are now' in l:
            b_first = i+5
        if 'NION' in l:
            ion_line = l
        if 'NKPTS' in l:
            kpt_line =l

    NKPTS,NKDIMS,NBANDS = [int(v) for v in re.findall(r"\d+",kpt_line)]
    NEDOS,NIONS = [int(v) for v in re.findall(r"\d+",ion_line)]
    n_kbi = (NKPTS,NBANDS,NIONS)
    # Data manipulation
    # Potential
    data = lines[start_index:stop_index]
    initial = np.loadtxt(StringIO(''.join(data[:-1]))).reshape((-1))
    last = np.loadtxt(StringIO(data[-1]))
    pot_arr = np.hstack([initial,last]).reshape((-1,2))
    pot_arr[:,0] = pot_arr[:,0]-1 # Ion index fixing
    # Nearest neighbors
    pos = lines[first:first+N]
    pos_arr = np.loadtxt(StringIO('\n'.join(pos)))
    pos_arr[pos_arr>0.98] = pos_arr[pos_arr>0.98]-1 # Fixing outer layers
    # positions and potential
    pos_pot = np.hstack([pos_arr,pot_arr[:,1:]])
    basis = np.loadtxt(StringIO(''.join(lines[b_first:b_first+3])))
    final_dict = {'ion_pot':pot_arr,'positions':pos_arr,'site_pot':pos_pot,'basis':basis[:,:3],'rec_basis':basis[:,3:],'n_kbi':n_kbi}
    return serializer.OutcarData(final_dict)

def export_locpot(locpot=None,e = True,m = False):
    """
    - Returns Data from LOCPOT and similar structure files like CHG. Loads only single set out of 2/4 magnetization data to avoid performance/memory cost while can load electrostatic and one set of magnetization together.
    - **Parameters**
        - locpot: path/to/LOCPOT or similar stuructured file like CHG. LOCPOT is auto picked in CWD.
        - e     : Electric potential/charge density. Default is True.
        - m     : Magnetization density m. Default is False. If True, picks `m` for spin polarized case, and `m_x` for non-colinear case. Additionally it can take 'x','y' and 'z' in case of non-colinear calculations.
    - **Exceptions**
        - Would raise index error if magnetization density set is not present in LOCPOT/CHG in case `m` is not False.
    """
    if locpot is None:
        if os.path.isfile('LOCPOT'):
            locpot = 'LOCPOT'
        else:
            raise FileNotFoundError('./LOCPOT not found.')
    else:
        if not os.path.isfile(locpot):
            raise FileNotFoundError("File {!r} does not exist!".format(locpot))
    if m not in [True,False,'x','y','z']:
        raise ValueError("m expects one of [True,False,'x','y','z'], got {}".format(e))
    # data fixing after reading islice from file.
    def fix_data(islice_gen,shape):
        new_gen = (float(l) for line in islice_gen for l in line.split())
        COUNT = np.prod(shape).astype(int)
        data = np.fromiter(new_gen,dtype=float,count=COUNT) # Count is must for performance
        # data written on LOCPOT is in shape of (NGz,NGy,NGx)
        N_reshape = [shape[2],shape[1],shape[0]]
        data = data.reshape(N_reshape).transpose([2,1,0])
        return data
    # Reading File
    with open(locpot,'r') as f:
        lines = []
        f.seek(0)
        for i in range(8):
            lines.append(f.readline())
        N = sum([int(v) for v in lines[6].split()])
        f.seek(0)
        poscar = []
        for i in range(N+8):
            poscar.append(f.readline())
        f.readline() # Empty one
        Nxyz = [int(v) for v in f.readline().split()] # Grid line read
        nlines = np.ceil(np.prod(Nxyz)/5).astype(int)
        #islice is faster generator for reading potential
        pot_dict = {}
        if e == True:
            pot_dict.update({'e':fix_data(islice(f, nlines),Nxyz)})
            ignore_set = 0 # Pointer already ahead.
        else:
            ignore_set = nlines # Needs to move pointer to magnetization
        #reading Magnetization if True
        ignore_n = np.ceil(N/5).astype(int)+1 #Some kind of useless data
        if m == True:
            print("m = True would pick m_x for non-colinear case, and m for ISPIN=2.\nUse m='x' for non-colinear or keep in mind that m will refer to m_x.")
            start = ignore_n+ignore_set
            pot_dict.update({'m': fix_data(islice(f, start,start+nlines),Nxyz)})
        elif m == 'x':
            start = ignore_n+ignore_set
            pot_dict.update({'m_x': fix_data(islice(f, start,start+nlines),Nxyz)})
        elif m == 'y':
            start = 2*ignore_n+nlines+ignore_set
            pot_dict.update({'m_y': fix_data(islice(f, start,start+nlines),Nxyz)})
        elif m == 'z':
            start = 3*ignore_n+2*nlines+ignore_set
            pot_dict.update({'m_z': fix_data(islice(f, start,start+nlines),Nxyz)})

    # Read Info
    from .sio import export_poscar # Keep inside to avoid import loop
    poscar_data = export_poscar(text_plain = '\n'.join(p.strip() for p in poscar))
    basis = poscar_data.basis
    system = poscar_data.SYSTEM
    ElemName = list(poscar_data.unique.keys())
    ElemIndex = [0,*[v[-1]+ 1 for v in poscar_data.unique.values()]]
    positions = poscar_data.positions

    final_dict = dict(SYSTEM=system,ElemName=ElemName,ElemIndex=ElemIndex,basis=basis,positions=positions)
    final_dict = {**final_dict,**pot_dict}
    return serializer.LocpotData(final_dict)
