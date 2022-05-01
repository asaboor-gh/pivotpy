from . import __version__, __doc__, docs
from .sio import str2kpath, InvokeMaterialsProject, _save_mp_API
import argparse
from argparse import RawTextHelpFormatter
def main():
    print(f'Pivotpy\n=======\nVersion: {__version__}')
    print(__doc__)
    print('Loading Online DOCS...')
    docs()

    
def get_kpath():
    parser = argparse.ArgumentParser(description='Process KPATH string.', formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('kpath_str', help='''kpath multiline string e.g.:\n\t0 0 0 !$\Gamma$ 5\n\t1/2 1/2 1/2 !L 3\n\t0 0 1/2 !A\n\
will add 5 kpoints in the interval 1 and 3 in interval 2.\
Empty lines are taken as breaks in path.''')
    parser.add_argument('-w', '--weight',type=int, default=None)
    parser.add_argument('-n', '--n',type=int, default=10, help='Number of kpoints per unit length.')
    parser.add_argument('-z','--ibzkpt', type=str, help='IBZKPT file path, useful for HSE06 calculations.')
    
    args = parser.parse_args()
    
    return str2kpath(args.kpath_str,weight=args.weight,n=args.n, ibzkpt=args.ibzkpt)

def poscar():
    parser = argparse.ArgumentParser(description='Download POSCAR from Materials Project Website.')
    
    parser.add_argument('formula',type=str, help='Formula like "GaAs"')
    parser.add_argument('-i','--mp_id',type=str, help='ID of crystal as on Materials Project Website (optional).')
    parser.add_argument('-M', '--max_sites',type=int, help='Number of max sites in POSCAR, do not need if --mp_id given.')
    parser.add_argument('-m', '--min_sites',type=int, help='Number of min sites in POSCAR, do not need if --mp_id given.')
    parser.add_argument('-k','--api_key',type=str, help= 'Materials project API key. Save it using `pivotpy_save_mp_key` just once.')
    
    args = parser.parse_args()
    mp = InvokeMaterialsProject(api_key=args.api_key)
    output = mp.request(args.formula,mp_id=args.mp_id,max_sites=args.max_sites, min_sites=args.min_sites)
    
    if output:
        for car in output:
            print("\ncar\n{('='*66}")
            print(car.export_poscar().text_plain)
    else:
        print('No POSCAR found with provided input, try increasing range')
        
def api_key_save():
    parser = argparse.ArgumentParser(description='Saves API key from Materials Project website on your local machine for every time use.')
    parser.add_argument('api_key',type=str, help='Get API key from Materilas Project website and enter here.')
    args = parser.parse_args()
    return _save_mp_API(api_key=args.api_key)
