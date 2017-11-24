import numpy as np
import pandas as pd
import glob

from vald_to_pkl import vald_convert


class Elements(object):
    
    Z2X = ['', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
           'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 
           'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
           
    X2Z = {}
    for i, X in enumerate(Z2X):
        X2Z[X] = i


def create_input(file):
    """
    Crée le fichier d'input à partir du fichier .pkl de l'ion.
    """
    
    ion = file[file.rindex('\\')+1:file.rindex('.')].split('_')
    
    df = pd.read_pickle(file)
    columns = ['wav1', 'loggf', 'el', 'eu', 'rad', 'stark', 'waals']
    format = ['%12.5f', '% 7.3f'] +  2*['%8.4f'] +3*['% 6.2f']
    header = ' Wav (A)    log gf    E (l, u) (eV)    rad    stark  waals'
    
    # print(df[columns].shape)
    
    Z = Elements.X2Z[ion[0]]
    name = 'Z%02dE%02d' % (Z, Z-int(ion[1])+1)
    
    np.savetxt('input/'+name+'.dat', df[columns].values, fmt = format, header = header)


if __name__ == '__main__':

    for file in glob.glob('RAW/*'):
        Elem = file[file.rindex('\\')+1:file.rindex('.')]
        vald_convert(file, Elem)
        
    for file in glob.glob('PKL/*'):
        create_input(file)