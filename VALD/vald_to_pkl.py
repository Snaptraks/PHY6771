import pandas as pd
import numpy as np


def vald_convert(file, Elem, wlim=None):
    """
    Fonction pour convertir les fichiers de VALD en pandas.DataFrame
    (puis les sauvegarder en .pkl pour plus tard)
    
    format:
    DataFrame: ['wav1','loggf','confl','terml','jl','el',
                'confu','termu','ju','eu','rad','stark','waals']
                
    wav1: longueur d'onde en Angstrom
    loggf: log de la force d'oscillateur g_f
    confx, termx, jx, ex: configuration, terme quantique, j et énergie du niveau
      x (l/u, lower/upper)
    rad, stark, waals: damping parameters
    
    Input: file (str) le nom (et path si nécessaire) du fichier de VALD
           wlim (float) la longueur d'onde maximale à traiter
    """

    # dict élément -> numéro atomique
    Elem_dict ={'C':  6,
                'N':  7,
                'O':  8,
                'F':  9,
                'Ne':10,
                'Na':11,
                'Mg':12,
                'Al':13,
                'Si':14,
                'P': 15,
                'S': 16,
                'Cl':17,
                'Ar':18,
                'K': 19,
                'Ca':20,
                'Sc':21,
                'Ti':22,
                'V': 23,
                'Cr':24,
                'Mn':25,
                'Fe':26,
                'Co':27,
                'Ni':28,
                'Cu':29,
                'Zn':30,
                'Ga':31,
                'Ge':32,
                'As':33,
                'Se':34}
                
    L_dict = {'S':0,
              'P':2,
              'D':4,
              'F':6,
              'G':8,
              'H':10,
              'I':12,
              'K':14,
              'L':16,
              'M':18,
              'N':20,
              'O':22}
    
        
    # ouvre le fichier de VALD
    valdlinelist = open(file,'r')
    ev2cm = 1.23984245e-4 # conversion eV -> cm**-1

    # pd.DataFrame() for future reference
    # Data = [[] for _ in range(Elem_dict[Elem])]
    nIonMax = 3
    Data = [[] for _ in range(nIonMax)]

    # read two first lines
    valdlinelist.readline()
        # Lande factors      Damping parameters
    valdlinelist.readline()
        # Elm Ion      WL_air(A)   log gf* E_low(eV) J lo  E_up(eV) J up  lower  upper   mean   Rad.  Stark  Waals

    # on itère sur les lignes. chaque raie est décrite sous 4 lignes,
    # donc les lignes sont traitées différemment
    for i,line in enumerate(valdlinelist.readlines()):
    
        # la première ligne
        if i%4 == 0:
            read_err = False
            # on sépare selon les ','
            Line = line.strip('\n').split(',')
            
            # on enlève les "'" du terme de l'ion, et on sépare élément et ionisation
            Ion = Line[0].replace("'",'').split()
            
            # l'état d'ionisation devient un entier
            try:
                Ion[1] = int(Ion[1])
            except ValueError:
                print(Elem,i)
                print(line)
                raise
            
            # on assigne les valeurs restantes à leur variable
            wav1, loggf, el, jl, eu, ju, lanl, lanu, lanm, rad, stark, waals = \
            [float(x) for x in Line[1:-1]]
            # print Line[1:-1]
            
        # deuxième ligne
        elif i%4 == 1:
            try:
                # configuration et terme quantique LOWER (niveau inférieur)
                confl, terml = line.replace("'",'').replace('\\ ','_').split()[1:]
            except ValueError:
                # print Elem, i
                # print line.replace('\\ ','_')
                # raise
                read_err = True
            
        # troisième ligne
        elif i%4 == 2:
            try:
                # configuration et terme quantique UPPER (niveau supérieur)
                confu, termu = line.replace("'",'').replace('\\ ','_').split()[1:]
            except ValueError:
                # print Elem, i
                # print line.replace('\\ ','_')
                # raise
                read_err = True
            
        # quatrième ligne
        elif i%4 == 3:        
            # cette ligne ne contient que les références (que l'on n'a pas
            # besoin) donc on ne la traite pas.
            # Par contre, on écrit l'information dans le DataFrame
        
            # print [type(c) for c in newline]
            # print Ion[1]
            
            # Si on a pas wlim dans l'input
            if wlim is None:
                boolwlim = True
            else: # sinon on vérifie si la raie est dans le range
                boolwlim = wav1 < wlim
                
            # si pas de limite sur wav1 OU si wav1 plus petit que wlim
            # ET si loggf plus grand que -10 (si <-10 il y a un problème de
            # format avec le code d'atmosphère)
            
            # On détermine les spins et L
            spins = {}
            ls = {}
            for l_u,term in zip(('l','u'),(terml,termu)):
                try:
                    spins[l_u] = (int(term[0]) -1)/2
                except ValueError:
                    spins[l_u] = (int(term[1]) -1)/2
                    
                if '+' in term or '(' in term:
                    ls[l_u] = np.nan
                    
                elif '[' in term:
                    if '/' in term:
                        # si de la forme N[X/2]
                        ls[l_u] = int(term[term.index('[')+1:term.index('/')])/2
                    else: # si de la forme N[X]
                        ls[l_u] = float(term[term.index('[')+1:term.index(']')])
            
                else:
                    try:
                        ls[l_u] = L_dict[term[1]]
                    except KeyError:
                        ls[l_u] = L_dict[term[2]]
                        
            sl, su = spins['l'], spins['u']
            ll, lu = ls['l'], ls['u']
                        
            if boolwlim and Ion[1]<=nIonMax and not read_err:
                # on ajoute l'information dans une liste
                Data[Ion[1]-1].append([wav1, loggf, confl, terml, sl, ll, jl, el,
                         confu, termu, su, lu, ju, eu, lanl, lanu, lanm,
                         rad, stark, waals])
    
    # on itère sur les états d'ionisation
    for i in range(nIonMax):
        # converti la liste en DataFrame
        Dataf = pd.DataFrame(Data[i],columns=['wav1','loggf','confl','terml',
                'sl','ll','jl','el','confu','termu','su','lu','ju','eu',
                'lanl','lanu','lanm','rad','stark','waals'])
        
        # Enregistre le DataFrame sous .pkl pour plus loin
        Dataf.to_pickle('PKL/{0}_{1}.pkl'.format(Elem, i+1))
            
    valdlinelist.close()
    

if __name__ == '__main__':
    vald_convert('./RAW/O.vald', 'O')

