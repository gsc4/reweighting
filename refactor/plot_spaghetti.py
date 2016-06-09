import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import fnmatch

from Tkinter import Tk
from tkFileDialog import askopenfilename
from tkFileDialog import askdirectory

if __name__ == "__main__":

    path_start = '/home/gsc4/scratch/SH3/'
    os.chdir(path_start)

    os.system('clear')
    print('Go into the directory containing our data')
    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    path_next =  askdirectory()
    os.chdir(path_next)

    os.system('clear')
    print 'Select loops file'
    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    loops_name = askopenfilename() # show an "Open" dialog box and return the path to the selected file
    loops = np.load(loops_name)    
    max_loop = np.max(loops)

    os.system('clear')
    print 'Select bin mid file'
    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    bin_mid_name = askopenfilename() # show an "Open" dialog box and return the path to the selected file
    bin_mid = np.load(bin_mid_name)
       
    os.system('clear')
    print 'Select spaghetti data file'
    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    Qi_vs_rxn_coord_name = askopenfilename() # show an "Open" dialog box and return the path to the selected file
    Qi_vs_rxn_coord = np.load(Qi_vs_rxn_coord_name)

    for i in range(len(Qi_vs_rxn_coord)): 
        plt.plot(bin_mid, Qi_vs_rxn_coord[i,:], c=cm.jet(loops[i]/max_loop))

    plt.ylabel('Percent Formation of Contact')

    print
    rxn_coord_name = raw_input('Reaction coordinate name (with units): ')
    plt.xlabel(rxn_coord_name)
    
    print
    title_name = raw_input('Graph title : ')
    plt.title(title_name)

    print
    save_name = raw_input('Save filename: ')
    
    if fnmatch.fnmatch(save_name, '*.pdf'):
        plt.savefig(save_name)
    else:
        plt.savefig(save_name + '.pdf')

