import os
import numpy as np
import mdtraj as md
import mbar_setup as mbsu

if __name__ == "__main__":
    
    r0 = mbsu.read_filenames('short_temps', '~~~~PATH~~~~')

    #This is specific to SH3 - 1st bead index is 0, last bead index is 57
    atom_pair = np.array([[0, 57]])

    energies, Etot = read_data(centers, '~~~~PATH~~~~')

    os.chdir('~~~~Path to get me to mostly the right place')

    for i in r0:

        os.chdir(r0[i])

        if not os.path.exists('r1N.npy'):
        
            traj = md.load('traj.xtc', top='conf.gro')
            r1N = md.compute_distances(traj, np.array
            np.save('r1N.npy', r1N)
        else:
            r1N = np.load('r1N.npy')

        
        
        

        os.chdir('..')
    
