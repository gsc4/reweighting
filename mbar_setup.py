import os
import numpy as np
import model_builder as mdb
import simulation.calc.observables as observables
import pymbar
import mdtraj as md

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

global kb_kJ_mol
kb_KJ_mol = .0083145

def read_filenames(names_file, path):
    
    """ 
    Reads the data found in names_file to get the names we'll need later to read in desired data

    Inputs:          

        names_file - (string) name of file storing the beginning of the names for all the files from which data is to be read

    Outputs:
    
        names - (list) contains the filename beginnings for the data to be read in later

    """

    # Get to the right place
    cwd = os.getcwd()
    os.chdir(path)

    # The file storing names has to be in a very specific format for this to work (newline after every name)
    names  = [ x.strip("\n") for x in open(str(names_file), "r").readlines() ]  
    
    # Take me home
    os.chdir(cwd)    

    return names

def read_data(names, path, stride=1):

    """
    Reads the data found in the files defined in <names-path>

    Inputs:

        names - (list) contains the filename beginnings for data files to be read in

        path - (string) contains the filename end and path for the data files to be read in

        stride (optional) - (int) considers only the data every stride distance

    Outputs:

        data - (numpy.ndarray) - matrix with all the data not yet put into a single vector

        data_long - (numpy.ndarray) - a 1xn (n is number of datapoints read in)  array of all the data read in 


    """ 

    names_strided = names[::stride]

    data = np.array([ np.loadtxt( str(path) + ("{}_0/Etot.xvg").format(x), usecols=[1]) for x in names_strided ] )
    
    data_long = np.concatenate(data)

    return (data, data_long)

def temp_series_u_kn_N_k(temps, energies, Etot, step_size):

    """
    Creates u_kn for data taken at a range of different temperatures        

    Inputs:
        
        temps - (list of strings) holds all the temps data was taken from

        energies - (numpy.ndarray) - holds 

    """
    
    T = np.arange(int(temps[0]), int(temps[-1]) + step_size, step_size)
    
    # Construct N_k
    N_k = np.zeros(len(T))
    count = 0
    for i in range(len(T)):
        if np.isclose(T[i], np.array(temps, float)).any():
            N_k[i] = len(energies[count])
            count += 1
        else:
            None #since N_k is initialzed to 0
    
    # Construct u_kn
    u_kn = np.zeros((len(T), len(Etot)), float)
    for i in range(len(T)):
        B = 1/(kb_KJ_mol*T[i])
        u_kn[i] = Etot*B
    
    return u_kn, N_k, T

def calc_Cv(mbar, Etot, T):
    
    """
    Calculates Cv at a temperature close to the the temps where the data for Etot was taken

    Inputs:
    
        mbar - (mbar object) need this for mbar stuff to work

        Etot - (numpy.ndarray) Energies from all trials concatenated together

    Outputs:

        Cv - (numpy.ndarray) array containing Cv as a function of T

    """

    B = 1/(kb_KJ_mol*T)    

    # Calculate expectations for E
    (Eavg, dEavg) = mbar.computeExpectations(Etot)

    # Calculate expectations for E^2
    (E2avg, dE2avg) = mbar.computeExpectations(Etot**2)

    # Cv = (1/kb*T^2)*( <E^2> - <E>^2)
    Cv = kb_KJ_mol*(B**2)*(E2avg - Eavg**2)

    return Cv

def make_umbrella_u_kn(r0, path, step_size, energies, Etot, atom_pair): 
   
    # I'm not sure if we can actually interpolate points this way or not 
    r0_interped = np.arange(float(r0[0]), float(r0[-1]) + step_size, step_size)
 
    # Construct N_k
    N_k = np.zeros(len(r0))
    count = 0
    for i in range(len(r0)):
        if np.isclose(r0_interped[i], np.array(r0, float)).any():
            N_k[i] = len(r0[count])
            count += 1
        else:
            None #since N_k is initialzed to 0
    
    r0 = read_filenames(r0, path)

    # This is specific to SH3 - 1st bead index is 0, last bead index is 464
    # atom_pair = np.array([[0, 464]])

    energies, Etot = read_data(centers, '~~~~PATH~~~~')

    os.chdir('~~~~Path to directory to contiain r1N')

    # Initialize array for u_kn and counter
    u_kn = np.zeros((len(r0_interped),len(Etot)), float)
    count = 0

    for k in r0:

        os.chdir(k)

        # Read end-end distances from file, if file doesn't exist, make one

        if not os.path.exists('r1N.npy'):

            traj = md.load('traj.xtc', top='conf.gro')
            r1N = md.compute_distances(traj, np.array)
            np.save('r1N.npy', r1N)

        else:
            r1N = np.load('r1N.npy')
        
        # Calculate the potential bias from umbrella and remove it
        u_bias = .5*kb_KJ_mol*(r1N - k)**2

        energy_unbias = energies - u_bias

        Etot_unbias = np.concatenate(energy_unbias)

        u_kn[count] = Etot_unbias_temp

        count += 1

        # Hop back out to parent directory so we can do it again
        os.chdir('..')
 
    return N_k, u_kn, r0_interped
    
  
### NO WAY THIS WORKS ###    

def make_me_some_spaghetti(protein_name, path_to_contacts, Q_long, numbins, thermo_states,  mbar ):
    
    cwd = os.getcwd

    # Take me to the contacts file, my body is ready
    os.chdir(path_to_contacts)

    contacts = np.loadtxt(protein_name + ".contacts")

    model, fitopts = mdb.inputs.load_model(protein_name)

    ri_0 = model.pairwise_distances
    
    # Make histogram    
    (bin_n, bin_edges) = np.histogram(Q , bins = numbins)
    bin_p = bin_n/len(Q_long)
    min_mid = .5*(bin_edges[1:] + bin_edges[:-1])

    # Get expectation of h function
    h_expectations = np.zeros((len(thermo_states), numbins), float) 
    for i in range(numbins):
        h = (Q_long > bin_edges[i]) & (Q_long < bin_edges[i+1])
        h_expectation[:,i] = np.sum(h)*1./len(h)
    
    gamma = 5
    width = 2./gamma

    for i in range(len(contacts)):
        qi = observables(top, contacts[i], ri_0, width)
        qtanh = observables.calculate_observable(trajfiles, qtanhsum_obs, saveas="Q.npy")
        for j in range(len(numbins)):
            Qij, dQij = mbar.computeExpectations(qi*h_expectation[j])/h_expectation[j]
            # Save Qij for all thermo states somehow, need to ask klubes about all of this first
    
# Test functions by playing in the sandbox

if __name__ == "__main__":

    """

    centers = read_filenames(umbrella_last, '~/scratch/SH3.constant.temp/test_melt/6-12-15_umbrella/kumb_0.05/')

    Q, Q_long = read_data(centers,  '~/scratch/SH3.constant.temp/test_melt/6-12-15_umbrella/kumb_0.05/Q.dat'
    
    """

    temps = read_filenames('short_temps', '/home/gsc4/scratch/SH3.constant.temp/test_melt')

    energies, Etot = read_data(temps, '/home/gsc4/scratch/SH3.constant.temp/test_melt/')

    u_kn, N_k, T = temp_series_u_kn_N_k(temps, energies, Etot, 1)   

    mbar = pymbar.MBAR(u_kn, N_k)

    Cv = calc_Cv(mbar, Etot, T)    

    # Generate plot of Cv vs T
    plt.figure()
    plt.plot(T,Cv)
    plt.title('Cv vs T')
    plt.xlabel('T (K)')
    plt.ylabel('Cv')

    plt.savefig("Cv_vs_T.pdf") 
