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
    Reads names_file to make a list of names

    Keyword arguments:
    names_file -- (string) name of file containing beginnings of other filenames

    Returns:
    names -- (list) contains the filename beginnings for the data to be read in later
    """

    # Get to the right place
    cwd = os.getcwd()
    os.chdir(path)

    # The file storing names has to be in a very specific format for this to work (newline after every name)
    names  = [ x.strip("\n") for x in open(str(names_file), "r").readlines() ]  
    
    # Take me home
    os.chdir(cwd)    

    return names

def read_data(name_start, name_end, path, stride=1):
    """
    Reads the data found in the named files on the given path

    Keyword arguments:
    names_start -- (list) start of names of files to read
    name_end -- (string) end of names of files to read (file type included)
    path -- (string) path to files to be read
    stride (optional) -- (int) only take  data every stride distance

    Returnss:
    data -- (numpy.ndarray) array of read in data
    data_long -- (numpy.ndarray) 1xn array of all the data read in (a vector)
    """ 

    names_strided = name_start[::stride]

    data = np.array([ np.loadtxt( str(path) + ("{}{}").format(x, name_end), usecols=[1]) for x in names_strided ] )
    
    data_long = np.concatenate(data)

    return (data, data_long)

def temp_series_u_kn_N_k(temps, energies, Etot, num_interp = 2):
    """
    Creates u_kn and N_k, interpolates 2 points by default        

    Keyword arguments:
    temps -- (list of strings) temperatures at which data was taken
    energies -- (numpy.ndarray)  kxn array of energy data
    Etot -- (numpy.ndarray) 1x(k*n) array of energy data (vector)
    num_interp (optional) -- (int) number interpolating points

    Returns:
    u_kn -- (numpy.ndarray) every frame evaluated at every thermodynamic state
    N_k -- (numpy.ndarray) number of actual data frames for corresponding T
    T -- (numpy.ndarray) all temperatures (including interpolated ones)

    """
    
    flt_temps = [ float(x) for x in temps ]
 
    for i in range(len(flt_temps) - 1):
        temp_diff = flt_temps[i+1] - flt_temps[i]
        for j in range(num_interp): 
            flt_temps.append(flt_temps[i] + (j+1)*temp_diff/(num_interp + 1))

    T = np.array(flt_temps)

    # Construct N_k
    N_k = np.zeros(len(T))
    for i in range(len(temps)):
            N_k[i] = len(energies[i])

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

def umbrella_u_kn_N_k(r0, energies, T, atom_pair, path, num_interp=2):
    
    cwd = os.getcwd()
      
    E_unbias = []
    r1N_all = []

    flt_r0 = [ float(x) for x in r0 ]
 
    for i in range(len(r0) - 1):
        r0_diff = flt_r0[i+1] - flt_r0[i]
        for j in range(num_interp): 
            flt_r0.append(flt_r0[i] + (j+1)*r0_diff/(num_interp + 1))

    r0_interp = np.array(flt_r0)
 
    # Construct N_k
    N_k = np.zeros(len(r0_interp) + 1)
    for i in range(len(r0)):
        N_k[i] = len(energies[i])
     
    os.chdir(path)    

    for i in range(len(r0)):
        
        # Hop into directory w/ data    
        os.chdir(r0[i])
        
        ###### NEED TO READ THIS VALUE FROM SOMEWHERE ###### 
        k_umb = 2   
        
        # Read end-end distances from file, if file doesn't exist, make one
        if not os.path.exists('r1N.npy'):
            traj = md.load('traj.xtc', top='conf.gro')
            r1N = md.compute_distances(traj, atom_pair)
            np.save('r1N.npy', r1N)
        else:
            r1N = np.load('r1N.npy')        
        
        r1N_all.extend([x for y in r1N.tolist() for x in y])

        Ebias = .5*k_umb*(( r1N - float(r0[i]) )**2)
        E_unbias_temp_pre = energies[i] - np.transpose(Ebias)
        E_unbias_temp=[x for y in E_unbias_temp_pre.tolist() for x in y]
        E_unbias.extend( E_unbias_temp )
        os.chdir('..')

    r1N_all = np.array(r1N_all)
    E_unbias = np.array(E_unbias)

    beta = 1./(kb_KJ_mol*T)

    # Initialize array for u_kn and value for indexing
    u_kn = np.zeros((len(r0_interp) + 1 ,len(E_unbias)), float)

    for k in range(len(r0_interp)):
 
        ###### NEED TO READ THIS VALUE FROM SOMEWHERE ###### 
        k_umb = 2   
        
        E_bias_k = .5*k_umb*((r1N_all - r0_interp[k])**2)
        
        u_kn[k,:] = beta*(E_unbias + E_bias_k)

    u_kn[-1,:] = beta*E_unbias

    # Hop back out to parent directory so we can do it again
    os.chdir(cwd)
        
    return u_kn, N_k, r0_interp 

def free_eng_profile(mbar, Q_long,  numbins=100):

    bin_n, bin_edges = np.histogram(Q_long, bins=numbins)
    bin_p = bin_n/len(Q_long)
    bin_mid = .5*(bin_edges[1:] + bin_edges[:-1])

    p_of_Q = np.zeros( (len(r0_interp), numbins), float)
    p_of_Q_min_0 =  np.zeros( (len(r0_interp), numbins), float)
    dp_of_Q = np.zeros( (len(r0_interp), numbins), float)

    for i in range(numbins):
        h = (Q_long > bin_edges[i]) & (Q_long < bin_edges[i+1])
        p_of_Q[:,i], dp_of_Q[:,i] = mbar.computeExpectations(h.astype(int))
        p_of_Q_min_0[:,i] = p_of_Q[:,i] - min(p_of_Q[:,i])

    A_over_RT = -np.log(p_of_Q)
    A_over_RT_min_0 = np.zeros(A_over_RT.shape)
    for i in range(numbins):
        A_over_RT_min_0[:,i] = A_over_RT[:,i] - min(A_over_RT_min_0[:,i])

    return A_over_RT, A_over_RT_min_0



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
   
###############################################################################
#                                    TEST                                     #
###############################################################################
 
# Test functions by playing in the sandbox

if __name__ == "__main__":

    ##### Test umbrella stuff #####

    path = '/home/gsc4/scratch/SH3.constant.temp/test_melt/umbrella_sampling/4-7-2016_end_end_umb/'

    T_umb = 129.0

    atom_pair = np.array([0,57])

    r0 = read_filenames('umbrella_last', path)
    
    energies, poop = read_data(r0, '/Etot.dat', path)
    
    u_kn, N_k, r0_interp = umbrella_u_kn_N_k(r0, energies, T_umb, atom_pair, path)
     
    mbar = pymbar.MBAR(u_kn, N_k)
    
    #Q, Q_long = read_data(r0, '/Q.dat', path)

    
 
    ##### Test temperature stuff #####
    """
    temps = read_filenames('short_temps', '/home/gsc4/scratch/SH3.constant.temp/test_melt')

    energies, Etot = read_data(temps, '_0/Etot.xvg', '/home/gsc4/scratch/SH3.constant.temp/test_melt/')

    u_kn, N_k, T = temp_series_u_kn_N_k(temps, energies, Etot, 1)   

    mbar = pymbar.MBAR(u_kn, N_k)

    Cv = calc_Cv(mbar, Etot, T)    
    """

    # Plot stuff
    """
######## This seems like a super circuitious route to do this... ####
    T_Cv_zip = zip(T,Cv)
    T_Cv_sort = sorted(T_Cv_zip, key=lambda x:x[0])
    T = [ x[0] for x in T_Cv_sort ] 
    Cv = [ x[1] for x in T_Cv_sort ]
#####################################################################
 
    plt.figure()
    plt.plot(T,Cv)
    plt.title('Cv vs T')
    plt.xlabel('T (K)')
    plt.ylabel('Cv')

    plt.savefig("Cv_vs_T.pdf") 
    """
