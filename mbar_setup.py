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
    Reads file to make a list of names

    Parameters
    ----------
    names_file : string
        Name of file
    path : string
        Path from to file location

    Returns
    -------
    names : list of strings
        Identifying name parts of files to be read in later
    """
    
    cwd = os.getcwd()
    os.chdir(path)
    
    names = [ x.strip("\n") for x in open(str(names_file), "r").readlines() ]  
    
    os.chdir(cwd)    

    return names

def read_data(name_start, name_end, path, stride=1):
    """
    Reads the data found in the named files on the given path

    Parameters
    ----------
    names_start : list of strings
        Identifying name parts (the beginning) of files to read
    name_end : string
        End of files to read including filetype (common to all files)
    path : string
        Path to file location
    stride : int, optional
        Only consider data every stride distance

    Returns
    -------
    data : numpy.ndarray
        Array of read-in data
    data_long : numpy.ndarray
        1xn array of all the read-in data (a vector)
    """ 

    names_strided = name_start[::stride]
    
    if name_end == '/Q.dat':
        data = np.array([np.loadtxt( str(path) + ("{}{}").format(x, name_end))\
             for x in names_strided ])
    else:
        data = np.array([np.loadtxt( str(path) + ("{}{}").format(x, name_end),
             usecols=[1]) for x in names_strided ] )
    
    data_long = np.concatenate(data)

    return (data, data_long)

def temp_series_u_kn_N_k(temps, energies, Etot, num_interp = 2):
    """
    Creates inputs for mbar on constant temperature simulations

    Parameters
    ----------
    temps : list of strings
        Temperatures at which data was taken
    energies : numpy.ndarray
        Array of energy data
        Dimensions: (num thermo states) x (num frames)     
    Etot : numpy.ndarray
        Vector of energy data
        Dimensions: 1 x ((num thermo states)*(num frames)) 
    num_interp : int, optional
        Number interpolating points

    Returns
    -------
    u_kn : numpy.ndarray
        Every frame evaluated at every thermodynamic state
        Dimensions: (num thermo states) x ((num thermo states)*(num frames))
    N_k : numpy.ndarray
        Number of data frames for each temperature in `T`
    T : numpy.ndarray
        All temperatures (including interpolated ones)
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
    Calculates heat capacity for given temperatures

    Parameters
    ----------
    mbar : mbar object
        Needed for calcualtions using pymbar      
    Etot : numpy.ndarray
        Vector of energy data
        Dimensions: 1 x ((num thermo states)*(num frames)) 
    T : numpy.ndarray
        All temperatures (including interpolated ones)
    
    Returns
    -------
    Cv : numpy.ndarray
        Array of heat capacity as a function of T
        Dimensions: 1 x (num temps in `T`)
    """

    beta = 1/(kb_KJ_mol*T)    

    # Calculate expectations for E
    (Eavg, dEavg) = mbar.computeExpectations(Etot)

    # Calculate expectations for E^2
    (E2avg, dE2avg) = mbar.computeExpectations(Etot**2)

    # Cv = (1/kb*T^2)*( <E^2> - <E>^2)
    Cv = kb_KJ_mol*(beta**2)*(E2avg - Eavg**2)

    return Cv

def umbrella_u_kn_N_k(r0, energies, T, atom_pair, path, num_interp=2):
    """
    Creates inputs for mbar on umbrella sampled simulations

    Parameters
    ----------
    r0 : list of strings
        Umbrella centers for which data was taken            
    energies : numpy.ndarray
        Array of energy data
        Dimensions: (num thermo states) x (num frames)
    atom_pair : np.ndarray   
        Numbers of the 2 atoms that were pulled by umbrella potential
        Dimensions: 1 x 2
    num_interp : int, optional
        Number of interpolating points

    Returns
    -------
    u_kn : numpy.ndarray
        Every frame evaluated at every thermodynamic state
        Dimensions: (num thermo states) x ((num thermo states)*(num frames))
    N_k : numpy.ndarray
        Number of data frames for each temperature in `r0_interp`
    r0_interp : numpy.ndarray
        All umbrella centers (including interpolated ones)
    r1N_all : numpy.ndarray
        Vector of all pairwise distances for the 2 atoms for all thermo states
        Dimensions: 1 x ((num thermo states)*(num frames))
 
    """

    # Check if num_interp is a number
    if not isinstance(num_interp, int): 
        try:
            print "Rouding interpolations to nearest non-negative integer"
            num_interp = int(round(float(num_interp)))
            if num_interp < 0:
                num_interp = 0
        except TypeError:
            raise Exception("num_interp has to be a number (preferably int)")

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

    # Initialize array for u_kn
    u_kn = np.zeros((len(r0_interp) + 1 ,len(E_unbias)), float)

    for k in range(len(r0_interp)):
 
        ###### NEED TO READ THIS VALUE FROM SOMEWHERE ###### 
        k_umb = 2   
        
        E_bias_k = .5*k_umb*((r1N_all - r0_interp[k])**2)
        
        u_kn[k,:] = beta*(E_unbias + E_bias_k)

    u_kn[-1,:] = beta*E_unbias

    # Hop back out to parent directory so we can do it again
    os.chdir(cwd)
        
    return u_kn, N_k, r0_interp, r1N_all 

def free_eng_profile(mbar, rxn_coord, numbins=100):
    """
    Calculates free energy as a function of a reaction coordinate

    Parameters
    ----------
    mbar : mbar object
        Needed for calculations using pymbar
    rxn_coord :
        Vector of reaction coordinate data for all frames in all thermo states
    numbins : int, optional
        Number of bins used in the indicator function `h`

    """

    bin_n, bin_edges = np.histogram(rxn_coord, bins=numbins)
    bin_mid = .5*(bin_edges[1:] + bin_edges[:-1])

    p_of_rxn_coord = np.zeros( (len(r0_interp) + 1, numbins), float)
    p_of_rxn_coord_min_0 =  np.zeros( (len(r0_interp) + 1, numbins), float)
    dp_of_rxn_coord = np.zeros( (len(r0_interp) + 1, numbins), float)

    for i in range(numbins):
        h = (rxn_coord > bin_edges[i]) & (rxn_coord < bin_edges[i+1])
        p_of_rxn_coord[:,i], dp_of_rxn_coord[:,i] = mbar.computeExpectations(h.astype(int))
        p_of_rxn_coord_min_0[:,i] = p_of_rxn_coord[:,i] - min(p_of_rxn_coord[:,i])

    A_over_RT = -np.log(p_of_rxn_coord)
    A_over_RT_min_0 = np.zeros(A_over_RT.shape)
    for i in range(A_over_RT.shape[0]):
        A_over_RT_min_0[i,:] = A_over_RT[i,:] - np.nanmin(A_over_RT[i,:])

    return A_over_RT, A_over_RT_min_0, bin_mid

### NO WAY THIS WORKS ###     
def make_me_some_spaghetti(protein_name, path_to_contacts, Q_long, numbins, thermo_states,  mbar ):
    """    

    HOT GARBAGE
    
    """

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
        qi = observables.TanhContacts(top, contacts[i], ri_0[i], width)
        qtanh = observables.calculate_observable(trajfiles, qtanhsum_obs)
        for j in range(numbins):
            h = 'which frames are in bin j' 'can move above loop 246 to here'
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
    
    u_kn, N_k, r0_interp, r1N_all = umbrella_u_kn_N_k(r0, energies, T_umb, atom_pair, path, num_interp=0)
    """ 
    mbar = pymbar.MBAR(u_kn, N_k)
    
    # Turns out I don't need this at all...  
    #Q, Q_long = read_data(r0, '/Q.dat', path)

    A_over_RT, A_over_RT_min_0, bin_mid = free_eng_profile(mbar, r1N_all)

    ##### Plot stuff #####
    
    plt.figure()
    for i in range(len(r0_interp) + 1):
        if i == len(r0_interp):
            plt.plot(bin_mid,A_over_RT_min_0[i,:],linewidth=5)

        else:
            plt.plot(bin_mid, A_over_RT_min_0[i,:])  
 
    centers = [ 'r0 = ' + str(x) + ' nm'  for x in r0_interp ]  
    
    centers.append('unbiased')
 
    plt.title('A/RT vs r1N')
    plt.xlabel('r1N (nm)')
    plt.ylabel('A/RT')
    plt.legend(centers, loc=9, prop = {'size':9})
    
    plt.savefig("A_over_RT_vs_r1N_min_0.pdf") 
    """

    ##### Test temperature stuff #####
    """
    temps = read_filenames('short_temps', '/home/gsc4/scratch/SH3.constant.temp/test_melt')

    energies, Etot = read_data(temps, '_0/Etot.xvg', '/home/gsc4/scratch/SH3.constant.temp/test_melt/')

    u_kn, N_k, T = temp_series_u_kn_N_k(temps, energies, Etot, 1)   

    mbar = pymbar.MBAR(u_kn, N_k)

    Cv = calc_Cv(mbar, Etot, T)    
    """

    ##### Plot stuff #####
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
