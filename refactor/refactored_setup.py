import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import itertools
import time
import fnmatch

import mdtraj as md
import pymbar

import simulation.calc.observables as observables

global kb_kJ_mol
kb_KJ_mol = .0083145

def read_data_const_temp(name_format, name_data, n_trials=1):
    """
    Reads data from files and concatenates them
    Parameters
    ----------
    name_format : string
        Blank model for data filenames
    name_data : list of list of string
        Fileame parts to insert into the blank model (state names, trial #, etc.)
    n_trials : int, optional
        Number of trials run at each thermodynamic state
    Returns
    -------
    data_long : numpy.ndarray
        Data vector
    n_frames : long
        Numer of frames considered at each thermodynamic state
    """

    name_data_comb = [ i for i in itertools.product(*name_data) ]   
    
    if name_data[-1][0] in ('Q.dat', 'qtanh.dat'):
        data = np.array([ np.loadtxt(name_format.format(*i)) for i  in name_data_comb ])
    else:
        data = np.array([ np.loadtxt(name_format.format(*i), usecols=[1]) for i in name_data_comb ])

    data_concat = np.concatenate(data)

    n_frames = data.shape[1]*n_trials

    return data_concat, n_frames

def read_data_umb(dir_name_format, dir_name_data, prot_name, rxn_coord_filename, umb_centers, atom_pair, n_trials=1):
    """
    Reads data from files and removes biasing potential
    Parameters
    ----------
    dir_name_format : string
        Blank model for directory names
    dir_name_data
        Directory name parts to insert into blank model 
    rxn_coord_filename : string
        Name of file containing rxn coordinate data
    umb_centers : list of strings
        Centers for umbrella biasing potentials
    atom_pair : numpy.ndarray
        For atom pulling umbrellas the indices of the 2 residues being pulled
    Returns
    -------
    E_unbias_concat : numpy.ndarray
        Vector of unbiased energy data
    rxn_coord_concat : numpy.ndarray
        Vector of all reaction coord data
    k_umb : numpy.ndarray
        Umbrella constants for each thermodynamic state
    n_frames : long
        Number of frames considered at each thermodynamic state
    """
    
    dir_name_data_comb = [ i for i in itertools.product(*dir_name_data) ]
    dir_names = [ dir_name_format.format(*i) for i in dir_name_data_comb ]

    E_unbias_concat = []
    rxn_coord_concat = []
    k_umb_all = []

    frame_toss = 0

    for i in range(len(dir_names)):
        # Hop into folder with the data
        os.chdir(dir_names[i])
        
        if os.path.exists('umbrella_params'):
            with open('umbrella_params') as fin:
                params = fin.readline().split()
            k_umb = float(params[1])
            frame_toss = float(params[4])/1000
 
        else:
            # Read in k_umb (assuming that there's an empty line at the end of run.mdp) 
            with open('run.mdp','r') as fin:
                k_umb_line = list(fin)[-2]
                for word in k_umb_line.split():
                    try:
                        k_umb = float(word)
                    except ValueError:
                        raise IOError('Couldn\'t find umbmrella params or get k_umb from run.mdp')
    
        k_umb_all.append(k_umb)

        if rxn_coord_filename == 'r1N.npy':
            # Read distances from file, if file doesn't exist, make one
            if not os.path.exists(rxn_coord_filename):
                traj = md.load('traj.xtc', top='conf.gro')
                rxn_coord_data = md.compute_distances(traj, atom_pair)
                np.save(rxn_coord_filename, r1N)
            else: 
                rxn_coord_data = np.load(rxn_coord_filename)
        
        else:
            # Read Q from file, if file doesn't exist, make one
            if not os.path.exists(rxn_coord_filename):
                pairs = np.loadtxt('../' + prot_name + '.contacts', dtype=int) - 1
                n_native_pairs = pairs.shape[0]

                top = 'ca.pdb'
                ref = md.load(top)
                r0 = md.compute_distances(ref, pairs)[0]
                widths = .05* np.ones(n_native_pairs, float)
                
                r0_cot = 1.2*r0
                qtanhsum_obs = observables.TanhContactSum(top, pairs, r0_cont, widths)
                rxn_coord_data = observables.calculate_observables('traj.xtc', qtanhsum_obs)
                np.save(rxn_coord_filename, r1N)
            elif fnmatch.fnmatch(rxn_coord_filename, '*.npy'):
                rxn_coord_data = np.load(rxn_coord_filename)
            else:
                rxn_coord_data = np.loadtxt(rxn_coord_filename)

        # Add all rxn coord data to list
        if isinstance(rxn_coord_data, np.ndarray) and (rxn_coord_data.ndim == 1):
            rxn_coord_concat.extend([ x for x in rxn_coord_data.tolist()])
        else:
            rxn_coord_concat.extend([x for y in rxn_coord_data.tolist() for x in y])

        # Read in energies
        if os.path.exists('Etot.dat'):    
            energies = np.loadtxt('Etot.dat', usecols=[1])[frame_toss:]   
        elif os.path.exists('Epot.dat'):
            energies = np.loadtxt('Epot.dat', usecols=[1])[frame_toss:]
        else:
            raise IOError('Couldn\'t locate an engergy file')

        # Remove biasing potential from total energy
        Ebias = .5*k_umb*( (rxn_coord_data - float(umb_centers[i/n_trials]))**2 )
        if (Ebias.ndim == 1) and (energies.ndim == 1):
            E_unbias = (energies - Ebias).tolist()
        else: 
            E_unbias = [x for y in (energies - Ebias.transpose()).tolist() for x in y]
        E_unbias_concat.extend(E_unbias)
        
        os.chdir('..')

    n_frames = len(rxn_coord_data)*n_trials

    k_umb = np.array(k_umb_all[0::n_trials])

    rxn_coord_concat = np.array(rxn_coord_concat)
    E_unbias_concat = np.array(E_unbias_concat)    

    return E_unbias_concat, rxn_coord_concat, k_umb,  n_frames, frame_toss

def create_mbar_const_temp(temps, E_concat, n_frames, n_interp=2):
    """
    Constructs mbar for constant temp simulations
    
    Parameters
    ----------
    temps : list of strings
        Temperatures at which data was taken
    E_concat : numpy.ndarray
        Vector of energy data
    n_frames : int
        Number of frames considered for each simulation
    n_interp : int, optional
        Number of interpolating points between temperateres in `temps`    
    Returns
    -------
    mbar : mbar object
        Needed for computations with mbar
    T_interp : numpy.ndarray
        All temperatures (including interpolated ones)
    """
    
    flt_temps = [ float(x) for x in temps ]
    
    for i in range(len(flt_temps) - 1):
        temp_diff = flt_temps[i+1] - flt_temps[i]
        for j in range(n_interp):
            flt_temps.append(flt_temps[i] + (j+1)*temp_diff/(n_interp +1))

    T_interp = np.array(flt_temps)

    # Construct N_k and u_kn
    N_k = np.zeros(len(T_interp))
    u_kn = np.zeros((len(T_interp), len(E_concat)), float)
    for i in range(len(T_interp)):

        if i < len(temps):
            N_k[i] = n_frames
        else:
            None

        beta = 1/(kb_KJ_mol*T_interp[i])
        u_kn[i] = E_concat*beta

    mbar = pymbar.MBAR(u_kn, N_k)

    return mbar, T_interp

def create_mbar_umb(umb_centers, E_unbias_concat, rxn_coord_concat, n_frames, k_umb, umb_T, n_interp=2):
    """
    Constructs mbar for umbrella simulations
    Parameters
    ----------
    umb_centers : list of strings
        Umbrella center values for which data was taken
    E_unbias_concat : numpy.ndarray
        Vector of unbiased energy data
    rxn_coord_concat : numpy.ndarray
        Vector of reaction coordinate values 
    n_frames : int
        Number of frames considered for each simulation    
    k_umb : numpy.ndarray
        Umbrella constants at each thermodynamic state
    umb_T : float
        Temperature at which umbrella simulation was run
    n_interp : int, optional
        Number of interpolating points between the ceneters in `umb_centers`
 
    Returns
    -------
    mbar : mbar object
        Needed for computations using mbar    
    umb_centers_interp : numpy.ndarray    
        All umbrella centers (including interpolated ones)    
    """

    flt_umb_centers = [ float(x) for x in umb_centers ] 
    
    for i in range(len(flt_umb_centers) - 1):
        cent_diff = flt_umb_centers[i+1] - flt_umb_centers[i]
        for j in range(n_interp):
            flt_umb_centers.append(flt_umb_centers[i] + (j+1)*cent_diff/(n_interp + 1))

    umb_centers_interp = np.array(flt_umb_centers)

    # Construct N_k and u_kn
    N_k = np.zeros(len(umb_centers_interp)+1)
    u_kn = np.zeros((len(umb_centers_interp) + 1, len(E_unbias_concat)), float)  
    
    beta = 1./(kb_KJ_mol*umb_T)

    for i in range(len(umb_centers_interp)):
        
        if i < len(umb_centers):
            N_k[i] = n_frames
            E_bias_i = .5*k_umb[i]*( (rxn_coord_concat - umb_centers_interp[i])**2 )
        else:
            # For interp points use same k_umb as closest umbrella center with known k_umb
            k_umb_index = min(range(len(umb_centers)), key=lambda x: abs(float(umb_centers[x]) - umb_centers_interp[i]))
            E_bias_i = .5*k_umb[k_umb_index]*( (rxn_coord_concat - umb_centers_interp[i])**2 )
    
        u_kn[i, :] = beta*(E_unbias_concat + E_bias_i)
   
    u_kn[-1, :] = beta*E_unbias_concat
 
    mbar = pymbar.MBAR(u_kn, N_k)

    return mbar, umb_centers_interp

def calc_Cv(mbar, E_concat, T_interp):
    """
    Calculates heat capacity
    Parameters
    ----------
    mbar : mbar object
        Needed for computations using mbar
    E_concat : numpy.ndarray
        Vector of energy data
    T_interp : numpy.ndarray
        All temperatures (including interpolated ones)
 
    Returns
    -------
    Cv : numpy.ndarray
        Vector of heat capacities corresponding to the temperatures in `T`
    """
    
    beta = 1/(kb_KJ_mol*T_interp)

    Eavg, dEavg = mbar.computeExpectations(E_concat)

    E2avg, dE2avg = mbar.computeExpectations(E_concat**2)
    
    Cv = kb_KJ_mol*(beta**2)*(E2avg - Eavg**2)

    return Cv

def empirical_spaghetti(dir_name_format, dir_name_data, prot_name, rxn_coord_concat, numbins=100 ):

    """
    Lukewarm garbage
    """
    
    name_data_comb = [ i for i in itertools.product(*dir_name_data) ]
    dirnames = [ dir_name_format.format(*i) for i in name_data_comb ] 
    
    pairs = np.loadtxg(prot_name + '.contacts', dtype=int) - 1
 
    trajfiles = []
    for i in dirnames:
        trajfiles.append(i + '/traj.xtc')   

    # Calculate contact distances
    if os.path.exists(dirnames[0] + '/ref.pdb'):
        top = dirnames[0] + '/ref.pdb'
    elif os.path.exists(dirnames[0] + '/ca.pdb'):
        top = dirnames[0] + '/ca.pdb'
    else:
        raise IOError('Couldn\'t locate .pdb file') 
    ref = md.load(top)
    r0 = md.compute_distances(ref, pairs)[0]
    r0_cont = 1.2*r0

    bin_n, bin_edges = np.histogram(rxn_coord_concat, bins=numbins)
    bin_mid = .5*(bin_edges[1:] + bin_edges[:-1])

    gamma = 5
    width = 2./gamma

    Qi_vs_rxn_coord = np.zeros((len(pairs), numbins), float) 
    loops = np.zeros(len(pairs), float)

    for i in range(len(pairs)):
        qi = observables.TanhContacts(ref, np.array([pairs[i]]), r0_cont[i], width)
        qi_tanh = np.concatenate(observables.calculate_observable(trajfiles, qi))
        loops[i] = pairs[i][1] - pairs[i][0]
        for j in range(numbins):
            h = (rxn_coord_concat > bin_edges[j]) & (rxn_coord_concat <= bin_edges[j+1])
            if np.any(h):
                Qi_vs_rxn_coord[i,j] = np.mean(qi_tanh[h])
        
    return Qi_vs_rxn_coord, bin_mid, loops

def mbar_spaghetti(dir_name_format, dir_name_data, prot_name, rxn_coord_concat, n_thermo_states, mbar, numbins=100, frame_toss=0):
    
    name_data_comb = [ i for i in itertools.product(*dir_name_data) ]
    dirnames = [ dir_name_format.format(*i) for i in name_data_comb ] 
    
    pairs = np.loadtxt(prot_name + '.contacts', dtype=int) - 1
 
    trajfiles = []
    for i in dirnames:
        trajfiles.append(i + '/traj.xtc')  
 
    # Calculate contact distances
    if os.path.exists(dirnames[0] + '/ref.pdb'):
        top = dirnames[0] + '/ref.pdb'
    elif os.path.exists(dirnames[0] + '/ca.pdb'):
        top = dirnames[0] + '/ca.pdb'
    else:
        raise IOError('Couldn\'t locate pdb file')
    
    ref = md.load(top)
    r0 = md.compute_distances(ref, pairs)[0]
    r0_cont = 1.2*r0

    bin_n, bin_edges = np.histogram(rxn_coord_concat, bins=numbins)
    bin_mid = .5*(bin_edges[1:] + bin_edges[:-1])

    gamma = 5
    width = 2./gamma

    Qi_vs_rxn_coord = np.zeros((n_thermo_states, len(pairs), numbins), float)     
    loops = np.zeros(len(pairs), float)
        
    for i in range(len(pairs)):
    #for i in [0]:
        qi = observables.TanhContacts(ref, np.array([pairs[i]]), r0_cont[i], width)
        qi_tanh = np.concatenate(np.array(observables.calculate_observable(trajfiles, qi))[:,frame_toss:,:])
        loops[i] = pairs[i][1] - pairs[i][0]
        time_start = time.time()
    
        for j in range(numbins):       
            h = (rxn_coord_concat > bin_edges[j]) & (rxn_coord_concat <= bin_edges[j+1])
            if np.any(h):
                h = h.astype(float)
                h_avg, dh_avg = mbar.computeExpectations(h, compute_uncertainty=False)
                numerator_avg, dnumerator_avg = mbar.computeExpectations((qi_tanh.transpose())[0]*h, compute_uncertainty=False) 
                Qi_vs_rxn_coord[:,i,j] = numerator_avg/h_avg
            else:
                Qi_vs_rxn_coord[:,i,j] = np.nan

            print j, time.time() - time_start
            time_start = time.time()

    return Qi_vs_rxn_coord, bin_mid, loops
