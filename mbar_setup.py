def read_filenames(names_file, path):
    
    import os
    
    """ 
    Reads the data found in names_file to get the names we'll need later to read in desired data

    Inputs:          

        names_file - (string) name of file storing the beginning of the names for all the files from which data is to be read

    Outputs:
    
        names - (list) contains the filename beginnings for the data to be read in later

    """

    # Get to the right place
    os.chdir(path)

    # The file storing names has to be in a very specific format for this to work (newline after every name)
    names  = [ x.strip("\n") for x in open(str(names_file), "r").readlines() ]  

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

    data = np.array([np.loadtxt(("{}" + str(path)).format(x.rstrip("\n")[::stride] for x in names ))])

    data_long = np.concatenate(data)

    return (data, data_long)

def temp_series_u_kn_N_k(temps, energies, Etot, step_size):

    """
    Creates u_kn for data taken at a range of different temperatures        

    Inputs:
        
        temps - (list of strings) holds all the temps data was taken from

        energies - (numpy.ndarray) - holds 

    """
    
    kb = .0083145 #kJ/mol
    
    T = np.arange(int(temps[0], int(temps[len(temps)-1])+step_size), step_size)
    
    # Construct N_k
    N_k = np.zeros(len(T))
    iter = 0
    for i in range(len(T)):
        if np.isclose(T[i], np.array(temps, float)).any():
            N_k[i] = len(energies[iter])
            iter += 1
        else:
            None #since N_k is initialzed to 0
    
    # Construct u_kn
    for i in range(len(T)):
        B = 1/(kb*T[i])
        u_kn[i] = Etot*B
    
    return (N_k, u_kn)

def calc_Cv(mbar, Etot, T):
    
    """
    Calculates Cv at a temperature close to the the temps where the data for Etot was taken

    Inputs:
    
        mbar - (mbar object) need this for mbar stuff to work

        Etot - (numpy.ndarray) Energies from all trials concatenated together

    Outputs:

        Cv - (numpy.ndarray) array containing Cv as a function of T

    """

    kb = .0083145
    B = 1/(kb*T)    

    # Calculate expectations for E
    (Eavg, dEavg) = mbar.computeExpectations(Etot)

    # Calculate expectations for E^2
    (E2avg, dE2avg) = mbar.computeExpectations(Etot**2)

    # Cv = (1/kb*T^2)*( <E^2> - <E>^2)
    Cv = kb*(B**2)*(E2avg - Eavg**2)

    return Cv

#def make_umbrella_u_kn():
