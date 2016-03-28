import numpy as np
import pymbar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

if __name__ == "__main__":
    
    #### READ FILES, SET CONSTANTS, INITIALIZE MBAR ####

    # Set boltzmann constant
    kb = .0083145 #kJ/mol*K

    # Temps being considered
    temps = [ x.rstrip("\n") for x in open("short_temps_about_melt", "r").readlines() ]
    
    stride = 100

    # Read Q files
    Q = np.array([np.loadtxt("{}_0/Q.dat".format(x.rstrip("\n")))[::stride] for x in temps ])
    Q_long = np.concatenate(Q)
   
    # Make my histogram thingy 
    numbins = 100

    (bin_n , bin_edges) = np.histogram(Q_long, bins = numbins)
    bin_p = bin_n/len(Q_long)
    bin_mid= .5*(bin_edges[1:] + bin_edges[:-1])
 
     # Turn list of temperatures into array of floats with a set step size
    step_size = .01
    T = np.arange(float(temps[0]),float(temps[len(temps)-1])+step_size,step_size)

    # Get energies from Etot.xvg files use these to create vector N_k
    energies = [ np.loadtxt("{}_0/Etot.xvg".format(x), usecols=[1])[::stride] for x in temps ]
    N_k = np.zeros(len(T))
    iter = 0
    for i in range(len(T)):
        if i % (.1/step_size) == 0:
            N_k[i] = len(energies[iter])
            iter += 1
        else:
            None #since it's already initialized to 0
    
    Etot = np.concatenate(energies)

    # Beta for calculating reduced potentials
    B = 1/(kb*T)

    u_kn = np.zeros((len(T),len(Etot)))

    # for loop to fill in u_kn
    for i in range(len(T)):
        B = 1/(kb*T[i])
        u_kn[i] = Etot*B
 
    # Initialize mbar
    mbar = pymbar.MBAR(u_kn, N_k)

    Q_avg = np.zeros((len(T),numbins),float)
    dQ_avg = np.zeros((len(T),numbins),float)

    #for i in range(len(bin_edges)-1):
    for i in [10]:
        h = (Q_long > bin_edges[i]) & (Q_long < bin_edges[i+1])
        Q_avg[:,i], dQ_avg[:,i] = mbar.computeExpectations(h.astype(int))
     



