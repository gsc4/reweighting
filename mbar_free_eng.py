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
    temps = [ x.rstrip("\n") for x in open("short_temps_about_melt", "r").rea

    # Read Q files
    Q = [np.loadtxt("{}_0/Q.dat".format(x.rstrip("\n"))) for x in temps ]
    
    #this is acceptable because all Q calculations have same length
    num_frames = len(Q[0])
   
    # Make my histogram thingy 
    numbins = 100

    bin_n = np.zeros((len(Q), numbins),float)
    bin_p = np.zeros((len(Q), numbins),float)
    bin_edges = np.zeros((len(Q), numbins+1))
    bin_mid = np.zeros((len(Q), numbins))

    for i in range(len(Q)):
        (bin_n[i] , bin_edges[i]) = np.histogram(Q[i], bins = numbins)
        bin_p[i] = bin_n[i]/len(Q[0])
        bin_mid[i] = .5*(bin_edges[i][1:] + bin_edges[i][:-1])



     # Turn list of temperatures into array of floats with a set step size
    step_size = .01
    T = np.arange(int(temps[0]),int(temps[len(temps)-1])+step_size,step_size)

    # Get energies from Etot.xvg files use these to create vector N_k
    energies = [ np.loadtxt("{}_0/Etot.xvg".format(x), usecols=[1]) for x in
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

    # for loop to fill in u_kn
    for i in range(len(T)):
        B = 1/(kb*T[i])
        u_kn[i] = Etot*B
    
    # Initialize mbar
    mbar = pymbar.MBAR(u_kn, N_k)




