import numpy as np
import pymbar
import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt

if __name__ == "__main__":

    #### READ FILES, SET CONSTANTS, AND GET THINGS READY ####    

    # My boy boltzmann has a very special number
    kb = 0.0083145  # kJ/mol K

    # Collect list of temperatures from file
    temps = [ x.rstrip("\n") for x in open("short_temps", "r").readlines() ]
 
    # Turn list of temperatures into array of floats with a set step size
    step_size = .1
    T = np.arange(int(temps[0]),int(temps[len(temps)-1])+step_size,step_size)

    #print len(T)

    # Get energies from Etot.xvg files use these to create vector N_k
    energies = [ np.loadtxt("{}_0/Etot.xvg".format(x), usecols=[1]) for x in temps ]
    N_k = np.zeros(len(T))
    iter = 0
    for i in range(len(T)):
        if i % (1/step_size) == 0:
            N_k[i] = len(energies[iter])
            iter += 1        
        else:
            None #since it's already initialized to 0
    
    #print N_k
    """
    # Put energies in array named Etot 
    Etot = np.array(energies, float)
    """
    
    Etot = np.concatenate(energies) 
    print len(Etot)
    
    # Read Q values for calculations later
    Q = [np.loadtxt("{}_0/Q.dat".format(x.rstrip("\n"))) for x in temps] 
    """
    Q_kln = np.empty((len(T), Etot.shape[0], Etot.shape[1]), float)
    Q_kln[:] = np.array(Q, float)    
    """
    Q = np.concatenate(Q)
    #Q_kn = np.empty((len(T),len(Etot)), float)
    #Q_kn[:] = Q 

    # Beta for calculating reduced potentials
    B = 1/(kb*T) 
    #print len(B)
    
    """
    # Initialize vector to hold reduced potentials
    u_kln = np.zeros((len(T), Etot.shape[0], Etot.shape[1]), float)    
    
    # Loop to fill in u_kln
    for i in range(len(T)):
        u_kln[i] = Etot*B[i]    
    
    
    # It's gonna be unhappy because this does not have K == L
    print u_kln.shape
    """

    u_kn = np.zeros((len(T), len(Etot)), float)

    #for loop to fill in u_kn
    for i in range(len(T)):
        B = 1/(kb*T[i])
        u_kn[i] = Etot*B
    #print u_kn.shape    
    
    # Initialize mbar
    #mbar = pymbar.mbar.MBAR(u_kln, N_k)
    mbar = pymbar.MBAR(u_kn, N_k)


    #### COMPUTE OBSERVABLES #### 
    # For some reason we need to use unreduced potentials, aka Etot over and over
    """
    E_kln = np.empty((len(T), Etot.shape[0], Etot.shape[1]), float)
    E_kln[:] = Etot
    """

    #E_kn = np.empty((len(T), len(Etot)),float)
    #E_kn[:] = Etot
    
    # Calculate expectations for E
    (Eavg, dEavg) = mbar.computeExpectations(Etot)


    # Calculate expectations for E^2
    (E2avg, dE2avg) = mbar.computeExpectations(Etot**2)


    # Cv = (1/kb*T^2)*( <E^2> - <E>^2)
    Cv = kb*(B**2)*(E2avg - Eavg**2)   


    # <Q>
    (Qavg, dQavg) = mbar.computeExpectations(Q)


    # <A> = A(T)
    # didn't do this one
   
    #### MAKE GRAPHS ####
    
    # Generate plot of Cv vs T
    plt.figure()
    plt.plot(T,Cv)
    plt.title('Cv vs T')
    plt.xlabel('T (K)')
    plt.ylabel('Cv')

    plt.savefig("Cv_vs_T.pdf")
 
    # Generate plot of Q vs T
    plt.figure()
    plt.plot(T,Qavg)
    plt.title('Q vs T')
    plt.xlabel('T (K)')
    plt.ylabel('Q')
    plt.savefig("Q_vs_T.pdf")
    
