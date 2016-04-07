import numpy as np
import pymbar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

if __name__ == "__main__":
    
    #### READ FILES, SET CONSTANTS, INITIALIZE MBAR ####

    # Set boltzmann constant
    kb = .0083145 #kJ/mol*K
    
    
    print 'reading in temps'
    # Temps being considered
    temps = [ x.rstrip("\n") for x in open("short_temps_about_melt", "r").readlines() ]
   
    # Only take every 100th point so python can handle the data 
    stride = 100
    
    print 'reading in Q.dat files'
    # Read Q files
    Q = np.array([np.loadtxt("{}_0/Q.dat".format(x.rstrip("\n")))[::stride] for x in temps ])
    Q_long = np.concatenate(Q)
   
    # Make my histogram thingy 
    numbins = 100

    (bin_n , bin_edges) = np.histogram(Q_long, bins = numbins)
    bin_p = bin_n/len(Q_long)
    bin_mid= .5*(bin_edges[1:] + bin_edges[:-1])
 
     # Turn list of temperatures into array of floats with a set step size
    step_size = .05
    T = np.arange(float(temps[0]),float(temps[len(temps)-1])+step_size,step_size)

    print 'reading in energies'
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

    u_kn = np.zeros((len(T),len(Etot)))

    # for loop to fill in u_kn
    for i in range(len(T)):
        B = 1/(kb*T[i])
        u_kn[i] = Etot*B
    
    print 'initializing mbar'
    # Initialize mbar
    mbar = pymbar.MBAR(u_kn, N_k)

    p_of_Q = np.zeros((len(T),numbins),float)
    p_of_Q_min_0 = np.zeros((len(T),numbins),float)
    dp_of_Q = np.zeros((len(T),numbins),float)

    print 'calculating Qavg'
    for i in range(numbins):
    #for i in [10]:
        h = (Q_long > bin_edges[i]) & (Q_long < bin_edges[i+1])
        p_of_Q[:,i], dp_of_Q[:,i] = mbar.computeExpectations(h.astype(int))
        p_of_Q_min_0[:,i] = p_of_Q[:,i] - min(p_of_Q[:,i])
    
    A_over_RT = -np.log(p_of_Q)
    A_over_RT_min_0 = np.zeros(A_over_RT.shape)
    for i in range(numbins):
        A_over_RT_min_0[:,i] = A_over_RT[:,i] - min(A_over_RT_min_0[:,i])    

    blue = np.array([0, 0, 1])
    red = np.array([1, 0, 0])
    num_colors = len(T)/(5)+1
    color_step_size = 1./num_colors
    rgb_step = color_step_size*(red-blue)
    
    color = np.zeros((9,3))
    for i in range(num_colors):
        color[i] = blue + i*rgb_step

    plt.figure()
    temps_considered = []
    for i in range(0,len(T), 5 ):
        plt.plot(bin_mid, p_of_Q[i,:],c=color[i/5])
        temps_considered.append('T = ' + str(T[i]) + ' K')        

    plt.legend(temps_considered,loc=0)
    plt.title('p(Q) vs Q')
    plt.ylabel('p(Q)')
    plt.xlabel('Q')
    plt.savefig("p_of_Q_vs_Q_shifted_down.pdf") 

    plt.figure()
    for i in range(0,len(T),5):
        plt.plot(bin_mid, A_over_RT_min_0[i,:],c=color[i/5])

    plt.ylim(0, 10)
    plt.legend(temps_considered,loc=0)
    plt.title('A/RT vs Q')
    plt.xlabel('Q')
    plt.ylabel('A/RT')
    plt.savefig("A_over_RT_vs_Q_shifted_down.pdf")

    A_over_RT = -np.log(p_of_Q)
    A_over_RT_min_0 = np.zeros(A_over_RT.shape)
    for i in range(numbins):
        A_over_RT_min_0[:,i] = A_over_RT[:,i] - min(A_over_RT_min_0[:,i])

    blue = np.array([0, 0, 1])
    red = np.array([1, 0, 0])
    num_colors = len(T)/(5)+1
    color_step_size = 1./num_colors
    rgb_step = color_step_size*(red-blue)

    color = np.zeros((9,3))
    for i in range(num_colors):
        color[i] = blue + i*rgb_step

    plt.figure()
    temps_considered = []
    for i in range(0,len(T), 5 ):
        plt.plot(bin_mid, p_of_Q[i,:],c=color[i/5])
        temps_considered.append('T = ' + str(T[i]) + ' K')

    plt.legend(temps_considered,loc=9,prop={'size':10})
    plt.title('p(Q) vs Q')
    plt.ylabel('p(Q)')
    plt.xlabel('Q')
    plt.savefig("p_of_Q_vs_Q_shifted_down.pdf")

    plt.figure()
    for i in range(0,len(T),5):
        plt.plot(bin_mid, A_over_RT_min_0[i,:],c=color[i/5])

    plt.ylim(0, 10)
    plt.legend(temps_considered,loc=8,prop ={'size':10})
    plt.title('A/RT vs Q')
    plt.xlabel('Q')
    plt.ylabel('A/RT')
    plt.savefig("A_over_RT_vs_Q_shifted_down.pdf")

