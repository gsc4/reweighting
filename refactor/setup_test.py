import os
import matplotlib.cm as cm
import time

from refactored_setup import *

if __name__ == '__main__':

#    t1 = time.time()

    path = '/home/gsc4/scratch/SH3.constant.temp/test_melt/umbrella_sampling/5-27-2016_sh3_ca_sbm/'

    cwd = os.getcwd()

    os.chdir(path)

##### CONSTANT TEMP #####
    """
    name_format = 'T_{}_{}/{}'
    name_data = [[128.5, 129.0, 129.5], [0,1,2],['Etot.xvg']]
    E_concat, n_frames = read_data_const_temp(name_format, name_data, n_trials=3)

    name_data = [[128.5,129.0,129.5], [0,1,2],['qtanh.dat']] 
    name_format = 'T_{}_{}/{}'
    Q_concat, dumbo = read_data_const_temp(name_format, name_data, n_trials=3)

    temps = [128.5, 129.0, 129.5]
    mbar, T_interp = create_mbar_const_temp(temps, E_concat, n_frames, n_interp=0)
    #Cv = calc_Cv(mbar, E_concat, T_interp)
    """
    
#    print time.time() - t1 

#### MBAR SPAGHETTI ####
    """
    dir_name_format = 'T_{}_{}'
    dir_name_data = [temps, range(3)]
    mbar_qivsQ, bin_mid, loops = mbar_spaghetti(dir_name_format, dir_name_data, 'SH3', Q_concat, len(T_interp), mbar, numbins = 50) 

    if not os.path.exists('mbar_QivsQ'):
        os.mkdir('mbar_QivsQ')
    os.chdir('mbar_QivsQ')
    for i in range(len(T_interp)):
        np.save('mbar_QivsQ_' + str(T_interp[i]) + '.npy', mbar_qivsQ[i,:,:])
        
    np.savetxt('temps', T_interp)
    np.savetxt('loop', loops)
    np.savetxt('bin_mid', bin_mid)
    """

#### EMPIRICAL SPAGHETTI    

    #### Qi vs Q ####
    """
    #dir_name_format = 'T_{}_{}'
    #dir_name_data = [[129.0], [0,1,2]]
    #qivsQ, bin_mid, loops = empirical_spaghetti(dir_name_format, dir_name_data, 'SH3', Q_concat)
    
    max_loop = np.max(loops)
    
    for i in range(len(qivsQ)):
        plt.plot(bin_mid, qivsQ[i,:], c=cm.jet(loops[i]/max_loop)) 
 
    plt.title('qi vs Qtot')
    plt.xlabel('Qtot')
    plt.ylabel('Percent formation of contact')
    
    plt.savefig('qivsQ.pdf')
    """

#### MBAR SPAGHETTI ####
    
    #### Qi vs Q ####
    """
    max_loop = np.max(loops)    

    for i in range(len(T_interp)):
        
        plt.figure()
        
        for j in range(len(mbar_qivsQ)):
            plt.plot(bin_mid, mbar_qivsQ[i,j,:], c=cm.jet(loops[i]/max_loop))
        
        plt.title('qi vs Qtot')
        plt.xlabel('Qtot')
        plt.ylabel('Percent formation of contact')
        plt.savefig('QivsQ_' + str(T_interp[i]))
        plt.close()       
    """

    #### Qi vs r1N ####
    path = '/home/gsc4/scratch/SH3.constant.temp/test_melt/umbrella_sampling/4-7-2016_end_end_umb/'
    
    os.chdir(path)

    t1 = time.time()

    dir_name_format = '{}'
    umb_centers = [x.strip('\n') for x in open(str('umbrella_last'), 'r').readlines() ]        
    dir_name_data = [umb_centers]
    atom_pair = np.array([0,57])
    umb_T = 129.

    print 'reading in umb data'
    E_unbias_concat, r1N_concat, k_umb, n_frames = read_data_umb(dir_name_format, dir_name_data, 'SH3', 'r1N.npy', umb_centers, atom_pair)
    
    print 'making mbar'    
    mbar, umb_centers_interp = create_mbar_umb(umb_centers, E_unbias_concat, r1N_concat, n_frames, k_umb, umb_T, n_interp=1)

    print time.time()-t1

    print 'spaghetti math'
    Qi_vs_r1N, bin_mid, loops = mbar_spaghetti(dir_name_format, dir_name_data, 'SH3', r1N_concat, len(umb_centers_interp)+1, mbar, numbins=50)

    max_loop = np.max(loops)
    
    if not os.path.exists('mbar_Qivsr1N'):
        os.mkdir('mbar_Qivsr1N')
    os.chdir('mbar_Qivsr1N')
    
    for i in range(len(umb_centers_interp) + 1):
        if i < len(umb_centers_interp):
            np.save('mbar_Qivsr1N_' + str(umb_centers_interp[i]) + '.npy', Qi_vs_r1N[i,:,:])
        else:
            np.save('mbar_Qivsr1N_unbiased.npy', Qi_vs_r1N[i,:,:]) 
    np.savetxt('centers', umb_centers_interp)
    np.savetxt('loop', loops)
    np.savetxt('bin_mid', bin_mid)
    #for i in range(len(qivsr1N)

#    t2 = time.time()

#    t_total = t2-t1

#    print t_total

##### PLOT CONST TEMP  #####
    """    
    ### argsort to reorder
    order = np.argsort(T_interp)

    plt.figure()
    plt.plot(T_interp[order],Cv[order])
    plt.title('Cv vs T')
    plt.xlabel('T (K)')
    plt.ylabel('Cv')

    plt.savefig("Cv_vs_T.pdf")

    """
