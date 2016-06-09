import os
import matplotlib.cm as cm
import time

from refactored_setup import *

if __name__ == '__main__':

    path = '/home/gsc4/scratch/SH3.constant.temp/test_melt/umbrella_sampling/5-27-2016_sh3_ca_sbm/'

    cwd = os.getcwd()

    os.chdir(path)

###################################################################################################

    # Read in energy data
    print 'reading in energy data'    
    name_format = 'T_{}_{}/{}'
    name_data = [[128.5, 129.0, 129.5], [0,1,2],['Etot.xvg']]
    E_concat, n_frames = read_data_const_temp(name_format, name_data, n_trials=3)

    # Initialize MBAR
    print 'Initializing mbar'
    temps = [128.5, 129.0, 129.5]
    mbar, T_interp = create_mbar_const_temp(temps, E_concat, n_frames, n_interp=0)
    
    # Calculating Cv
    print 'Calculating Cv with MBAR'
    Cv = calc_Cv(mbar, E_concat, T_interp)

###################################################################################################

    # Empirical Spaghetti    
    print 'reading in Q data'
    name_data = [[129.0], [0,1,2],['qtanh.dat']]
    name_format = 'T_{}_{}/{}'
    Q_concat, dummy = read_data_const_temp(name_format, name_data, n_trials=3)

    print 'Calculating Empirical Spaghetti data'
    dir_name_format = 'T_{}_{}'
    dir_name_data = [[129.0], [0,1,2]]
    qivsQ, bin_mid, loops = empirical_spaghetti(dir_name_format, dir_name_data, 'SH3', Q_concat)

    print 'Plotting Empirical Spaghetti'
    max_loop = np.max(loops)
    for i in range(len(qivsQ)):
        plt.plot(bin_mid, qivsQ[i,:], c=cm.jet(loops[i]/max_loop))
    plt.title('qi vs Qtot')
    plt.xlabel('Qtot')
    plt.ylabel('Percent formation of contact')
    plt.savefig('test_empirical_qivsQ.pdf')

###################################################################################################

    # MBAR Spaghetti
    
    print 'reading in energy data'    
    name_format = 'T_{}_{}/{}'
    name_data = [[128.5, 129.0, 129.5], [0,1,2],['Etot.xvg']]
    E_concat, n_frames = read_data_const_temp(name_format, name_data, n_trials=3)

    print 'Initializing mbar'
    temps = [128.5, 129.0, 129.5]
    mbar, T_interp = create_mbar_const_temp(temps, E_concat, n_frames, n_interp=0)
    
    print 'reading in Q data'
    name_data = [[128.5,129.0,129.5], [0,1,2],['qtanh.dat']]
    name_format = 'T_{}_{}/{}'
    Q_concat, dumbo = read_data_const_temp(name_format, name_data, n_trials=3)

    print 'Calculating MBAR Spaghetti Data'
    dir_name_format = 'T_{}_{}'
    dir_name_data = [temps, range(3)]
    mbar_qivsQ, bin_mid, loops = mbar_spaghetti(dir_name_format, dir_name_data, 'SH3', Q_concat, len(T_interp), mbar, numbins = 50)

    print 'Plotting MBAR Spaghetti'
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


