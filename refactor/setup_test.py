import os

from refactored_setup import *

if __name__ == '__main__':

    path = '/home/gsc4/scratch/SH3.constant.temp/test_melt/umbrella_sampling/5-27-2016_sh3_ca_sbm/'

    cwd = os.getcwd()

    os.chdir(path)

    ##### Constant temp test #####

    name_format = 'T_{}_{}/{}'
    name_data = [[128.5, 129.0, 129.5], [0,1,2],['Etot.xvg']]

    #E_concat, n_frames = read_data_const_temp(name_format, name_data, n_trials=3)

    name_data = [[129.0], [0,1,2],['qtanh.dat']]

    Q_concat, dumbo = read_data_const_temp(name_format, name_data, n_trials=3)

    temps = [128.5, 129.0, 129.5]

    #mbar, T_interp = create_mbar_const_temp(temps, E_concat, n_frames)

    #Cv = calc_Cv(mbar, E_concat, T_interp)

    dir_name_format = 'T_{}_{}'
    dir_name_data = [[129.0], [0,1,2]]

    empirical_spaghetti(dir_name_format, dir_name_data, 'SH3', Q_concat, 3)

    os.chdir(cwd)
    
    ##### Plot temperature stuff #####

######## This seems like a super circuitious route to do this... ####
    T_Cv_zip = zip(T_interp,Cv)
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



