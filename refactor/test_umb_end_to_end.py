import os
import matplotlib.cm as cm
import time

from refactored_setup import *


if __name__ == '__main__':

    path = '/home/gsc4/scratch/SH3.constant.temp/test_melt/umbrella_sampling/4-7-2016_end_end_umb/'

    cwd = os.getcwd()

    os.chdir(path)
    
    t1 = time.time()

    # Read umbrella center names
    dir_name_format = '{}'
    umb_centers = [x.strip('\n') for x in open(str('umbrella_last'), 'r').readlines() ]
    dir_name_data = [umb_centers]
    atom_pair = np.array([0,57])
    umb_T = 129.0

    print 'reading in umb data'
    E_unbias_concat, r1N_concat, k_umb, n_frames = read_data_umb(dir_name_format, dir_name_data, 'SH3', 'r1N.npy', umb_centers, atom_pair)

    print 'making mbar'
    mbar, umb_centers_interp = create_mbar_umb(umb_centers, E_unbias_concat, r1N_concat, n_frames, k_umb, umb_T, n_interp=0)

    print 'initialization time', time.time()-t1

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
    np.save('centers.npy', umb_centers_interp)
    np.save('loop.npy', loops)
    np.save('bin_mid.npy', bin_mid)

    os.chdir(cwd)

