import os
import matplotlib.cm as cm
import time

from refactored_setup import *


if __name__ == '__main__':

    path = '/home/gsc4/scratch/SH3.constant.temp/test_melt/umbrella_sampling/4-12-2016_q_umbrella/'

    cwd = os.getcwd()

    os.chdir(path)

    t1 = time.time()

    # Read umbrella center names
    dir_name_format = '{}'
    umb_centers = [x.strip('\n') for x in open(str('umbrella_last'), 'r').readlines() ]
    dir_name_data = [umb_centers]
    umb_T = 130.95
    atom_pair = None

    print 'reading in umb data'
    E_unbias_concat, Q_concat, k_umb, n_frames, n_frame_toss = read_data_umb(dir_name_format, dir_name_data, 'SH3', 'qtanh.npy', umb_centers, atom_pair)
    
    print 'making mbar'
    mbar, umb_centers_interp = create_mbar_umb(umb_centers, E_unbias_concat, Q_concat, n_frames, k_umb, umb_T, n_interp=0)

    print 'initialization time', time.time() - t1

    print 'spaghetti math'
    Qi_vs_Q, bin_mid, loops = mbar_spaghetti(dir_name_format, dir_name_data, 'SH3', Q_concat, len(umb_centers_interp)+1, mbar, numbins=50, frame_toss=n_frame_toss)
    
    max_loop = np.max(loops)

    if not os.path.exists('mbar_QivsQ_umb'):
        os.mkdir('mbar_QivsQ_umb')
    os.chdir('mbar_QivsQ_umb')

    for i in range(len(umb_centers_interp) + 1):
        if i < len(umb_centers_interp):
            np.save('mbar_QivsQ_umb_' + str(umb_centers_interp[i]) + '.npy', Qi_vs_Q[i,:,:])
        else:
            np.save('mbar_QivsQ_umb_unbiased.npy', Qi_vs_Q[i,:,:])
    np.save('centers.npy', umb_centers_interp)
    np.save('loop.npy', loops)
    np.save('bin_mid.npy', bin_mid)

    os.chdir(cwd)




