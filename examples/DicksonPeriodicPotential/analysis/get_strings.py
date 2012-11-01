import numpy as np
import h5py
import argparse
import sys

def run(fin_name, fout_name):
    fin = h5py.File(fin_name,'r')
    fout = h5py.File(fout_name,'a')

    # Get last iteration of string update
    niters = fin['stringmethod'].attrs['last_update'] + 2

    data_grp = fout.require_group('string')

    for iiter in xrange(2,niters):
        if iiter % 100 == 0:
            print 'Processing {} of {}'.format(iiter,niters-1)
            fout.flush()

        avgpos_flag = False

        try:
            iter_grp = fin['iter_{:08d}'.format(iiter)]
            
            centers = iter_grp['stringmethod']['centers'][:]

            if 'avgpos' in iter_grp['stringmethod']:
                avgpos = iter_grp['stringmethod']['avgpos'][:]
                avgpos_flag = True

        except:
            print 'Error in processing iteration: {}'.format(iiter)
            print sys.exc_info()
            break

        data_iter_grp = data_grp.require_group('iter_{:08d}'.format(iiter))

        centers_ds = data_iter_grp.require_dataset('centers',shape=centers.shape,dtype=np.float64)
        centers_ds[...] = centers

        if avgpos_flag:
            avg_pos_ds = data_iter_grp.require_dataset('avgpos',shape=avgpos.shape,dtype=np.float64)
            avg_pos_ds[...] = avgpos

    fin.close()
    fout.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get string coordinates and average positions if present')
    parser.add_argument('-f', dest='h5in', help='input h5 file')
    parser.add_argument('-o', dest='h5out', help='name of output file')

    args = parser.parse_args() 

    run(args.h5in,args.h5out)

