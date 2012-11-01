import numpy as np
import h5py
import argparse
import sys

def run(fin_name, fout_name):
    fin = h5py.File(fin_name,'r')
    fout = h5py.File(fout_name,'a')
    # Get last iteration in input file
    niters = fin.attrs.get('wemd_current_iteration') - 1

    print 'Getting pcoord data ({} iterations)'.format(niters)

    data_grp = fout.require_group('pcoord')
    for iiter in xrange(2,niters,500):
        print 'Processing {} of {}'.format(iiter,niters-1)
        fout.flush()

        try:
            iter_grp = fin['iter_{:08d}'.format(iiter)]
            crd = iter_grp['pcoord'][:,-1,:]

        except:
            print 'Error in processing iteration: {}'.format(iiter)
            print sys.exc_info()
            break

        data_iter_grp = data_grp.require_group('iter_{:08d}'.format(iiter))

        crd_ds = data_iter_grp.require_dataset('pcoord',shape=crd.shape,dtype=np.float64)
        crd_ds[...] = crd

    fin.close()
    fout.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get pcoords for representative snapshots')
    parser.add_argument('-f', dest='h5in', help='input h5 file')
    parser.add_argument('-o', dest='h5out', help='name of output file')

    args = parser.parse_args() 

    run(args.h5in,args.h5out)

