import numpy as np
import h5py
import argparse
import sys

def run(fin_name, fout_name):
    fin = h5py.File(fin_name,'r')
    fout = h5py.File(fout_name,'a')
    nbins = 100

    # Get last iteration in input file
    niters = fin.attrs.get('wemd_current_iteration') - 1

    data_grp = fout.require_group('dist')
    if 'data' in data_grp:
        dshape = data_grp['data'].shape
        if dshape[0] < niters:
            data_grp['data'].resize((niters-2,nbins))
        
        data = data_grp['data']
        start_iter = data_grp.attrs['last_completed_iter']
    else:
        data = data_grp.require_dataset('data',(niters-2,nbins),np.float,exact=False,maxshape=(None,nbins))
        start_iter = 2
        data_grp.attrs['last_completed_iter'] = 2

    print 'starting iteration: {}'.format(start_iter)

    for iiter in xrange(start_iter,niters):
        if iiter % 1000 == 0:
            print 'Processing {} of {}'.format(iiter,niters-1)
            fout.flush()

        try:
            iter_grp = fin['iter_{:08d}'.format(iiter)]

            weight = iter_grp['seg_index']['weight']
            crd = iter_grp['pcoord'][:,-1,1]

            assert weight.shape[0] == crd.shape[0]

            # simulation time spent on iteration
            h,edges = np.histogram(crd,weights=weight,range=(0.0,1.0),bins=nbins)

            data[iiter-2,:] = h
            data_grp.attrs['last_completed_iter'] = iiter

        except:
            print 'Error in processing iteration: {}'.format(iiter)
            print sys.exc_info()
            break

    fin.close()
    fout.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='WEMD distribution analysis script')
    parser.add_argument('-f', dest='h5in', help='input h5 file')
    parser.add_argument('-o', dest='h5out', help='name of output file')

    args = parser.parse_args() 

    run(args.h5in,args.h5out)
