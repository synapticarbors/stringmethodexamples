import numpy as np
import h5py
import os
import sys

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

try:
    import timeseries # from pymbar
except:
    sys.exit('This script requires the timeseries module from pymbar')

file_dir = os.path.dirname(os.path.abspath(__file__))
basedir = os.path.split(file_dir)[0]

# Font settings
### rcParams are the default parameters for matplotlib
mpl.rcParams['font.size'] = 12.
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['axes.labelsize'] = 10.
mpl.rcParams['xtick.labelsize'] = 10.
mpl.rcParams['ytick.labelsize'] = 10.

legfont = fm.FontProperties(size=10)
fsize = (3.375, 3.0)

data_file = os.path.join(basedir, 'we_freestr', 'bin_populations_all.h5')
stat_ineff_file = os.path.join(basedir, 'we_freestr', 'stat_inefficiency.np')
skip = 500


def get_probs(fname):
    f = h5py.File(data_file,'r')
    bp = f['data'][...]
    f.close()

    return bp

if __name__ == '__main__':

    bp = get_probs(data_file)[skip:,:]
    n_samples, n_bins = bp.shape
    bi = np.arange(n_bins)

    pmf = -0.6*np.log(bp)
    pmf_mean = pmf.mean(axis=0)
    pmf_mean -= np.min(pmf_mean)

    pmf_std = pmf.std(axis=0)

    # Calculate statistical inefficiency
    try:
        g = np.load(stat_ineff_file)
    except:
        g = np.zeros((n_bins,))
        for k in xrange(n_bins):
            g[k] = timeseries.statisticalInefficiency(pmf[:,k])

        np.save(stat_ineff_file, g)

    pmf_err = 1.96*g*pmf_std/np.sqrt(n_samples)

    offset = np.min(pmf_mean - pmf_err)
    pmf_mean -= offset

    fig = plt.figure(1, figsize=fsize)
    ax = fig.add_subplot(111)
    ax.plot(bi, pmf_mean, color='black', lw=2, zorder=1000)
    ax.fill_between(bi, pmf_mean + pmf_err, pmf_mean - pmf_err, alpha=.4, facecolor='gray', edgecolor='black', lw=1, zorder=900)
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$G_{\alpha}$ (kcal/mol)')

    fig.set_size_inches(fsize)
    plt.tight_layout()
    eps_file = os.path.join(basedir, 'generate_figures', 'ntrc_pmf.eps')
    plt.savefig(eps_file, dpi=600, format='eps', bbox_inches='tight')
    #plt.show()


