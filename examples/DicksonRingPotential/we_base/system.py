from __future__ import division, print_function
__metaclass__ = type
import time
import os
import numpy as np

import west
from west.propagators import WESTPropagator
from west import Segment, WESTSystem
from westpa.binning import VoronoiBinMapper

import cIntegratorSimple
import ForceFields
from utils import dfunc

import logging
log = logging.getLogger(__name__)
log.debug('loading module %r' % __name__)

pcoord_dtype = np.float32


def genrandint():
    """Generates a random integer between 0 and (2^32)-1"""
    x = 0
    for i in range(4):
        x = (x << 8) + ord(os.urandom(1))

    return x


class SimpleLangevinPropagator(WESTPropagator):

    def __init__(self, rc=None):
        super(SimpleLangevinPropagator, self).__init__(rc)

        rc = self.rc.config['west','simplelangevin']
        self.ndim = rc.get('ndim', 2)
        self.nsteps = rc.get('blocks_per_iteration', 2)
        self.nsubsteps = rc.get('steps_per_block')
        self.beta = rc.get('beta')

        ff = ForceFields.Dickson2dRingForce()
        MASS = 1.0
        XI = 1.5
        BETA = self.beta
        NDIMS = 2
        DT = 0.005
        ISPERIODIC = np.array([0, 0], dtype=np.int)
        BOXSIZE = np.array([1.0E8, 1.0E8], dtype=pcoord_dtype)

        self.integrator = cIntegratorSimple.SimpleIntegrator(ff, MASS, XI, BETA, DT, NDIMS, ISPERIODIC, BOXSIZE, genrandint())

    def get_label(self, x, last_state):
        if (x[0] + 3.0)**2 + x[1]**2 < 1.0:
            state = 0
        elif (x[0] - 3.0)**2 + x[1]**2 < 1.0:
            state = 1
        else:
            state = last_state

        return state

    def get_pcoord(self, state):
        pcoord = None
        if state.label == 'initA':
            pcoord = [-3.0, 0.0, 0]
        elif state.label == 'initI1':
            pcoord = [0.0, 3.0, 0]
        elif state.label == 'initB':
            pcoord = [3.0, 0.0, 1]
        elif state.label == 'initI2':
            pcoord = [0.0, -3.0, 1]

        state.pcoord = pcoord

    def propagate(self, segments):

        for segment in segments:
            starttime = time.time()
            new_pcoords = np.empty((self.nsteps, self.ndim+1), dtype=pcoord_dtype)
            new_pcoords[0,:] = segment.pcoord[0,:]
            x = new_pcoords[0,:2].copy()

            for istep in xrange(1, self.nsteps):
                self.integrator.step(x, self.nsubsteps)
                new_pcoords[istep,:2] = x
                new_pcoords[istep,2] = self.get_label(x, new_pcoords[istep-1,2])

            segment.pcoord = new_pcoords[...]
            segment.status = Segment.SEG_STATUS_COMPLETE
            segment.walltime = time.time() - starttime

        return segments


class System(WESTSystem):

    def initialize(self):

        rc = self.rc.config['west', 'system']

        self.pcoord_ndim = 3
        self.pcoord_len = 2
        self.pcoord_dtype = pcoord_dtype
        self.target_count = rc.get('target_count')
        self.nbins = rc.get('nbins')

        slen = self.nbins // 2
        x = np.linspace(-3.0, 3.0, slen)
        centers = np.zeros((self.nbins, self.pcoord_ndim), dtype=self.pcoord_dtype)
        centers[:slen, 0] = x
        centers[slen:, 0] = x[::-1]

        centers[:slen, 2] = 0.0
        centers[slen:, 2] = 1.0

        self.bin_mapper = VoronoiBinMapper(dfunc, centers)
        self.bin_target_counts = np.zeros((self.bin_mapper.nbins,), dtype=np.int_)
        self.bin_target_counts[...] = self.target_count

        slen = self.nbins // 2
        self.sm_params = {'slen': [slen, slen],
                          'kappa': 0.001,
                          'dtau': 0.15,
                          'fixed_ends': False,
                          'sciflag': True,
                          'mpairs': [[0, self.nbins - 1], [slen - 1, slen]],
                          'slabels': [2],
                          'fourierflag': True,
                          'fourier_P': 2}
