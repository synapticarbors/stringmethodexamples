all: shared examples

examples: check-env DicksonRingPotential DicksonPeriodicPotential ElasticNetworkModel

shared:
	cd shared/pylangevin-integrator && python setup.py build_ext --inplace
	cd shared/pyqcprot && python setup.py build_ext --inplace
	cd shared/elasticnetwork-langevin && python setup.py build_ext --inplace

DicksonPeriodicPotential:
	$(MAKE) -C examples/DicksonPeriodicPotential

DicksonRingPotential:
	$(MAKE) -C examples/DicksonRingPotential

ElasticNetworkModel:
	$(MAKE) -C examples/ElasticNetworkModel

clean:
	cd shared/pylangevin-integrator && rm -rf *.so build
	cd shared/pyqcprot && rm -rf *.so build
	cd shared/elasticnetwork-langevin && rm -rf *.so build
	$(MAKE) -C examples/DicksonRingPotential clean
	$(MAKE) -C examples/DicksonPeriodicPotential clean
	$(MAKE) -C examples/ElasticNetworkModel clean

check-env:
	@if [ -z "$${WEST_ROOT}" ]; then echo "The env variable WEST_ROOT must be specified"  && exit 1; fi

.PHONY: examples shared check-env
