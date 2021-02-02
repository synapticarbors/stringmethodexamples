all: shared examples

examples: DicksonRingPotential DicksonPeriodicPotential ElasticNetworkModel

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

.PHONY: examples shared
