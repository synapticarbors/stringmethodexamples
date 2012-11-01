all: shared examples

examples: check-env DicksonRingPotential DicksonPeriodicPotential

shared:
	cd shared/pylangevin-integrator && python setup.py build_ext --inplace

DicksonPeriodicPotential:
	$(MAKE) -C examples/DicksonPeriodicPotential

DicksonRingPotential:
	$(MAKE) -C examples/DicksonRingPotential

clean:
	cd shared/pylangevin-integrator && rm -rf *.so build
	$(MAKE) -C examples/DicksonRingPotential clean

check-env:
	@if [ -z "$${WEST_ROOT}" ]; then echo "The env variable WEST_ROOT must be specified"  && exit 1; fi

.PHONY: examples shared check-env
