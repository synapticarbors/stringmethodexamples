PYLANG_TARGET := $(shell python -c 'import os; import sys; print(os.path.realpath("../../shared/pylangevin-integrator"))')
SRC_TARGET := $(shell python -c 'import os; import sys; print(os.path.realpath("src"))')

all: clean setup

setup: 
	cd src && python setup.py build_ext --inplace
	for file in $(PYLANG_TARGET)/*.so; do ln -s $${file} we_base/$$(basename -- $${file}); done
	for file in $(PYLANG_TARGET)/*.so; do ln -s $${file} bruteforce_base/$$(basename -- $${file}); done
	for file in $(SRC_TARGET)/*.so; do ln -s $${file} we_base/$$(basename -- $${file}); done
	for file in $(SRC_TARGET)/*.so; do ln -s $${file} analysis/$$(basename -- $${file}); done
	for file in $(SRC_TARGET)/*.so; do ln -s $${file} generate_figures/$$(basename -- $${file}); done
	touch we_base/env.sh
	echo 'export WM_WORK_MANAGER=serial' >> we_base/env.sh

clean:
	rm -f we_base/env.sh we_base/*.so bruteforce_base/*.so analysis/*.so generate_figures/*.so src/*.so

.PHONY: all clean setup 

