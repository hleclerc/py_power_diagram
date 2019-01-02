comp:
	test -e ext/power_diagram || git clone git@github.com:hleclerc/power_diagram.git ext/power_diagram
	test -e ext/xsimd || git clone https://github.com/quantstack/xsimd.git ext/xsimd
	python3 setup.py build
	python3 setup.py install --user 

	
test:
	python3 -m unittest discover tests

.PHONY: test comp
