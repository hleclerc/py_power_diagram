comp:
	git submodule update --init --recursive
	python3 setup.py install --user 

	
test:
	python3 -m unittest discover tests

.PHONY: test comp
