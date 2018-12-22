all:
	git submodule update --init --recursive
	python3 setup.py install --user 

	