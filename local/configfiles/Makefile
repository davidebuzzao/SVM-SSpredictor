######################
##	CONFIGURATION	##
######################

SHELL := /bin/bash

all: configuration
	echo 'Done'

configuration:
	echo 'Writing symbolic links'
	for i in $$(ls ../src/); do ln -s ../src/$$i ../bin/$$(echo $$i | cut -d '.' -f '1'); done
	touch ../bin/__init__.py 