all: parafly

OS := $(shell uname)

CXX = g++
CC = gcc

parafly:
	cd plugins/ParaFly && sh ./configure --prefix=`pwd` CXX=$(CXX) CC=$(CC) CFLAGS="-fopenmp" CXXFLAGS="-fopenmp" && $(MAKE) install 



large_sample_data:
	git clone https://github.com/EVidenceModeler/EVM_sample_data.git

