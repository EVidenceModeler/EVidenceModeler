all: parafly


parafly:
	cd plugins/ParaFly && ./configure --prefix=`pwd` && make && make install
