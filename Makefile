all: parafly


parafly:
	cd plugins/ParaFly && ./configure --prefix=`pwd` && make && make install



large_sample_data:
	git clone https://github.com/EVidenceModeler/EVM_sample_data.git

