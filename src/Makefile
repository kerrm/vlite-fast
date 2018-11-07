LPGPLOT = -L/usr/local/pgplot -lpgplot -lcpgplot -L/usr/X11R6/lib -lpng -lz -lX11

INCDIFXMESSAGE =  /home/vlite-master/frb-search/difxmessage-2.4.0
PARSE_FILES = eop.c  executor.c  options.c vlaant.c  vlite_xml.c multicast.c alert.c
PARSE_OBJ = eop.o  executor.o  options.o vlaant.o  vlite_xml.o multicast.o alert.o
PARSE_HEADERS = efaults.h  eop.h  executor.h  options.h  vlaant.h  vlite_xml.h multicast.h alert.h

#F77 = gfortran -std=legacy -ffixed-line-length-none 
CC = gcc
NVCC = nvcc
CUDA_ARCH = -gencode arch=compute_35,code=sm_35 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_61,code=sm_61
CUDA_DIR = /usr/cuda
PACKAGE_DIR = ../..

all: process_baseband genbase readbase messenger writer

process_baseband: process_baseband.o util.o cuda_util.o utils.o
	$(NVCC) -o process_baseband process_baseband.o utils.o util.o cuda_util.o -lcufft -lcurand -L${PACKAGE_DIR}/lib -lvdifio -lpsrdada

process_baseband.o: process_baseband.cu process_baseband.h
	$(NVCC) -g -O2 $(CUDA_ARCH) --compiler-options -Wall -c -o $@ -I${CUDA_DIR}/include -I${PACKAGE_DIR}/include process_baseband.cu 

genbase: genbase.cu util.o cuda_util.o
	$(NVCC) -g -O2 $(CUDA_ARCH) --compiler-options -Wall -I${CUDA_DIR}/include -I${CUDA_DIR}/samples/common/inc -I${PACKAGE_DIR}/include -c genbase.cu
	$(NVCC) -o genbase util.o cuda_util.o genbase.o -lcufft -lcurand -L${PACKAGE_DIR}/lib -lvdifio -lpsrdada

readbase:	readbase.c util.o
	g++ -g -O2 -Wall -c readbase.c -I${PACKAGE_DIR}/include
	g++ -o $@ readbase.o util.o -L${PACKAGE_DIR}/lib -lpsrdada
	chmod a+rx $@

messenger:	messenger.c def.h utils.o $(PARSE_FILES)
	gcc -g -O2 -Wall -std=gnu99 -c messenger.c $(PARSE_FILES) -I${PACKAGE_DIR}/include -I$(INCDIFXMESSAGE)
	gcc -o $@ messenger.o utils.o $(PARSE_OBJ) -L${PACKAGE_DIR}/lib -lexpat -lpsrdada
	chmod a+rx $@

writer:	writer.c def.h utils.o
	gcc -g -O2 -Wall -std=gnu99 -pthread -c writer.c -I${PACKAGE_DIR}/include -I$(INCDIFXMESSAGE)
	gcc -o $@ writer.o utils.o -pthread -L${PACKAGE_DIR}/lib -lpsrdada -lvdifio
	chmod a+rx $@

util.o: util.h util.c
	g++ -g -O2 -Wall -c util.c -I${PACKAGE_DIR}/include

cuda_util.o: cuda_util.h cuda_util.cu
	$(NVCC) -g -O2 --compiler-options -Wall -c -o $@ -I${CUDA_DIR}/include -I{${PACKAGE_DIR}/include cuda_util.cu

utils.o: utils.h utils.c Connection.h
	gcc -g -O2 -Wall -std=gnu99 -c -o $@ utils.c -I${PACKAGE_DIR}/include

clean:
	rm -f *.o core* process_baseband genbase readbase messenger writer


