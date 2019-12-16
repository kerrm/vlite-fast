LPGPLOT = -L/usr/local/pgplot -lpgplot -lcpgplot -L/usr/X11R6/lib -lpng -lz -lX11

INCDIFXMESSAGE =  /home/vlite-master/frb-search/difxmessage-2.4.0
PARSE_FILES = eop.c  executor.c  options.c vlaant.c  vlite_xml.c multicast.c alert.c
PARSE_OBJ = eop.o  executor.o  options.o vlaant.o  vlite_xml.o multicast.o alert.o
PARSE_HEADERS = defaults.h  eop.h  executor.h  options.h  vlaant.h  vlite_xml.h multicast.h alert.h

#F77 = gfortran -std=legacy -ffixed-line-length-none 
CC = gcc
#NVCC = nvcc
CUDA_ARCH = -gencode arch=compute_35,code=sm_35 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_61,code=sm_61
#CUDA_DIR = /usr/cuda
CUDA_DIR = /usr/cuda-10.1.243.418.87.00
NVCC = $(CUDA_DIR)/bin/nvcc
MAKE_PACKAGE_DIR = /home/vlite-master/mtk
#PACKAGE_DIR = /home/vlite-master/surya
PACKAGE_DIR = /home/vlite-master/mtk

all: process_baseband genbase readbase messenger writer

install: all
	cp -p process_baseband ${PACKAGE_DIR}/bin
	cp -p genbase ${PACKAGE_DIR}/bin
	cp -p readbase ${PACKAGE_DIR}/bin
	cp -p messenger ${PACKAGE_DIR}/bin
	cp -p writer ${PACKAGE_DIR}/bin
	${PACKAGE_DIR}/vlite-fast/scripts/install_binaries

process_baseband: process_baseband.o util.o cuda_util.o utils.o multicast.o pb_kernels.o
	$(NVCC) -o $@ $^ -L${CUDA_DIR}/lib64  -lcufft -lcurand -L${MAKE_PACKAGE_DIR}/lib -lvdifio -lpsrdada
	chmod a+rx $@

process_baseband.o: process_baseband.cu pb_kernels.cu process_baseband.h def.h
	$(NVCC) -g -O2 $(CUDA_ARCH) --compiler-options -Wall -c -o $@ -I${CUDA_DIR}/include -I${MAKE_PACKAGE_DIR}/include process_baseband.cu

pb_kernels.o: pb_kernels.cu process_baseband.h def.h
	$(NVCC) -g -O2 $(CUDA_ARCH) --compiler-options -Wall -c -o $@ -I${CUDA_DIR}/include -I${MAKE_PACKAGE_DIR}/include pb_kernels.cu

genbase: genbase.cu util.o cuda_util.o
	$(NVCC) -g -O2 $(CUDA_ARCH) --compiler-options -Wall -I${CUDA_DIR}/include -I${CUDA_DIR}/samples/common/inc -I${MAKE_PACKAGE_DIR}/include -c genbase.cu
	$(NVCC) -o genbase util.o cuda_util.o genbase.o -lcufft -lcurand -L${MAKE_PACKAGE_DIR}/lib -lvdifio -lpsrdada
	chmod a+rx $@

readbase:	readbase.c util.o
	g++ -g -O2 -Wall -c readbase.c -I${MAKE_PACKAGE_DIR}/include
	g++ -o $@ readbase.o util.o -L${MAKE_PACKAGE_DIR}/lib -lpsrdada
	chmod a+rx $@

messenger:	messenger.c def.h utils.o $(PARSE_FILES)
	gcc -g -O2 -Wall -std=gnu99 -c messenger.c $(PARSE_FILES) -I${MAKE_PACKAGE_DIR}/include -I$(INCDIFXMESSAGE)
	gcc -o $@ messenger.o utils.o $(PARSE_OBJ) -L${MAKE_PACKAGE_DIR}/lib -lexpat -lpsrdada
	chmod a+rx $@

writer:	writer.c def.h utils.o $(PARSE_FILES)
	gcc -g -O2 -Wall -std=gnu99 -pthread -c writer.c $(PARSE_FILES) -I${MAKE_PACKAGE_DIR}/include -I$(INCDIFXMESSAGE)
	gcc -o $@ writer.o utils.o $(PARSE_OBJ) -pthread -L${MAKE_PACKAGE_DIR}/lib -lpsrdada -lvdifio -lexpat -Wl,-rpath /home/vlite-master/mtk/lib
	chmod a+rx $@

util.o: util.h util.c
	g++ -g -O2 -Wall -c util.c -I${MAKE_PACKAGE_DIR}/include

cuda_util.o: cuda_util.h cuda_util.cu
	$(NVCC) -g -O2 --compiler-options -Wall -c -o $@ -I${CUDA_DIR}/include -I{${MAKE_PACKAGE_DIR}/include cuda_util.cu

utils.o: utils.h utils.c Connection.h
	gcc -g -O2 -Wall -std=gnu99 -c -o $@ utils.c -I${MAKE_PACKAGE_DIR}/include

clean:
	rm -f *.o core* process_baseband genbase readbase messenger writer


