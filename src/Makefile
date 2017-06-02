LPGPLOT = -L/usr/local/pgplot -lpgplot -lcpgplot -L/usr/X11R6/lib -lpng -lz -lX11

INCDIFXMESSAGE =  /home/vlite-master/frb-search/difxmessage-2.4.0
PARSE_FILES = eop.c  executor.c  options.c vlaant.c  vlite_xml.c multicast.c alert.c
PARSE_OBJ = eop.o  executor.o  options.o vlaant.o  vlite_xml.o multicast.o alert.o
PARSE_HEADERS = efaults.h  eop.h  executor.h  options.h  vlaant.h  vlite_xml.h multicast.h alert.h

#F77 = gfortran -std=legacy -ffixed-line-length-none 
CC = gcc
NVCC = nvcc

all: process_baseband genbase readbase messenger writer

process_baseband: process_baseband.o util.o cuda_util.o utils.o
	$(NVCC) -o process_baseband process_baseband.o utils.o util.o cuda_util.o -lcufft -lcurand -L../lib -lvdifio -lpsrdada

process_baseband.o: process_baseband.cu
	$(NVCC) -g -O2 --compiler-options -Wall -c -o $@ -I/usr/cuda-6.5/include -I../include process_baseband.cu 

genbase: genbase.cu util.o cuda_util.o
	$(NVCC) -g -O2 --compiler-options -Wall -I/usr/cuda-6.5/include -I/usr/cuda-6.5/samples/common/inc -I/users/mkerr/mtk/include -c genbase.cu
	$(NVCC) -o genbase util.o cuda_util.o genbase.o -lcufft -lcurand -L../lib -lvdifio -lpsrdada

readbase:	readbase.c util.o
	g++ -g -O2 -Wall -c readbase.c -I../include
	g++ -o $@ readbase.o util.o -L../lib -lpsrdada
	chmod a+rx $@

messenger:	messenger.c def.h utils.o $(PARSE_FILES)
	gcc -g -O2 -Wall -std=gnu99 -c messenger.c $(PARSE_FILES) -I../include -I$(INCDIFXMESSAGE)
	gcc -o $@ messenger.o utils.o $(PARSE_OBJ) -L../lib -lexpat -lpsrdada
	chmod a+rx $@

writer:	writer.c def.h utils.o
	gcc -g -O2 -Wall -std=gnu99 -c writer.c -I../include -I$(INCDIFXMESSAGE)
	gcc -o $@ writer.o utils.o -L../lib -lpsrdada
	chmod a+rx $@

util.o: util.h util.c
	g++ -g -O2 -Wall -c util.c -I../include

cuda_util.o: cuda_util.h cuda_util.cu
	$(NVCC) -g -O2 --compiler-options -Wall -c -o $@ -I/usr/cuda-6.5/include -I../include cuda_util.cu

utils.o: utils.h utils.c Connection.h
	gcc -g -O2 -Wall -std=gnu99 -c -o $@ utils.c -I../include

clean:
	rm -f *.o core* process_baseband genbase readbase messenger writer

dataclean:
	rm -f ../data/2017*.fil

