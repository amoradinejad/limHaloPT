# NOTE: Depending on your installation paths for gsl and cuba libraries, you should add the links 
# for the INCLUDES and LIBS below, replacing -I/opt/local/include and  -L/opt/local/lib which are the 
# paths on my machine. 

CC        = gcc -c
CLINK     = gcc 

INCLUDES  =  -I../Include/  
INCLUDES +=  -I/opt/local/include
INCLUDES +=  -I../Class/include  

LIBS      =  -L/opt/local/lib  
LIBS     +=  -L../Class/lib
LIBS     +=  -L../lib

CFLAGS    =  -Wall
CFLAGS   += -fopenmp 
CFLAGS   += -g 

LIBFLAGS  = -lgsl 
LIBFLAGS += -lgslcblas 
LIBFLAGS += -lcuba  
LIBFLAGS += -lclass 
LIBFLAGS += -lm
LIBFLAGS += -llimHaloPT


test: test.o
	$(CLINK) test.o $(INCLUDES) $(LIBS) $(CFLAGS) $(LIBFLAGS) -o test 

test.o: test.c 
	$(CC) $(INCLUDES) $(CFLAGS) test.c -o test.o


clean:
	rm test.o test

