# NOTE: Depending on your installation paths for gsl and cuba libraries, you should add the links 
# for the INCLUDES and LIBS below, replacing -I/opt/local/include and  -L/opt/local/lib which are the 
# paths on my machine. 

CC        = gcc -c
CLINK     = gcc 

INCLUDES  =  -I../Include/  
INCLUDES += -I/opt/local/include
INCLUDES +=  -I../Class/include  

LIBS      =  -L/opt/local/lib  
LIBS     +=  -L../Class/lib

CFLAGS   = -Wall
CFLAGS   += -fopenmp 
CFLAGS   += -g 

LIBFLAGS  =  -lgsl 
LIBFLAGS += -lgslcblas 
LIBFLAGS += -lcuba  
LIBFLAGS += -lclass 
LIBFLAGS += -lm

SOURCE_FILES = $(wildcard *.c)
$(info $$SOURCE_FILES = $(SOURCE_FILES))

OBJECT_FILES = $(patsubst %.c, %.o, $(notdir $(SOURCE_FILES)))
$(info $$OBJECT_FILES = $(OBJECT_FILES))

limHaloPT: ${OBJECT_FILES}
	$(CLINK) ${OBJECT_FILES} $(INCLUDES) $(LIBS) $(CFLAGS) $(LIBFLAGS) -o limHaloPT

%.o: %.c
	$(CC) $< -o $@ $(INCLUDES)

clean:
	rm *.o limHaloPT
