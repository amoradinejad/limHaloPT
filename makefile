# NOTE: Depending on your installation paths for gsl and cuba libraries, you should add the links 
# for the INCLUDES and LIBS below, replacing -I/opt/local/include and  -L/opt/local/lib which are the 
# paths on my machine. 

CC        = gcc -c
CLINK     = gcc 
MAKELIB   = ar -rcs

INCLUDES  =  -I./Include/  
INCLUDES +=  -I/opt/local/include
INCLUDES +=  -I./Class/include  

LIBS      =  -L/opt/local/lib  
LIBS     +=  -L./Class/lib
LIBS     +=  -L./lib

CFLAGS    =  -Wall
CFLAGS   += -fopenmp 
CFLAGS   += -g 

LIBFLAGS += -llimHaloPT
LIBFLAGS += -lcuba  
LIBFLAGS += -lclass 
LIBFLAGS += -lgsl 
LIBFLAGS += -lgslcblas 
LIBFLAGS += -lm

SOURCE_DIR = ./Source
BUILD_DIR  = ./Build

MAIN_FILE = $(SOURCE_DIR)/main.c
MAIN_OBJECT = $(BUILD_DIR)/main.o

SOURCE_FILES = $(wildcard $(SOURCE_DIR)/*.c)
$(info $$SOURCE_FILES = $(SOURCE_FILES))

LIB_FILES 	 = $(filter-out $(MAIN_FILE), $(SOURCE_FILES))
$(info $$LIB_FILES = $(LIB_FILES))

LIB_OBJECT_FILES = $(addprefix $(BUILD_DIR)/, $(patsubst %.c, %.o, $(notdir $(LIB_FILES))))
$(info $$LIB_OBJECT_FILES = $(LIB_OBJECT_FILES))

limHaloPT: ${MAIN_OBJECT} lib/liblimHaloPT.a
	$(CLINK) ${MAIN_OBJECT} $(INCLUDES) $(LIBS) $(CFLAGS) $(LIBFLAGS) -o limHaloPT
	rm -rf limHaloPT.dSYM

lib/liblimHaloPT.a: $(LIB_OBJECT_FILES)
	$(MAKELIB) lib/liblimHaloPT.a $(LIB_OBJECT_FILES)

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.c
	$(CC) $< -o $@ $(INCLUDES)

clean:
	rm ${LIB_OBJECT_FILES} $(MAIN_OBJECT) limHaloPT lib/liblimHaloPT.a
