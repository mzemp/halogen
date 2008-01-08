# Makefile for halogen

# Executable

BASE    = halogen
VERSION = 1.1
EXT     = 64
EXE     = $(BASE)$(EXT)

# Compiler stuff

CC	= gcc
CFLAGS	= -O3 -Wall -I$(LOCAL_LIB_PATH)/include
LIBS	= -L$(LOCAL_LIB_PATH)/lib -lm -liof

# Object definition

OBJ	= $(BASE).o functions.o routines.o

# Rules

$(EXE):	$(OBJ) Makefile
	$(CC) $(CFLAGS) $(OBJ) -o $(EXE) $(LIBS)

clean:
	-rm -f *.o *~ $(EXE)

tar:
	cd ..; tar cvf - $(BASE)/*.c $(BASE)/*.h $(BASE)/Makefile > $(BASE)-$(VERSION).tar

# Dependencies

halogen.o: definitions.h functions.h routines.h
functions.o: definitions.h functions.h routines.h
routines.o: definitions.h functions.h routines.h

