# Makefile for halogen

# Executable

BASE    = halogen
VERSION = 1.1
EXT     = 64
EXE     = $(BASE)$(EXT)

# Compiler stuff

CC	= gcc -DDPP
CFLAGS	= -O3 -Wall
LIBS	= -lm

# Object definition

OBJ	= $(BASE).o functions.o routines.o IOfunctions.o

# Rules

$(EXE):	$(OBJ) Makefile
	$(CC) $(CFLAGS) $(OBJ) -o $(EXE) $(LIBS)

clean:
	-rm -f *.o *~ $(EXE)

tar:
	cd ..; tar cvf - $(BASE)/*.c $(BASE)/*.h $(BASE)/Makefile > $(BASE)-$(VERSION).tar

# Dependencies

halogen.o: definitions.h functions.h routines.h IOfunctions.h
functions.o: definitions.h functions.h routines.h IOfunctions.h
routines.o: definitions.h functions.h routines.h IOfunctions.h

