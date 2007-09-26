# Makefile for halogen

# Executable

BASE    = halogen
EXT     = 64
EXE     = $(BASE)$(EXT)

# Compiler stuff

CC	= gcc
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
	cd ..; tar cvf - $(BASE)/*.c $(BASE)/*.h $(BASE)/Makefile > $(BASE).tar

# Dependencies

halogen.o: definitions.h functions.h routines.h IOfunctions.h
functions.o: definitions.h functions.h routines.h IOfunctions.h
routines.o: definitions.h functions.h routines.h IOfunctions.h

