# Makefile for halogen

# Executable

BASE    = halogen
VERSION = 1.1
EXT     = 
EXE     = $(BASE)$(EXT)

# Compiler stuff

CC	= gcc
CFLAGS	= -O3 -Wall -I$(LOCAL_LIB_PATH)/include
LIBS	= -L$(LOCAL_LIB_PATH)/lib -lm -liof

# Object definition

OBJ	= $(BASE).o functions.o routines.o \
	arguments.o allocate.o check.o initialise.o usage.o write.o
INCL	= definitions.h functions.h routines.h \
	arguments.h allocate.h check.h initialise.h usage.h write.h

# Rules

$(EXE):	$(OBJ) Makefile
	$(CC) $(CFLAGS) $(OBJ) -o $(EXE) $(LIBS)

clean:
	-rm -f *.o *~ $(EXE)

tar:
	cd ..; tar cvf - $(BASE)/*.c $(BASE)/*.h $(BASE)/Makefile > $(BASE)-$(VERSION).tar

# Dependencies

halogen.o: $(INCL)
functions.o: $(INCL)
routines.o: $(INCL)
arguments.o: $(INCL)
allocate.o: $(INCL)
check.o: $(INCL)
initialise.o: $(INCL)
usage.o: $(INCL)
write.o: $(INCL)


