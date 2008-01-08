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

OBJ	= $(BASE).o arguments.o functions.o routines.o write.o check.o usage.o
INCL	= definitions.h arguments.h functions.h routines.h write.h check.h usage.h

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
write.o: $(INCL)
check.o: $(INCL)
usage.o: $(INCL)

