# Makefile for halogen

NAME	= halogen
VERSION	= $(shell git describe --tags --long)

CC		= gcc
CFLAGS	= -O3 -mcmodel=medium -Wall -pedantic -I$(LOCAL_LIB_PATH)/include -DVERSION=\"${VERSION}\"
LIBS	= -L$(LOCAL_LIB_PATH)/lib -lm -liof

SRCS	= $(wildcard *.c)
INCS	= $(wildcard *.h)
OBJS	= $(SRCS:.c=.o)

# Rules

$(NAME): $(OBJS) $(INCS) Makefile
	$(CC) $(CFLAGS) $(OBJS) -o $(NAME) $(LIBS)

clean:
	rm -f *~ *.o $(NAME)
