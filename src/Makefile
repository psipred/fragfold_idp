######################################################################
#	Makefile for fragfold_idp components                               #
######################################################################

CC		= cc
LIBS		= -lm

all:	psifilt

install:
	/bin/cp psifilt ../bin

clean:
	/bin/rm -f psifilt
	/bin/rm -f ../bin/psifilt

psifilt:
	$(CC) psifilt.c $(LIBS) -o psifilt
