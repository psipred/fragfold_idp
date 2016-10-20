######################################################################
#	Makefile for fragfold_idp components                               #
######################################################################

CC		= cc
LIBS		= -lm

all:	psifilt rmsdclust

install:
	/bin/cp psifilt ../bin
	/bin/cp rmsdclust ../bin

clean:
	/bin/rm -f psifilt
	/bin/rm -f ../bin/psifilt
	/bin/rm -f rmsdclust
	/bin/rm -f ../bin/rmsdclust

psifilt:
	$(CC) psifilt.c $(LIBS) -o psifilt

rmsdclust:
	$(CC) rmsdclust.c $(LIBS) -o rmsdclust
