######################################################################
#	Makefile for fragfold_idp components                               #
######################################################################

CC		= cc
LIBS		= -lm

all:	psifilt rmsdclust

install:
	/bin/cp psifilt ../bin
	/bin/cp rmsdclust_RMSDarray ../bin

clean:
	/bin/rm -f psifilt
	/bin/rm -f ../bin/psifilt
	/bin/rm -f rmsdclust_RMSDarray
	/bin/rm -f ../bin/rmsdclust_RMSDarray

psifilt:
	$(CC) psifilt.c $(LIBS) -o psifilt

rmsdclust:
	$(CC) rmsdclust_RMSDarray.c $(LIBS) -o rmsdclust_RMSDarray
