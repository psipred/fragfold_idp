######################################################################
#	Makefile for PSIpred4                                        #
######################################################################

#### CC		Command to execute C compiler
#### CFLAGS	Flags to pass to C compiler.

CC		= cc
LIBS		= -lm

all:	psifilt

install:
	/bin/cp psifilt ../bin

clean:
	/bin/rm -f psifilt
	/bin/rm -f ../bin/psifilt


chkparse:
	$(CC) psifilt.c $(LIBS) -o psifilt
