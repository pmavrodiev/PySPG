# put common definitions in here

TARGET          = rank
TARGET_SCRIPT   = script-rank
CONFIGFILE      = 

BASEDIR         = $(HOME)/opt
CPP             = g++
CPPFLAGS        = -O2
FORTRAN         = gfortran
CC              = gcc
LD              = g++
CTT             = ~/opt/bin/ctt.py
LDFLAGS         = -L${BASEDIR}/lib -lm -ldranxor -lgsl -lgslcblas
# -lwrappers  -lgauges -lgenerator -lnetworks -lstochastic -lutil -linteraction -lmeasures -lmcsys -ldynsys
RM              = rm
SHELL           = /bin/bash
INCLUDES        = -I${BASEDIR}/include
EXEEXT          = .p4
EXEEXT_SCRIPT   = .sh
EXEPREFIX       = ctx-
OBJS            = base.o main.o 

.PHONY: clean all install
.PRECIOUS: base.cxx #Pavlin: want to keep it

all:  ${EXEPREFIX}${TARGET}${EXEEXT}

install:
	install ${EXEPREFIX}${TARGET_SCRIPT}${EXEEXT_SCRIPT} ${BASEDIR}/bin
	install ${EXEPREFIX}${TARGET}${EXEEXT} ${BASEDIR}/bin
	install -T -m 0644 conf.in ${BASEDIR}/etc/spg-conf/${TARGET_SCRIPT}.in
	install -T -m 0644 conf.stdout ${BASEDIR}/etc/spg-conf/${TARGET_SCRIPT}.stdout
	install -T -m 0644 conf.params ${BASEDIR}/etc/spg-conf/${TARGET_SCRIPT}.params
	install -T -m 0644 base.ct ${BASEDIR}/etc/spg-conf/${TARGET_SCRIPT}.ct


${EXEPREFIX}${TARGET}${EXEEXT}: ${OBJS}
	${LD} ${OBJS} -o $@ ${CPPFLAGS} ${LDFLAGS}

%.o : %.cxx
	${CPP} ${INCLUDES} ${CPPFLAGS} -c $<

%.cxx : %.ct
	${CTT} -n -i $<

clean:
	${RM} -f *.o ${EXEPREFIX}*${EXEEXT}

