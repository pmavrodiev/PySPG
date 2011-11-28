# put common definitions in here

TARGET          = random-walk

CONFIGFILE      = 

BASEDIR         = $(HOME)/opt
CPP             = g++
CPPFLAGS        = -O2
FORTRAN         = gfortran
CC              = gcc
LD              = g++
CTT             = ~/opt/bin/ctt.py
LDFLAGS         = -L${BASEDIR}/lib -lm -ldranxor
<<<<<<< HEAD
# -lwrappers  -lgauges -lgenerator -lnetworks -lstochastic -lutil -linteraction -lmeasures -lmcsys -ldynsys
=======
>>>>>>> 7cf2fdd8ed287c32cb2b824f5cf1ccdde2782367
RM              = rm
SHELL           = /bin/bash
INCLUDES        = -I${BASEDIR}/include
EXEEXT          = .p4
EXEPREFIX       = ctx-
OBJS            = base.o main.o 

.PHONY: clean all install
<<<<<<< HEAD
.PRECIOUS: base.cxx #Pavlin: want to keep it
=======
.PRECIOUS: base.cxx


>>>>>>> 7cf2fdd8ed287c32cb2b824f5cf1ccdde2782367

all:  ${EXEPREFIX}${TARGET}${EXEEXT}

install:
	install ${EXEPREFIX}${TARGET}${EXEEXT} ${BASEDIR}/bin
	
	mkdir -p ${BASEDIR}/etc/spg-conf
	install -T -m 0644 conf.in ${BASEDIR}/etc/spg-conf/${TARGET}.in
	install -T -m 0644 conf.stdout ${BASEDIR}/etc/spg-conf/${TARGET}.stdout
	install -T -m 0644 conf.params ${BASEDIR}/etc/spg-conf/${TARGET}.params
	install -T -m 0644 base.ct ${BASEDIR}/etc/spg-conf/${TARGET}.ct


${EXEPREFIX}${TARGET}${EXEEXT}: ${OBJS}
	${LD} ${OBJS} -o $@ ${CPPFLAGS} ${LDFLAGS}

%.o : %.cxx
	${CPP} ${INCLUDES} ${CPPFLAGS} -c $<

%.cxx : %.ct
	${CTT} -n -i $<

clean:
	${RM} -f *.o ${EXEPREFIX}*${EXEEXT}

