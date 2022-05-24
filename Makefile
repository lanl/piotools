CC  = gcc-10
CXX = g++-10
FC  = gfortran-10

OPT ?= -O3
CFLAGS   = ${OPT}
CXXFLAGS = ${OPT}
FTNFLAGS = ${OPT}

.cpp.o: $< 
	$(CXX) -c  $(CXXFLAGS) $(INCLUDE)  $<

.F90.o: $< 
	echo	$(FC) -c  $(FTNFLAGS) $(INCLUDE)  $<


SRCS = pioTest.F90 pio_interface.F90 pioInterface.cpp
OBJS2 = $(OBJS1:.cpp=.o)
OBJS1 = $(SRCS:.F90=.o)
OBJS = $(OBJS2:.c=.o)

TFILE ?= triplepoint-dmp000119
all: pioctest pioftest

test: all
	pioctest ${TFILE}
	pioftest ${TFILE}

pioTest.o: pioTest.F90 pio_interface.o pioInterface.o
	$(FC) -c   $(FTNFLAGS) $(INCLUDE) $<

pio_interface.o: pio_interface.F90 
	$(FC) -c   $(FTNFLAGS) $(INCLUDE)  $<

pioInterface.o: pioInterface.cpp
	$(CXX) -c   $(CXXFLAGS) $(INCLUDE)  $<

pioftest: $(OBJS)
	$(FC) $(OBJS) -lstdc++ -o pioftest -flto

pioctest: pioInterface.cpp
	${CXX} -o pioctest ${CXXFLAGS} $< ${CXXFLAGS} -DDOPIOMAIN -flto

clean:
	/bin/rm *.o *.mod pioctest pioftest
