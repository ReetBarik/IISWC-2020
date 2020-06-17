TARGET = minla
LEDA = ${LEDAROOT}
LIBS = -L${LEDA} -lleda -lX11
#//-lblas -llapack -lm -lg2c -lG -lL 

CXXFLAGS  = -O3 -I${LEDA}/incl -Wno-unused-result 

SRC := $(wildcard *.cpp)
HDR := $(wildcard *.h)
OBJ := $(patsubst %.cpp, %.o, $(SRC))

.phony: ${TARGET}

.o: .cpp
	${CXX} -o $@ ${CXXFLAGS} $<


${TARGET}: ${OBJ} ${HDR}
	$(CXX) -o $@  ${OBJ} ${LIBS}

clean	:
	${RM} ${OBJ} core ${TARGET}



