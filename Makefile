CXX = g++
MSTOOLKIT = mstoolkit
override CXXFLAGS += -O3 -Wall -Wextra -static -Wno-char-subscripts -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D__LINUX__ -I$(MSTOOLKIT)/include
EXECNAME = xlinkx.exe
OBJS = xlinkx.o
DEPS = 

LIBPATHS = -L$(MSTOOLKIT)
LIBS = -lmstoolkitlite -lm -lpthread
ifdef MSYSTEM
   LIBS += -lws2_32
endif

xlinkx.exe: $(OBJS)
	git submodule init; git submodule update
	cd $(MSTOOLKIT) ; make lite 
	${CXX} $(CXXFLAGS) $(OBJS) $(LIBPATHS) $(LIBS) -o ${EXECNAME}

xlinkx.o: xlinkx.cpp $(DEPS)
	git submodule init; git submodule update
	${CXX} ${CXXFLAGS} xlinkx.cpp -c

clean:
	rm -f *.o ${EXECNAME}
	cd $(MSTOOLKIT) ; make clean
