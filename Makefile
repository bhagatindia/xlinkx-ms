CXX = g++
MSTOOLKIT = mstoolkit
HASH = hash
PROTOBUF = protobuf
HARDKLOR = hardklor
override CXXFLAGS +=  -g  -std=c++11 -Wall -Wextra -static -Wno-char-subscripts -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D__LINUX__ -I$(MSTOOLKIT)/include
EXECNAME = xlinkx.exe
OBJS = xlinkx.o xlinkx_Preprocess.o xlinkx_Search.o xlinkx_MassSpecUtils.o  $(HASH)/xlinkx-hash.o $(HASH)/protein_pep_hash.pb.o
DEPS = xlinkx.h Common.h xlinkx_Data.h xlinkx_DataInternal.h xlinkx_Preprocess.h xlinkx_MassSpecUtils.h

LIBS = -L$(MSTOOLKIT) -lmstoolkitlite -lm -pthread -L/usr/local/lib -lprotobuf 
ifdef MSYSTEM
   LIBS += -lws2_32
endif



xlinkx.exe: $(OBJS)
	git submodule init; git submodule update
	cd $(MSTOOLKIT) ; make lite 
	cd $(HASH) ; make; 
	rm MSToolkit; ln -s mstoolkit MSToolkit
	cd $(HARDKLOR) ; make; 
	${CXX} $(CXXFLAGS) $(OBJS) $(LIBS) -o ${EXECNAME}

xlinkx.o: xlinkx.cpp $(DEPS)
	cd $(HASH) ; make; 
	git submodule init; git submodule update
	${CXX} ${CXXFLAGS} xlinkx.cpp -c

xlinkx_Preprocess.o: xlinkx_Preprocess.cpp Common.h xlinkx_Preprocess.h xlinkx.h Common.h xlinkx_Data.h xlinkx_DataInternal.h
	git submodule init; git submodule update
	${CXX} ${CXXFLAGS} xlinkx_Preprocess.cpp -c

xlinkx_Search.o: xlinkx_Search.cpp Common.h xlinkx_Search.h xlinkx.h Common.h xlinkx_Data.h xlinkx_DataInternal.h
	git submodule init; git submodule update
	${CXX} ${CXXFLAGS} xlinkx_Search.cpp -c

xlinkx_MassSpecUtils.o: xlinkx_MassSpecUtils.cpp Common.h xlinkx_MassSpecUtils.h xlinkx.h Common.h xlinkx_Data.h xlinkx_DataInternal.h
	git submodule init; git submodule update
	${CXX} ${CXXFLAGS} xlinkx_MassSpecUtils.cpp -c

clean:
	rm -f *.o ${EXECNAME}
	cd $(MSTOOLKIT) ; make clean
	cd $(HASH) ; make clean
	cd $(PROTOBUF); make clean
