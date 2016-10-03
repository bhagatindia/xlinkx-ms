CXX = g++
override CXXFLAGS += -g  -std=c++11 -Wall -Wextra -static 

LIBS=lprotobuf

xlinkx-hash: protein_pep_hash.proto xlinkx-hash.cpp
	protoc --cpp_out=./ protein_pep_hash.proto	
	g++ -g -std=c++11 -lprotobuf xlinkx-hash.cpp protein_pep_hash.pb.cc  -o xlinkx-hash

clean:
	rm protein_pep_hash.pb.cc protein_pep_hash.pb.h

