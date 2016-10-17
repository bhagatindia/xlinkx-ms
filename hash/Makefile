CXX = g++
override CXXFLAGS += -g  -std=c++11 -Wall -Wextra -static 

LIBS=lprotobuf

xlinkx-hash: protein_pep_hash.proto xlinkx-hash.cpp
	protoc --cpp_out=./ protein_pep_hash.proto	
	g++ -g -c -std=c++11 xlinkx-hash.cpp protein_pep_hash.pb.cc
	ar rvs xlinkx-hash xlinkx-hash.o protein_pep_hash.pb.o

phd-main:
	g++   -g -std=c++11 phd_main.cpp xlinkx-hash -lprotobuf -o phd-main

clean:
	rm protein_pep_hash.pb.cc protein_pep_hash.pb.h xlinkx-hash.o xlinkx-hash protein_pep_hash.pb.o phd-main

