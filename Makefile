## compiler & flags
CXX = g++
DEBUG = -g
CXXFLAGS = -Wall -c -O3 $(DEBUG)
INC = -I/usr/local/include/htslib
LIBS = -Wall -lhts 

OBJECTS = fasta.o sam_util.o main.o
HEADERS = fasta.h sam_util.h


read_sam : $(OBJECTS)
	$(CXX) $(LIBS) $(OBJECTS) -o read_sam

main.o : main.cpp fasta.h sam_util.h
	$(CXX) $(CXXFLAGS) $(INC) main.cpp

sam_util.o : sam_util.cpp sam_util.h
	$(CXX) $(CXXFLAGS) $(INC) sam_util.cpp

fasta.o : fasta.cpp fasta.h
	$(CXX) $(CXXFLAGS) fasta.cpp


clean : 
	\rm *.o *~ read_sam
