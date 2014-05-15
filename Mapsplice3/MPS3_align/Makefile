all: alignAll alignHeadTail alignMid alignOneEndUnmapped compareJunction

OPTFLAGS = -O3
OUTPUTDIR = bin/

CFLAGS += $(OPTFLAGS) -fopenmp

alignAll: directories
	g++  $(CFLAGS) -o $(OUTPUTDIR)alignAll alignAll.cpp

alignHeadTail: directories
	g++  $(CFLAGS) -o $(OUTPUTDIR)alignHeadTail alignHeadTail.cpp

alignMid: directories
	g++  $(CFLAGS) -o $(OUTPUTDIR)alignMid alignMid_main.cpp

alignOneEndUnmapped: directories
	g++  $(CFLAGS) -o $(OUTPUTDIR)alignOneEndUnmapped alignOneEndUnmapped.cpp

compareJunction: directories
	g++  $(CFLAGS) -o $(OUTPUTDIR)compareJunction compareJunction.cpp
	
clean:
	rm $(OUTPUTDIR)*

directories:
	mkdir -p $(OUTPUTDIR)