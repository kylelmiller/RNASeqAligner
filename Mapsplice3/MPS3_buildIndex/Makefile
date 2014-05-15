all: buildWholeIndex buildLocalIndex

OPTFLAGS = -O3
OUTPUTDIR = bin/

CFLAGS += $(OPTFLAGS)

buildWholeIndex: directories
	g++ -o $(OUTPUTDIR)buildWholeGenomeIndex buildWholeGenomeIndex_final.cpp

buildLocalIndex: directories
	g++ -o $(OUTPUTDIR)buildLocalIndex buildLocalIndex.cpp

clean:
	rm $(OUTPUTDIR)*

directories:
	mkdir -p $(OUTPUTDIR)