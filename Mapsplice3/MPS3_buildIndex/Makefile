all: buildWholeIndex buildAll2ndLevelIndex

OPTFLAGS = -O3
OUTPUTDIR = bin/

CFLAGS += $(OPTFLAGS)

buildWholeIndex: directories
	g++ -o $(OUTPUTDIR)buildWholeGenomeIndex buildWholeGenomeIndex_final.cpp

buildAll2ndLevelIndex: directories
	g++ -o $(OUTPUTDIR)buildAll2ndLevelIndex buildAll2ndLevelIndex.cpp

clean:
	rm $(OUTPUTDIR)*

directories:
	mkdir -p $(OUTPUTDIR)