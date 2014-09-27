all: buildWholeIndex buildAll2ndLevelIndex

OPTFLAGS =  -g -O0 -rdynamic -Wall -std=c++0x
OUTPUTDIR = bin/

CFLAGS += $(OPTFLAGS)

buildWholeIndex: directories
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildWholeGenomeIndex buildWholeGenomeIndex_final.cpp

buildAll2ndLevelIndex: directories
	g++ $(CFLAGS) -o $(OUTPUTDIR)buildAll2ndLevelIndex buildAll2ndLevelIndex.cpp

clean:
	rm -f $(OUTPUTDIR)buildWholeGenomeIndex
	rm -f $(OUTPUTDIR)buildAll2ndLevelIndex

directories:
	mkdir -p $(OUTPUTDIR)