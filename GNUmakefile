
CXX=g++
CXXFLAGS=-g -fPIC -I.

all: scrape_anatree

scrape_anatree: scrape_anatree.cxx
	$(CXX) $(CXXFLAGS) `root-config --cflags` -c scrape_anatree.cxx -o scrape_anatree.o
	$(CXX) -o scrape_anatree scrape_anatree.o `root-config --libs`
	@rm scrape_anatree.o
clean:
	@rm scrape_anatree
