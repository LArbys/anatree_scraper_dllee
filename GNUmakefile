
CXX=g++
CXXFLAGS=-g -fPIC -I.

all: scrape_anatree

%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c -o $@ $^

scrape_anatree: scrape_anatree.cxx dwall.o
	$(CXX) $(CXXFLAGS) `root-config --cflags` -c scrape_anatree.cxx -o scrape_anatree.o
	$(CXX) -o scrape_anatree scrape_anatree.o dwall.o `root-config --libs`
	@rm scrape_anatree.o
clean:
	@rm scrape_anatree *.o
