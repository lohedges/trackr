CXX := g++
CXXFLAGS := -O3 -Isrc -I eigen

# Compile benchmark and test.
all:
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) src/benchmark.cpp src/TrackGenerator.cpp src/KalmanFilter.cpp -o trackr_benchmark
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) src/test.cpp src/TrackGenerator.cpp src/KalmanFilter.cpp -o trackr_test

# Remove binaries.
clean:
	rm trackr_benchmark
	rm trackr_test
