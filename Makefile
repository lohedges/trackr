CXX := g++
CXXFLAGS := -O3 -Isrc -I eigen
ifeq ($(CXX),g++)
	OMPFLAGS := -fopenmp
else ifeq ($(CXX),clang++)
	OMPFLAGS += -fopenmp=libomp
endif

# Compile benchmark and test. We could create and link against a re-usable
# library, but will leave for now.
all:
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) src/test.cpp src/TrackGenerator.cpp src/KalmanFilter.cpp -o trackr_test
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) src/benchmark.cpp src/TrackGenerator.cpp src/KalmanFilter.cpp -o trackr_benchmark
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(OPTFLAGS) src/benchmark_omp.cpp src/TrackGenerator.cpp src/KalmanFilter.cpp -o trackr_benchmark_omp

# Remove binaries.
clean:
	rm trackr_test
	rm trackr_benchmark
	rm trackr_benchmark_omp
