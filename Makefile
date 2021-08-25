.DEFAULT_GOAL := all

CXX := g++
CXXFLAGS := -O3 -Isrc -I eigen
ifeq ($(CXX),g++)
	OMPFLAGS := -fopenmp
else ifeq ($(CXX),clang++)
	OMPFLAGS += -fopenmp=libomp
endif

ABIFLAGS := -D_GLIBCXX_USE_CXX11_ABI=0
POPLARFLAGS :=-std=c++17 -L/opt/poplar/lib64 -lpoplar

# Compile benchmark and test. We could create and link against a re-usable
# library, but will leave for now.
all:
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) src/test.cpp src/TrackGenerator.cpp src/KalmanFilter.cpp -o trackr_test
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) src/benchmark.cpp src/TrackGenerator.cpp src/KalmanFilter.cpp -o trackr_benchmark
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(OPTFLAGS) src/benchmark_omp.cpp src/TrackGenerator.cpp src/KalmanFilter.cpp -o trackr_benchmark_omp

# Compile benchmark and test for IPU.
ipu:
	$(CXX) $(CXXFLAGS) $(ABIFLAGS) $(POPLARFLAGS) $(OPTFLAGS) src/benchmark_ipu.cpp src/TrackGenerator.cpp src/KalmanFilterIPU.cpp -o trackr_benchmark_ipu
	$(CXX) $(CXXFLAGS) $(ABIFLAGS) $(POPLARFLAGS) $(OPTFLAGS) src/test_ipu.cpp src/TrackGenerator.cpp src/KalmanFilterIPU.cpp -o trackr_test_ipu

# Remove binaries.
clean:
	rm -f trackr_test
	rm -f trackr_benchmark
	rm -f trackr_benchmark_omp
	rm -f trackr_test_ipu
	rm -f trackr_benchmark_ipu
