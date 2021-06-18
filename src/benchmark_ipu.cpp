/*
  Copyright (c) 2021 Lester Hedges <lester.hedges@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>

#include "KalmanFilterIPU.h"
#include "TrackGenerator.h"

#include <Eigen/Dense>
#include <Eigen/StdVector>

#include <poplar/DeviceManager.hpp>

/// Setup and connect to a Poplar IPU device.
poplar::Device setIpuDevice();

int main(int argc, char *argv[])
{
    // Set model parameters. This will ultimately be read from file, or passed
    // on the command-line.
    int   num_tiles = 1;
    int   hits_per_tile = 5;
    int   num_planes = 5;
    float distance = 1;
    float sigma = 10e-2;
    float theta0 = 1e-3;
    int   seed = 42;
    int   num_repeats = 100;
    int   num_hits;
    bool  profile = false;

    // Get the number of tiles from the command-line.
    if (argc > 1)
    {
        std::string arg = argv[1];

        try
        {
            std::size_t pos;
            num_tiles = std::stoi(arg, &pos);
            if (pos < arg.size())
            {
                std::cerr << "Trailing characters after number: " << arg << '\n';
                exit(-1);
            }
        }
        catch (std::invalid_argument const &ex)
        {
            std::cerr << "Invalid number: " << arg << '\n';
            exit(-1);
        }
        catch (std::out_of_range const &ex)
        {
            std::cerr << "Number out of range: " << arg << '\n';
            exit(-1);
        }
    }
    // Get the number of hits-per-tile from the command-line.
    if (argc > 2)
    {
        std::string arg = argv[2];

        try
        {
            std::size_t pos;
            hits_per_tile = std::stoi(arg, &pos);
            if (pos < arg.size())
            {
                std::cerr << "Trailing characters after number: " << arg << '\n';
                exit(-1);
            }
        }
        catch (std::invalid_argument const &ex)
        {
            std::cerr << "Invalid number: " << arg << '\n';
            exit(-1);
        }
        catch (std::out_of_range const &ex)
        {
            std::cerr << "Number out of range: " << arg << '\n';
            exit(-1);
        }
    }
    // Get the number of repeats from the command-line.
    if (argc > 3)
    {
        std::string arg = argv[3];

        try
        {
            std::size_t pos;
            num_repeats = std::stoi(arg, &pos);
            if (pos < arg.size())
            {
                std::cerr << "Trailing characters after number: " << arg << '\n';
                exit(-1);
            }
        }
        catch (std::invalid_argument const &ex)
        {
            std::cerr << "Invalid number: " << arg << '\n';
            exit(-1);
        }
        catch (std::out_of_range const &ex)
        {
            std::cerr << "Number out of range: " << arg << '\n';
            exit(-1);
        }
    }
    // Get the profiling flag from the command-line.
    if (argc > 4)
    {
        std::stringstream ss(argv[4]);

        if (!(ss >> std::boolalpha >> profile))
        {
            std::cerr << "Invalid profiling flag: options are 'true' or 'false'\n";
            exit(-1);
        }
    }

    // Work out the number of hits.
    num_hits = num_tiles * hits_per_tile;

    // Initialise the track generator.
    TrackGenerator trackGenerator(num_planes,
                                  distance,
                                  sigma,
                                  theta0,
                                  seed);

    // Summed throughput.
    double throughput_sum = 0;

    // Initialise a vector to hold the throughput for each repeat.
    std::vector<double> throughputs(num_repeats);

    for (int i=0; i<num_repeats; ++i)
    {
        // Initalise a vector to hold the tracks.
        std::vector<MatrixRowMajorXf> hits(num_hits);

        // Generate the tracks.
        for (int j=0; j<num_hits; ++j)
            hits[j] = trackGenerator.generateTrack().first;

        // Setup the IPU device handle.
        auto device = setIpuDevice();

        // Initalise the Kalman filter.
        KalmanFilterIPU kalmanFilter(std::move(device),
                                     hits,
                                     num_tiles,
                                     hits_per_tile,
                                     distance,
                                     sigma);

        // The time for the IPU engine to run in seconds.
        double secs;

        // Execute the Kalman filter and return the smoothed hits at each plane,
        // performing a "warmup" run prior to timing the benchmark.
        auto smoothed_hits = kalmanFilter.execute(secs, true, profile);

        double throughput = (num_hits / secs) / 1e6;
        throughput_sum += throughput;
        throughputs[i] = throughput;
    }

    // Calculate the mean throughput.
    double mean = throughput_sum / num_repeats;

    // Calculate the variance.
    double sum = 0;
    for (int i=0; i<num_repeats; ++i)
    {
        double tmp = throughputs[i] - mean;
        sum += tmp*tmp;
    }
    double variance = sum / num_repeats;

    // Output timing statistics.
    std::cout << std::setw(8) << num_tiles;
    std::cout << "\t" << std::setw(17) << hits_per_tile;
    std::cout << "\t" << std::setw(10) << std::setprecision(6) << std::fixed << mean;
    std::cout << "\t" << std::setw(8) << std::setprecision(6) << std::fixed << variance;
    std::cout << "\t" << std::setw(7) << num_repeats;
    std::cout << std::endl;

    return 0;
}

poplar::Device setIpuDevice()
{
    auto dm = poplar::DeviceManager::createDeviceManager();
    auto hwDevices = dm.getDevices(poplar::TargetType::IPU, 1);
    if (hwDevices.size() > 0)
    {
        for (auto &d : hwDevices)
        {
            if (d.attach())
            {
                return std::move(d);
            }
        }
    }

    std::cerr << "Unable to connect to IPU device!\n";
    exit(-1);
}
