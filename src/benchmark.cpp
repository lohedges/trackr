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

#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include "KalmanFilter.h"
#include "TrackGenerator.h"

#include <Eigen/Dense>
#include <Eigen/StdVector>

int main(int argc, char *argv[])
{
    // Set model parameters. This will ultimately be read from file, or passed
    // on the command-line.
    int   num_hits = 5;
    int   num_planes = 5;
    float distance = 1;
    float sigma = 10e-2;
    float theta0 = 1e-3;
    int   seed = 42;
    int   repeats = 100;

    // Get the number of hits from the command-line.
    if (argc > 1)
    {
        std::string arg = argv[1];

        try
        {
            std::size_t pos;
            num_hits = std::stoi(arg, &pos);
            if (pos < arg.size())
            {
                std::cerr << "Trailing characters after number: " << arg << '\n';
            }
        }
        catch (std::invalid_argument const &ex)
        {
            std::cerr << "Invalid number: " << arg << '\n';
        }
        catch (std::out_of_range const &ex)
        {
            std::cerr << "Number out of range: " << arg << '\n';
        }
    }
    // Get the number of repeats from the command-line.
    if (argc > 2)
    {
        std::string arg = argv[2];

        try
        {
            std::size_t pos;
            repeats = std::stoi(arg, &pos);
            if (pos < arg.size())
            {
                std::cerr << "Trailing characters after number: " << arg << '\n';
            }
        }
        catch (std::invalid_argument const &ex)
        {
            std::cerr << "Invalid number: " << arg << '\n';
        }
        catch (std::out_of_range const &ex)
        {
            std::cerr << "Number out of range: " << arg << '\n';
        }
    }

    // Initialise the track generator.
    TrackGenerator trackGenerator(num_planes,
                                  distance,
                                  sigma,
                                  theta0,
                                  seed);

    // Summed throughput.
    double throughput_sum = 0;

    // Initialise a vector to hold the throughput for each run.
    std::vector<double> throughputs(repeats);

    for (int i=0; i<repeats; ++i)
    {
        // Initalise a vector to hold the tracks.
        std::vector<Eigen::MatrixXf> hits(num_hits);

        // Generate the tracks.
        for (int j=0; j<num_hits; ++j)
            hits[j] = std::get<0>(trackGenerator.generateTrack());

        // Initalise the Kalman filter.
        KalmanFilter kalmanFilter(hits, distance, sigma);

        // Record start time.
        auto start = std::chrono::high_resolution_clock::now();

        // Execute the  Kalman filter.
        auto smoothed_hits = kalmanFilter.execute();

        // Record end time.
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;

        // Calculate and store the throughput.
        double secs = std::chrono::duration<double>(elapsed).count();
        double throughput = (num_hits / secs) / 1e6;
        throughput_sum += throughput;
        throughputs[i] = throughput;
    }

    // Calculate the mean throughput.
    double mean = throughput_sum / repeats;

    // Calculate the variance.
    double sum = 0;
    for (int i=0; i<repeats; ++i)
    {
        double tmp = throughputs[i] - mean;
        sum += tmp*tmp;
    }
    double variance = sum / repeats;

    // Output timing statistics.
    std::cout << std::setw(8) << num_hits;
    std::cout << "\t" << std::setw(10) << std::setprecision(6) << std::fixed << mean;
    std::cout << "\t" << std::setw(8) << std::setprecision(6) << std::fixed << variance;
    std::cout << "\t" << std::setw(7) << repeats;
    std::cout << std::endl;

    return 0;
}
