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
    int   num_repeats = 100;
    int   num_batches = 1;

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
                std::cerr << "Trailing characters after number of hits: " << arg << '\n';
                exit(-1);
            }
        }
        catch (std::invalid_argument const &ex)
        {
            std::cerr << "Invalid number of hits: " << arg << '\n';
            exit(-1);
        }
        catch (std::out_of_range const &ex)
        {
            std::cerr << "Number of hits out of range: " << arg << '\n';
            exit(-1);
        }
    }
    // Get the number of repeats from the command-line.
    if (argc > 2)
    {
        std::string arg = argv[2];

        try
        {
            std::size_t pos;
            num_repeats = std::stoi(arg, &pos);
            if (pos < arg.size())
            {
                std::cerr << "Trailing characters after number of repeats: " << arg << '\n';
                exit(-1);
            }
        }
        catch (std::invalid_argument const &ex)
        {
            std::cerr << "Invalid number of repeats: " << arg << '\n';
            exit(-1);
        }
        catch (std::out_of_range const &ex)
        {
            std::cerr << "Number of repeats out of range: " << arg << '\n';
            exit(-1);
        }
    }
    // Get the number of batches from the command-line.
    if (argc > 3)
    {
        std::string arg = argv[3];

        try
        {
            std::size_t pos;
            num_batches = std::stoi(arg, &pos);
            if (pos < arg.size())
            {
                std::cerr << "Trailing characters after number of batches: " << arg << '\n';
                exit(-1);
            }
        }
        catch (std::invalid_argument const &ex)
        {
            std::cerr << "Invalid number of batches: " << arg << '\n';
            exit(-1);
        }
        catch (std::out_of_range const &ex)
        {
            std::cerr << "Number of batches out of range: " << arg << '\n';
            exit(-1);
        }
    }

    // Initalise a track generator for each batch.
    std::vector<TrackGenerator> track_generators(num_batches);
    for (int i=0; i<num_batches; ++i)
    {
        track_generators[i] = TrackGenerator(num_planes,
                                             distance,
                                             sigma,
                                             theta0);
    }

    // Summed throughput.
    double throughput_sum = 0;

    // Initialise a vector to hold the throughput for each repeat.
    std::vector<double> throughputs(num_repeats);

    for (int i=0; i<num_repeats; ++i)
    {
        // Initalise a vector to hold the tracks.
        std::vector<std::vector<Eigen::MatrixXf> > hits(
                num_batches,
                std::vector<Eigen::MatrixXf>(num_hits));

        // Generate the tracks.
        # pragma omp parallel for
        for (int j=0; j<num_batches; ++j)
            for (int k=0; k<num_hits; ++k)
                hits[j][k] = track_generators[j].generateTrack().first;

        // Initalise a Kalman filter for each batch.
        std::vector<KalmanFilter> kalman_filters(num_batches);
        # pragma omp parallel for
        for (int j=0; j<num_batches; ++j)
            kalman_filters[j] = KalmanFilter(hits[j], distance, sigma);

        // Record start time.
        auto start = std::chrono::high_resolution_clock::now();

        // Initalise a vector to hold the smoothed track hits from each batch.
        std::vector<std::vector<Eigen::MatrixXf> > smoothed_hits(
                num_batches,
                std::vector<Eigen::MatrixXf>(num_planes));

        // Execute the Kalman filter for each batch and return the smoothed
        // hits at each plane.
        # pragma omp parallel for
        for (int j=0; j<num_batches; ++j)
            smoothed_hits[j] = kalman_filters[j].execute();

        // Record end time.
        auto finish = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;

        // Calculate and store the throughput.
        double secs = std::chrono::duration<double>(elapsed).count();
        double throughput = ((num_hits * num_batches) / secs) / 1e6;
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
    std::cout << std::setw(8) << num_hits;
    std::cout << "\t" << std::setw(10) << std::setprecision(6) << std::fixed << mean;
    std::cout << "\t" << std::setw(8) << std::setprecision(6) << std::fixed << variance;
    std::cout << "\t" << std::setw(7) << num_repeats;
    std::cout << "\t" << std::setw(7) << num_batches;
    std::cout << std::endl;

    return 0;
}
