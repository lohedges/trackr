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

#include <cassert>
#include <cmath>
#include <iostream>

#include "KalmanFilter.h"
#include "TrackGenerator.h"

#include <Eigen/Dense>
#include <Eigen/StdVector>

int main()
{
    // Set model parameters. This will ultimately be read from file, or passed
    // on the command-line.
    int   num_hits = 5;
    int   num_planes = 5;
    float distance = 1;
    float sigma = 10e-2;

    // Initialise a vector to hold the tracks.
    std::vector<Eigen::MatrixXf> hits;

    // Add some reference tracks from the Python implementation.

    Eigen::MatrixXf hit(num_planes, 2);
    hit <<  0.02937927, -0.12062073,
            0.02937927, -0.12062073,
           -0.12062073, -0.42062073,
           -0.12062073, -0.57062073,
           -0.12062073, -0.57062073;
    hits.push_back(hit);

    hit <<  0.17937927, -0.12062073,
            0.47937927, -0.27062073,
            0.62937927, -0.42062073,
            0.77937927, -0.72062073,
            0.92937927, -0.87062073;
    hits.push_back(hit);

    hit <<  0.17937927,  0.17937927,
            0.17937927,  0.32937927,
            0.32937927,  0.62937927,
            0.47937927,  0.62937927,
            0.47937927,  0.77937927;
    hits.push_back(hit);

    hit <<  0.02937927,  0.02937927,
            0.17937927,  0.17937927,
            0.17937927,  0.17937927,
            0.17937927,  0.17937927,
            0.32937927,  0.32937927;
    hits.push_back(hit);

    hit << -0.12062073,   0.17937927,
           -0.27062073,   0.17937927,
           -0.42062073,   0.32937927,
           -0.42062073,   0.47937927,
           -0.57062073,   0.47937927;
    hits.push_back(hit);

    // Initialise the Kalman filter.
    KalmanFilter kalmanFilter(hits, distance, sigma);

    // Execute the Kalman filter and return the smoothed hits at each plane.
    auto smoothed_hits = kalmanFilter.execute();

    // Store some reference eoutput form the Python implementation.

    std::vector<Eigen::MatrixXf> smoothed_hits_ref;

    Eigen::MatrixXf smooth_hit(4, num_hits);
    smooth_hit <<  0.01222589,  0.28513798,  0.19511343,  0.07082304, -0.17921787,
                  -0.03856441,  0.16282753,  0.07284389,  0.05570415, -0.09426857,
                  -0.15350827, -0.17637944,  0.2551298,   0.07082304,  0.19511343,
                  -0.11140831, -0.16282753,  0.13283298,  0.05570415,  0.07284389;
    smoothed_hits_ref.push_back(smooth_hit);

    smooth_hit << -0.02633852,  0.44796551,  0.26795733,  0.1265272,  -0.27348644,
                  -0.03856441,  0.16282753,  0.07284389,  0.05570415, -0.09426857,
                  -0.26491657, -0.33920696,  0.38796278,  0.1265272,   0.26795733,
                  -0.11140831, -0.16282753,  0.13283298,  0.05570415,  0.07284389;
    smoothed_hits_ref.push_back(smooth_hit);

    smooth_hit << -0.06490293,  0.61079304,  0.34080122,  0.18223135, -0.36775501,
                  -0.03856441,  0.16282753,  0.07284389,  0.05570415, -0.09426857,
                  -0.37632488, -0.50203449,  0.52079576,  0.18223135,  0.34080122,
                  -0.11140831, -0.16282753,  0.13283298,  0.05570415,  0.07284389;
    smoothed_hits_ref.push_back(smooth_hit);

    smooth_hit << -0.10346735,  0.77362056,  0.41364511,  0.2379355,  -0.46202358,
                  -0.03856441,  0.16282753,  0.07284389,  0.05570415, -0.09426857,
                  -0.48773319, -0.66486202,  0.65362875,  0.2379355,   0.41364511,
                  -0.11140831, -0.16282753,  0.13283298,  0.05570415,  0.07284389;
    smoothed_hits_ref.push_back(smooth_hit);

    smooth_hit << -0.14203176,  0.93644809,  0.48648901,  0.29363966, -0.55629215,
                  -0.03856441,  0.16282753,  0.07284389,  0.05570415, -0.09426857,
                  -0.59914149, -0.82768954,  0.78646173,  0.29363966,  0.48648901,
                  -0.11140831, -0.16282753,  0.13283298,  0.05570415,  0.07284389;
    smoothed_hits_ref.push_back(smooth_hit);

    // Assert that the smoothed tracks are the same.
    float delta = 1e-6;
    for (unsigned int i=0; i<smoothed_hits.size(); ++i)
    {
        auto a = smoothed_hits[i].reshaped();
        auto b = smoothed_hits_ref[i].reshaped();

        for (int j=0; j<smoothed_hits[i].size(); ++j)
            assert (abs(a[j] - b[j]) < delta);
    }

    std::cout << "Hooray, test passed!\n";

    return 0;
}
