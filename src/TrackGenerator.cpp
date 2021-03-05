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

#include <cmath>
#include <iostream>

#include <Eigen/StdVector>

#include "TrackGenerator.h"

TrackGenerator::TrackGenerator(
        int   num_planes,
        float distance,
        float sigma,
        float theta0,
        int   seed) :
        num_planes(num_planes),
        distance(distance),
        sigma(sigma),
        seed(seed)
{
    // Set the limits of theta and phi. Here they are the same, but they could
    // equally be independent.
    this-> theta_phi_limit = asin(1.0 / 5.0);

    // Work out the minimum x,y coordinate of each detector plane. Again, these
    // could be independent.
    this->min_xy = -(num_planes * distance * tan(this->theta_phi_limit)) - 2.0 * sigma;

    // Seed the random number generator and store the seed.
    if (seed < 0)
        seed = std::random_device{}();
    this->generator.seed(seed);
    this->seed = seed;
}

std::pair<Eigen::MatrixXf, Eigen::MatrixXf> TrackGenerator::generateTrack()
{
    // Initialise matrices to hold the digitised and true track hits.
    Eigen::MatrixXf digit_hits(this->num_planes, 2);
    Eigen::MatrixXf true_hits(this->num_planes, 2);

    // Loop over each detector plane.
    for (int i=0; i<this->num_planes; ++i)
    {
        // Approximate the scattering error for this plane.
        const float scattering_error = (sqrt(i+1) * this->theta0)
                                     * this->default_normal_distribution(this->generator);

        // Generate random theta and phi values.
        const auto theta = this->default_uniform_real_distribution(this->generator)
                         * this->theta_phi_limit;
        const auto phi   = this->default_uniform_real_distribution(this->generator)
                         * this->theta_phi_limit;

        // Work out the true position of the hit on this plane.
        auto x = this->distance * (i+1) * theta;
        auto y = this->distance * (i+1) * phi;

        // Append to the matrix.
        true_hits.row(i) << x, y;

        // Add error to the hit position.
        x += scattering_error;
        y += scattering_error;

        // Now digitize the hit, i.e. work out what detector pixel it corresponds to.

        // Work out the pixel indices.
        const int bin_x = int((x - this->min_xy) / this->sigma);
        const int bin_y = int((y - this->min_xy) / this->sigma);

        // Convert to pixel coordinates.
        x = this->min_xy + (this->sigma * bin_x);
        y = this->min_xy + (this->sigma * bin_y);

        // Append to the matrix.
        digit_hits.row(i) << x, y;
    }

    // Return the track.
    return std::pair<Eigen::MatrixXf, Eigen::MatrixXf>{digit_hits, true_hits};
}
