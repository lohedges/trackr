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

#ifndef _TRACKGENERATOR_H
#define _TRACKGENERATOR_H

#include <random>
#include <utility>

#include <Eigen/Dense>

class TrackGenerator
{
public:
    //! Constructor.
    /*! \param num_planes
            The number of detector planes.

        \param distance
            The distance between detector planes (in metres).

        \param sigma
            The resolution of the detector planes (in metres).

        \param theta0
            Multiple scattering uncertainty (in radians).

        \param seed
            The random number seed. If this is negative then
            a seed will generated using std::random_device.
     */
     TrackGenerator(int   num_planes=5,
                    float distance=1.0,
                    float sigma=10e-2,
                    float theta0=1e-3,
                    int   seed=-1
                   );

    //! Generate a track hit.
    /*! \return num_planes
            The digitised and true track hits at each detector plane.
     */
    std::pair<Eigen::MatrixXf, Eigen::MatrixXf> generateTrack();

private:
    /// The number of detector planes.
    int num_planes;

    /// The distance between detector planes.
    float distance;

    /// The resolution of the detector planes.
    float sigma;

    /// Multiple scattering uncertainty.
    float theta0;

    /// Random number seed.
    int seed;

    // The limit of theta and phi.
    float theta_phi_limit;

    /// Minimum x/y coordinate of detector planes.
    float min_xy;

    /// The Mersenne-Twister generator.
    std::mt19937 generator;

    /// Default uniform_real distribution [-1 - 1].
    std::uniform_real_distribution<float> default_uniform_real_distribution{-1.0, 1.0};

    /// Default normal distribution with zero mean and unit standard deviation.
    std::normal_distribution<float> default_normal_distribution{0.0, 1.0};
};

#endif	/* _TRACKGENERATOR_H */
