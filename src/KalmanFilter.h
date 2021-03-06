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

#ifndef _KALMANFILTER_H
#define _KALMANFILTER_H

#include <Eigen/Dense>
#include <Eigen/StdVector>

class KalmanFilter
{
public:
    //! Constructor.
    /*! \param hits
            A vector of x,y track hits.

        \param distance
            The distance between detector planes (in metres).

        \param sigma
            The resolution of the detector planes (in metres).
     */
     KalmanFilter(std::vector<Eigen::MatrixXf> hits,
                  float distance=1.0,
                  float sigma=10e-2
                );

    //! Execute the Kalman filter.
    /*! \return smoothed_hits
            The reconstructed and smoothed hits at each detector plane.
     */
    std::vector<Eigen::MatrixXf> execute();

private:
    /// The vector of track hits.
    std::vector<Eigen::MatrixXf> hits;

    /// The vector of hits for each detector plane.
    std::vector<Eigen::MatrixXf> hits_plane;

    /// The number of detector planes.
    int num_planes;

    /// The number of track hits.
    int num_hits;

    /// The distance between detector planes.
    float distance;

    /// The resolution of the detector planes.
    float sigma;

    /// Transfer matrix. (Linear operator.)
    Eigen::Matrix4f F;

    /// Weight matrix for measurement noise.
    Eigen::Matrix4f G;

    /// Matrix for relationship between measurement 'm' and state 'p'.
    Eigen::Matrix4f H;

    /// Covariance matrix.
    Eigen::Matrix4f C0;
};

#endif	/* _KALMANFILTER_H */
