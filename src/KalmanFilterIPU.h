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

#ifndef _KALMANFILTERIPU_H
#define _KALMANFILTERIPU_H

#include <Eigen/Dense>
#include <Eigen/StdVector>

#include <poplar/DeviceManager.hpp>
#include <poplar/Engine.hpp>
#include <poplar/Graph.hpp>

/// Typedef for poplar compatible RowMajor matrices.
using MatrixRowMajorXf =
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

/// Convenience typedef for DynamicStride.
using DynamicStride = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;

/// Configure how data is mapped from Poplar to Eigen.
using EigenPoplarMapXf = Eigen::Map<MatrixRowMajorXf, Eigen::Unaligned, DynamicStride>;

class KalmanFilterIPU
{
public:
    //! Constructor.
    /*! \param device
            The Poplar device.

        \param hits
            A vector of x,y track hits.

        \param distance
            The distance between detector planes (in metres).

        \param sigma
            The resolution of the detector planes (in metres).
     */
     KalmanFilterIPU(poplar::Device device,
                     std::vector<MatrixRowMajorXf> hits,
                     float distance=1.0,
                     float sigma=10e-2
                    );

    //! Execute the Kalman filter.
    /*! \param secs
            The time taken for the IPU engine to run in seconds.

        \return smoothed_hits
            The reconstructed and smoothed hits at each detector plane.
     */
    std::vector<MatrixRowMajorXf> execute(double&);

private:
    /// The vector of track hits.
    std::vector<MatrixRowMajorXf> hits;

    /// The vector of hits for each detector plane.
    std::vector<MatrixRowMajorXf> hits_plane;

    /// The number of detector planes.
    int num_planes;

    /// The number of track hits.
    int num_hits;

    /// The distance between detector planes.
    float distance;

    /// The resolution of the detector planes.
    float sigma;

    /// Transfer matrix. (Linear operator.)
    MatrixRowMajorXf F;

    /// Weight matrix for measurement noise.
    MatrixRowMajorXf G;

    /// Matrix for relationship between measurement 'm' and state 'p'.
    MatrixRowMajorXf H;

    /// Covariance matrix.
    MatrixRowMajorXf C0;

    /// The Poplar device.
    poplar::Device device;

    /// The Poplar graph.
    poplar::Graph graph;

    /// The Poplar program sequence.
    poplar::program::Sequence prog;

    /// Set up the graph program.
    void setupGraphProgram();
};

#endif	/* _KALMANFILTERIPU_H */
