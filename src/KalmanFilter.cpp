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
#include <vector>

#include "KalmanFilter.h"

#ifndef M_PI
    #define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406
#endif

KalmanFilter::KalmanFilter()
{
}

KalmanFilter::KalmanFilter(
        std::vector<Eigen::MatrixXf> hits,
        float distance,
        float sigma) :
        hits(hits)
{
    // Store the number of hits and the number of detector planes.
    this->num_hits = hits.size();
    this->num_planes = hits[0].rows();

    // Convert the hits to a vector containing the hits on each plane.

    // Loop over each of the planes.
    for (int i=0; i<this->num_planes; ++i)
    {
        // Initialise a matrix to hold the hits for this plane.
        Eigen::MatrixXf plane_hits(this->num_hits, 2);

        // Loop over each of the hits.
        for (int j=0; j<this->num_hits; j++)
        {
            plane_hits.row(j) = hits[j].row(i);
        }

        // Append to the vector.
        hits_plane.push_back(plane_hits);
    }

    // Initialise Kalman matrices.

    // Transfer matrix. (Linear operator.)
    this->F << 1, distance, 0,        0,
               0,        1, 0,        0,
               0,        0, 1, distance,
               0,        0, 0,        1;

    // Weight matrix for measurement noise.
    float sigma_sqd = sigma*sigma;
    float x = 1.0 / (sigma_sqd);
    this->G << x, 0, 0, 0,
               0, 0, 0, 0,
               0, 0, x, 0,
               0, 0, 0, 0;

    // Matrix for relationship between measurement 'm' and state 'p'.
    this->H << 1, 0, 0, 0,
               0, 0, 0, 0,
               0, 0, 1, 0,
               0, 0, 0, 0;

    // Initialise the covariance matrix.
    this->C0 << sigma_sqd, 0,    0,         0,
                0,         M_PI, 0,         0,
                0,         0,    sigma_sqd, 0,
                0,         0,    0,         M_PI;
}

std::vector<Eigen::MatrixXf> KalmanFilter::execute()
{
    // Initialise state of system: { x, tan(theta), y, tan(phi) }
    // Use {x, 0, y, 0} for hits from the first plane.
    Eigen::MatrixXf p0 = Eigen::MatrixXf::Zero(4, this->num_hits);
    p0.row(0) = this->hits_plane[0].col(0).transpose();
    p0.row(2) = this->hits_plane[0].col(1).transpose();

    std::vector<Eigen::MatrixXf> p_projs(this->num_planes);
    std::vector<Eigen::MatrixXf> C_projs(this->num_planes);
    std::vector<Eigen::MatrixXf> p_filts(this->num_planes);
    std::vector<Eigen::MatrixXf> C_filts(this->num_planes);

    // Initialise the state and covariance matrix.
    auto p = p0;
    auto C = this->C0;

    // Pre-compute re-used variables.
    auto FT = this->F.transpose();
    auto HTG = this->H.transpose() * this->G;
    auto HTGH = HTG * this->H;

    // Initialise a matrix for the measurements: ({ x, tan(theta), y, tan(phi) }
    Eigen::MatrixXf mm = Eigen::MatrixXf::Zero(4, this->num_hits);

    // Step forward through the planes to linearly project the state and
    // correct (filter) using the measured hits and covariance.
    for (int i=0; i<this->num_planes; ++i)
    {
        // Generate the projections.
        auto p_proj = this->F*p;
        auto C_proj = this->F*C*FT;

        // Store the projected state and covariance.
        p_projs[i] = p_proj;
        C_projs[i] = C_proj;

        // Now perform the filter.

        // Calculate the filtered covariance matrix.
        auto C_filt = (C_proj.inverse() + HTGH).inverse();

        // Extract the hits for the ith plane.
        auto m = this->hits_plane[i].transpose();
        mm.row(0) = m.row(0);
        mm.row(2) = m.row(1);

        // Calculate the filtered state. (N.B. the matrix operations have
        // been reordered to improve performance.)
        auto tmp = (p_proj.transpose() * C_proj.inverse()) + (mm.transpose() * HTG);
        auto p_filt = (tmp * C_filt).transpose();

        // Store the filtered state and covariance.
        p_filts[i] = p_filt;
        C_filts[i] = C_filt;

        // Update the state and covariance.
        p = p_filt;
        C = C_filt;
    }

    // Initialise a vector to hold the smoothed track hits at each plane.
    std::vector<Eigen::MatrixXf> p_smoothed(this->num_planes);

    // Initialise the smmoothed hits for the final plane.
    p_smoothed[this->num_planes-1] = p_filts[this->num_planes-1];

    // Now propagate state and covariance updates backwards through the states
    // to find the globally optimum configuration.
    for (int i=this->num_planes-2; i>=0; --i)
    {
        // Compute the backward transport operator.
        auto A = C_filts[i] * FT * C_projs[i+1].inverse();

        // Smooth in the backwards direction to update the state.
        auto p_smooth = p_filts[i] + A*(p_smoothed[i+1] - p_projs[i+1]);

        // Store the smoothed track.
        p_smoothed[i] = p_smooth;
    }

    // Return the smoothed track hits.
    return p_smoothed;
}
