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
#include <cmath>
#include <iostream>
#include <vector>

#include "KalmanFilterIPU.h"

#ifndef M_PI
    #define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406
#endif

KalmanFilterIPU::KalmanFilterIPU(
        poplar::Device device,
        std::vector<MatrixRowMajorXf> hits,
        float distance,
        float sigma) :
        hits(hits),
        distance(distance),
        sigma(sigma),
        device(std::move(device)),
        graph(this->device)
{
    // Store the number of hits and the number of detector planes.
    this->num_hits = hits.size();
    this->num_planes = hits[0].rows();

    // Convert the hits to a vector containing the hits on each plane.

    // Loop over each of the planes.
    for (int i=0; i<this->num_planes; ++i)
    {
        // Initalise a matrix to hold the hits for this plane.
        MatrixRowMajorXf plane_hits(this->num_hits, 2);

        // Loop over each of the hits.
        for (int j=0; j<this->num_hits; j++)
        {
            plane_hits.row(j) = hits[j].row(i);
        }

        // Append to the vector.
        hits_plane.push_back(plane_hits);
    }

    // Initalise Kalman matrices.

    // Transfer matrix. (Linear operator.)
    this->F = MatrixRowMajorXf(4, 4);
    this->F << 1, distance, 0,        0,
               0,        1, 0,        0,
               0,        0, 1, distance,
               0,        0, 0,        1;

    // Weight matrix for measurement noise.
    float sigma_sqd = sigma*sigma;
    float x = 1.0 / (sigma_sqd);
    this->G = MatrixRowMajorXf(4, 4);
    this->G << x, 0, 0, 0,
               0, 0, 0, 0,
               0, 0, x, 0,
               0, 0, 0, 0;

    // Matrix for relationship between measurement 'm' and state 'p'.
    this->H = MatrixRowMajorXf(4, 4);
    this->H << 1, 0, 0, 0,
               0, 0, 0, 0,
               0, 0, 1, 0,
               0, 0, 0, 0;

    // Initialise the covariance matrix.
    this->C0 = MatrixRowMajorXf(4, 4);
    this->C0 << sigma_sqd, 0,    0,         0,
                0,         M_PI, 0,         0,
                0,         0,    sigma_sqd, 0,
                0,         0,    0,         M_PI;

    // Initalise the graph and apply codelets.
    poplar::Graph graph(this->device);
    graph.addCodelets({"src/KalmanFilterCodelet.cpp"}, "-O3");
    this->graph = std::move(graph);

    // Setup the graph program.
    this->setupGraphProgram();
}

std::vector<MatrixRowMajorXf> KalmanFilterIPU::execute(double &secs)
{
    // Initialise the Poplar engine and load the IPU device.
    poplar::Engine engine(this->graph, this->prog);
    engine.load(this->device);

    // Record start time.
    auto start = std::chrono::high_resolution_clock::now();

    // Run the program.
    engine.run(0);

    // Record end time.
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;

    // Calculate run time in seconds..
    secs = std::chrono::duration<double>(elapsed).count();

    // Read the data back into a single matrix packed as hits along the columns
    // and planes along the rows. (Every fourth row is the start of a new plane.)
    MatrixRowMajorXf p_smooths(4*this->num_planes, this->num_hits);
    engine.readTensor("p_smooths", p_smooths.data());

    // Split the hits by plane. (For consistency with the CPU code.)
    std::vector<MatrixRowMajorXf> p_smoothed;
    for (int i=0; i<this->num_planes; ++i)
    {
        p_smoothed.push_back(p_smooths.block(4*i, 0, 4, this->num_hits));
    }

    return p_smoothed;
}

void KalmanFilterIPU::setupGraphProgram()
{
    // Pre-compute matrices needed for Kalman equations.

    // Create the transpose of the linear operator by mapping to regular
    // column-major Eigen matrix. (Taking the transpose of a row-major
    // matrix doesn't change the layout of data in memory.)
    Eigen::MatrixXf FT = Eigen::Map<
                            Eigen::Matrix<
                                float,
                                Eigen::Dynamic,
                                Eigen::Dynamic,
                                Eigen::RowMajor> >(this->F.data(), 4, 4);
    auto HTG = this->H.transpose() * this->G;
    auto HTGH = HTG * this->H;

    // Initalise state of system: { x, tan(theta), y, tan(phi) }
    // Use {x, 0, y, 0} column vectors for the initial hits. The matrix is
    // packed as hits along the columns and planes along the rows. (Every
    // fourth row is the start of a new plane.)
    MatrixRowMajorXf p0 = MatrixRowMajorXf::Zero(4*this->num_planes, this->num_hits);
    for (int i=0; i<this->num_planes; ++i)
    {
        p0.row(i*4 + 0) = this->hits_plane[i].col(0).transpose();
        p0.row(i*4 + 2) = this->hits_plane[i].col(1).transpose();
    }

    // Number of rows and columns for packed matrices.
    auto nRows = ulong(4*this->num_planes);
    auto nCols = ulong(this->num_hits);

    // Create the graph variables.
    poplar::Tensor v_p = this->graph.addVariable(poplar::FLOAT, {4, nCols}, "p");
    poplar::Tensor v_m = this->graph.addVariable(poplar::FLOAT, {4, nCols}, "m");
    poplar::Tensor v_p0 = this->graph.addVariable(poplar::FLOAT, {nRows, nCols}, "p0");
    poplar::Tensor v_p_proj = this->graph.addVariable(poplar::FLOAT, {4, nCols}, "p_proj");
    poplar::Tensor v_p_projs = this->graph.addVariable(poplar::FLOAT, {nRows, nCols}, "p_projs");
    poplar::Tensor v_p_filt = this->graph.addVariable(poplar::FLOAT, {4, nCols}, "p_filt");
    poplar::Tensor v_p_filts = this->graph.addVariable(poplar::FLOAT, {nRows, nCols}, "p_filts");
    poplar::Tensor v_p_smooth = this->graph.addVariable(poplar::FLOAT, {4, nCols}, "p_smooth");
    poplar::Tensor v_p_smooths = this->graph.addVariable(poplar::FLOAT, {nRows, nCols}, "p_smooths");
    poplar::Tensor v_F = this->graph.addVariable(poplar::FLOAT, {4, 4}, "F");
    poplar::Tensor v_C = this->graph.addVariable(poplar::FLOAT, {4, 4}, "C0");
    poplar::Tensor v_C0 = this->graph.addVariable(poplar::FLOAT, {4, 4}, "C0");
    poplar::Tensor v_C_proj = this->graph.addVariable(poplar::FLOAT, {4, 4}, "C_proj");
    poplar::Tensor v_C_projs = this->graph.addVariable(poplar::FLOAT, {nRows, 4}, "C_projs");
    poplar::Tensor v_C_filt = this->graph.addVariable(poplar::FLOAT, {4, 4}, "C_filt");
    poplar::Tensor v_C_filts = this->graph.addVariable(poplar::FLOAT, {nRows, 4}, "C_filts");
    poplar::Tensor v_FT = this->graph.addVariable(poplar::FLOAT, {4, 4}, "FT");
    poplar::Tensor v_HTG = this->graph.addVariable(poplar::FLOAT, {4, 4}, "HTG");
    poplar::Tensor v_HTGH = this->graph.addVariable(poplar::FLOAT, {4, 4}, "HTGH");
    poplar::Tensor v_tmp_4x4_0 = this->graph.addVariable(poplar::FLOAT, {4, 4}, "tmp_4x4_0");
    poplar::Tensor v_tmp_4x4_1 = this->graph.addVariable(poplar::FLOAT, {4, 4}, "tmp_4x4_1");
    poplar::Tensor v_tmp_4x4_2 = this->graph.addVariable(poplar::FLOAT, {4, 4}, "tmp_4x4_2");
    poplar::Tensor v_tmp_4xN_0 = this->graph.addVariable(poplar::FLOAT, {4, nCols}, "tmp_4xN_0");
    poplar::Tensor v_tmp_4xN_1 = this->graph.addVariable(poplar::FLOAT, {4, nCols}, "tmp_4xN_1");
    this->graph.setTileMapping(v_p, 0);
    this->graph.setTileMapping(v_m, 0);
    this->graph.setTileMapping(v_p0, 0);
    this->graph.setTileMapping(v_p_proj, 0);
    this->graph.setTileMapping(v_p_projs, 0);
    this->graph.setTileMapping(v_p_filt, 0);
    this->graph.setTileMapping(v_p_filts, 0);
    this->graph.setTileMapping(v_p_smooth, 0);
    this->graph.setTileMapping(v_p_smooths, 0);
    this->graph.setTileMapping(v_F, 0);
    this->graph.setTileMapping(v_C, 0);
    this->graph.setTileMapping(v_C0, 0);
    this->graph.setTileMapping(v_C_proj, 0);
    this->graph.setTileMapping(v_C_projs, 0);
    this->graph.setTileMapping(v_C_filt, 0);
    this->graph.setTileMapping(v_C_filts, 0);
    this->graph.setTileMapping(v_FT, 0);
    this->graph.setTileMapping(v_HTG, 0);
    this->graph.setTileMapping(v_HTGH, 0);
    this->graph.setTileMapping(v_tmp_4x4_0, 0);
    this->graph.setTileMapping(v_tmp_4x4_1, 0);
    this->graph.setTileMapping(v_tmp_4x4_2, 0);
    this->graph.setTileMapping(v_tmp_4xN_0, 0);
    this->graph.setTileMapping(v_tmp_4xN_1, 0);

    // Create the graph constants.
    poplar::Tensor c_p0 = this->graph.addConstant<float>(poplar::FLOAT, {nRows, nCols}, p0.data());
    poplar::Tensor c_F = this->graph.addConstant<float>(poplar::FLOAT, {4, 4}, F.data());
    poplar::Tensor c_C0 = this->graph.addConstant<float>(poplar::FLOAT, {4, 4}, C0.data());
    poplar::Tensor c_FT = this->graph.addConstant<float>(poplar::FLOAT, {4, 4}, FT.data());
    poplar::Tensor c_HTG = this->graph.addConstant<float>(poplar::FLOAT, {4, 4}, HTG.eval().data());
    poplar::Tensor c_HTGH = this->graph.addConstant<float>(poplar::FLOAT, {4, 4}, HTGH.eval().data());
    this->graph.setTileMapping(c_p0, 0);
    this->graph.setTileMapping(c_F, 0);
    this->graph.setTileMapping(c_C0, 0);
    this->graph.setTileMapping(c_FT, 0);
    this->graph.setTileMapping(c_HTG, 0);
    this->graph.setTileMapping(c_HTGH, 0);

    // Copy the inital data.
    this->prog.add(poplar::program::Copy(c_p0, v_p0));
    this->prog.add(poplar::program::Copy(c_F, v_F));
    this->prog.add(poplar::program::Copy(c_C0, v_C0));
    this->prog.add(poplar::program::Copy(c_FT, v_FT));
    this->prog.add(poplar::program::Copy(c_HTG, v_HTG));
    this->prog.add(poplar::program::Copy(c_HTGH, v_HTGH));

    // Create a compute set.
    poplar::ComputeSet computeSet = this->graph.addComputeSet("computeSet");

    // Add a vertex to the compute set.
    poplar::VertexRef vtx = this->graph.addVertex(computeSet, "KalmanFilter");

    // Connect graph variables to the vertex.
    this->graph.connect(vtx["p"], v_p);
    this->graph.connect(vtx["m"], v_m);
    this->graph.connect(vtx["p0"], v_p0);
    this->graph.connect(vtx["p_proj"], v_p_proj);
    this->graph.connect(vtx["p_projs"], v_p_projs);
    this->graph.connect(vtx["p_filt"], v_p_filt);
    this->graph.connect(vtx["p_filts"], v_p_filts);
    this->graph.connect(vtx["p_smooth"], v_p_smooth);
    this->graph.connect(vtx["p_smooths"], v_p_smooths);
    this->graph.connect(vtx["F"], v_F);
    this->graph.connect(vtx["FT"], v_FT);
    this->graph.connect(vtx["C"], v_C);
    this->graph.connect(vtx["C0"], v_C0);
    this->graph.connect(vtx["C_proj"], v_C_proj);
    this->graph.connect(vtx["C_projs"], v_C_projs);
    this->graph.connect(vtx["C_filt"], v_C_filt);
    this->graph.connect(vtx["C_filts"], v_C_filts);
    this->graph.connect(vtx["HTG"], v_HTG);
    this->graph.connect(vtx["HTGH"], v_HTGH);
    this->graph.connect(vtx["tmp_4x4_0"], v_tmp_4x4_0);
    this->graph.connect(vtx["tmp_4x4_1"], v_tmp_4x4_1);
    this->graph.connect(vtx["tmp_4x4_2"], v_tmp_4x4_2);
    this->graph.connect(vtx["tmp_4xN_0"], v_tmp_4xN_0);
    this->graph.connect(vtx["tmp_4xN_1"], v_tmp_4xN_1);

    // Map the vertex to a tile.
    this->graph.setTileMapping(vtx, 0);

    // Add the compute set to the program.
    this->prog.add(poplar::program::Execute(computeSet));

    // Create a read to the smoothed tracks.
    this->graph.createHostRead("p_smooths", v_p_smooths);
}
