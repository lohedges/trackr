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
#include <fstream>
#include <vector>

#include "KalmanFilterIPU.h"

#ifndef M_PI
    #define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406
#endif

KalmanFilterIPU::KalmanFilterIPU(
        poplar::Device device,
        std::vector<MatrixRowMajorXf> hits,
        int   num_tiles,
        int   hits_per_tile,
        float distance,
        float sigma,
        bool is_test) :
        hits(hits),
        num_tiles(num_tiles),
        hits_per_tile(hits_per_tile),
        device(std::move(device)),
        graph(this->device)
{
    // Work out the number of hardware threads per-tile.
    auto num_workers = this->device.getTarget().getNumWorkerContexts();

    if (not is_test)
    {
        // We need at least (4 x num_workers) hits per tile.
        if (hits_per_tile < 4*num_workers)
        {
            std::cerr << "The minimum number of hits per tile is 24!\n";
            exit(-1);
        }

        // Make sure the number of hits is a multiple of num_workers.
        if (hits_per_tile % num_workers != 0)
        {
            std::cerr << "The number of hits per tile must be divisible by "
                      << num_workers << "\n";
            exit(-1);
        }

        // Make sure the number of hits per worker is a multiple of 4. (For
        // certain operations we process 4 floats simulatneously by casting
        // to float4.) If loop unrolling, we will also need to ensure that
        // the hits per worker is a multiple of 4 times the unroll factor.
        auto hits_per_worker = hits_per_tile / num_workers;
        if (hits_per_worker % 4 != 0)
        {
            std::cerr << "There are " << num_workers << " workers per tile. "
                      << "The number of hits per worker must be equal to 4. "
                      << "Currently equal to " << hits_per_worker << ".\n";
            exit(-1);
        }
    }

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

    // Initalise the graph and apply codelets. Set optimisation level to -O3,
    // which gives the best performance at present. (The default when no
    // additional compiler flags are specified is -O2.)
    poplar::Graph graph(this->device);
    graph.addCodelets({"src/KalmanFilterCodelet.cpp"}, "-O3");
    this->graph = std::move(graph);

    // Setup the graph program.
    this->setupGraphProgram();
}

std::vector<MatrixRowMajorXf> KalmanFilterIPU::execute(double &secs,
                                                       bool warmup,
                                                       bool profile)
{
    // Initialise the Poplar engine and load the IPU device.

    auto optionFlags = poplar::OptionFlags{};
    if (profile)
    {
        // Taken from UoB-HPC IPU cookbook:
        // https://github.com/UoB-HPC/ipu-hpc-cookbook
        optionFlags = poplar::OptionFlags
        {
            {"target.saveArchive",                "archive.a"},
            {"debug.instrument",                  "true"},
            {"debug.instrumentCompute",           "true"},
            {"debug.loweredVarDumpFile",          "vars.capnp"},
            {"debug.instrumentControlFlow",       "true"},
            {"debug.computeInstrumentationLevel", "tile"},
            {"debug.outputAllSymbols",            "true"},
            {"autoReport.all",                    "true"},
            {"autoReport.outputSerializedGraph",  "true"},
            {"debug.retainDebugInformation",      "true"}
        };
    }

    poplar::Engine engine(this->graph, this->prog, optionFlags);
    engine.load(this->device);

    // Perform a warmup run.
    if (warmup)
        engine.run(0);

    // Record start time.
    auto start = std::chrono::high_resolution_clock::now();

    // Run the program.
    engine.run(0);

    // Record end time.
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;

    // Calculate run time in seconds..
    secs = std::chrono::duration<double>(elapsed).count();

    // Write profiling information to file.
    if (profile)
    {
        std::ofstream profile;
        profile.open("profile.txt");
        engine.printProfileSummary(profile, {{"showExecutionSteps", "true"}});
        profile.close();
    }

    std::vector<MatrixRowMajorXf> p_smoothed;
    if (this->num_tiles == 1)
    {
        // Read the data back into a single matrix packed as hits along the columns
        // and planes along the rows. (Every fourth row is the start of a new plane.)
        MatrixRowMajorXf p_smooths(4*this->num_planes, this->hits_per_tile);
        engine.readTensor("p_smooths", p_smooths.data());

        // Split the hits by plane. (For consistency with the CPU code.)
        for (int i=0; i<this->num_planes; ++i)
            p_smoothed.push_back(p_smooths.block(4*i, 0, 4, this->hits_per_tile));
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

    // Create a compute set.
    poplar::ComputeSet computeSet = this->graph.addComputeSet("computeSet");

    // Tally counter for hits that have been processed so far.
    int hit_tally = 0;

    // Work out the number of hardware threads per-tile.
    auto num_workers = this->device.getTarget().getNumWorkerContexts();

    // Loop over the number of IPU tiles.
    for (int i=0; i<this->num_tiles; ++i)
    {
        // Create getNumWorkerContexts threads per tile.
        for (int j=0; j<num_workers; ++j)
        {
            // Work out the number of hits per thread.
            int hits_per_thread = int(this->hits_per_tile / num_workers);

            // If this is the first thread, then batch any extra tracks from
            // the remainder of the division.
            if (j == 0)
            {
                int remainder = int(this->hits_per_tile - 6*hits_per_thread);
                hits_per_thread += remainder;
            }

            // Only continue if there are hits to process.
            if (hits_per_thread > 0)
            {
                // Initalise state of system: { x, tan(theta), y, tan(phi) }
                // Use {x, 0, y, 0} column vectors for the initial hits. The matrix is
                // packed as hits along the columns and planes along the rows. (Every
                // fourth row is the start of a new plane.)
                MatrixRowMajorXf p0 = MatrixRowMajorXf::Zero(4*this->num_planes, hits_per_thread);
                for (int k=0; k<this->num_planes; ++k)
                {
                    p0.row(k*4 + 0) = this->hits_plane[k].col(0).transpose().segment(hit_tally, hits_per_thread);
                    p0.row(k*4 + 2) = this->hits_plane[k].col(1).transpose().segment(hit_tally, hits_per_thread);
                }

                // Update the hit tally.
                hit_tally += hits_per_thread;

                // Number of rows and columns for packed matrices.
                auto nRows = ulong(4*this->num_planes);
                auto nCols = ulong(hits_per_thread);

                // Create the vertex identification string.
                std::string vertex_id = std::to_string(i) + std::to_string(j);

                // Create the graph variables.
                poplar::Tensor v_p = this->graph.addVariable(poplar::FLOAT, {4, nCols}, "p" + vertex_id);
                poplar::Tensor v_m = this->graph.addVariable(poplar::FLOAT, {4, nCols}, "m" + vertex_id);
                poplar::Tensor v_p0 = this->graph.addVariable(poplar::FLOAT, {nRows, nCols}, "p0" + vertex_id);
                poplar::Tensor v_p_proj = this->graph.addVariable(poplar::FLOAT, {4, nCols}, "p_proj" + vertex_id);
                poplar::Tensor v_p_projs = this->graph.addVariable(poplar::FLOAT, {nRows, nCols}, "p_projs" + vertex_id);
                poplar::Tensor v_p_filt = this->graph.addVariable(poplar::FLOAT, {4, nCols}, "p_filt" + vertex_id);
                poplar::Tensor v_p_filts = this->graph.addVariable(poplar::FLOAT, {nRows, nCols}, "p_filts" + vertex_id);
                poplar::Tensor v_p_smooth = this->graph.addVariable(poplar::FLOAT, {4, nCols}, "p_smooth" + vertex_id);
                poplar::Tensor v_p_smooths = this->graph.addVariable(poplar::FLOAT, {nRows, nCols}, "p_smooths" + vertex_id);
                poplar::Tensor v_F = this->graph.addVariable(poplar::FLOAT, {4, 4}, "F" + vertex_id);
                poplar::Tensor v_C = this->graph.addVariable(poplar::FLOAT, {4, 4}, "C0" + vertex_id);
                poplar::Tensor v_C0 = this->graph.addVariable(poplar::FLOAT, {4, 4}, "C0" + vertex_id);
                poplar::Tensor v_C_proj = this->graph.addVariable(poplar::FLOAT, {4, 4}, "C_proj" + vertex_id);
                poplar::Tensor v_C_projs = this->graph.addVariable(poplar::FLOAT, {nRows, 4}, "C_projs" + vertex_id);
                poplar::Tensor v_C_filt = this->graph.addVariable(poplar::FLOAT, {4, 4}, "C_filt" + vertex_id);
                poplar::Tensor v_C_filts = this->graph.addVariable(poplar::FLOAT, {nRows, 4}, "C_filts" + vertex_id);
                poplar::Tensor v_FT = this->graph.addVariable(poplar::FLOAT, {4, 4}, "FT" + vertex_id);
                poplar::Tensor v_HTG = this->graph.addVariable(poplar::FLOAT, {4, 4}, "HTG" + vertex_id);
                poplar::Tensor v_HTGH = this->graph.addVariable(poplar::FLOAT, {4, 4}, "HTGH" + vertex_id);
                poplar::Tensor v_tmp_4x4_0 = this->graph.addVariable(poplar::FLOAT, {4, 4}, "tmp_4x4_0" + vertex_id);
                poplar::Tensor v_tmp_4x4_1 = this->graph.addVariable(poplar::FLOAT, {4, 4}, "tmp_4x4_1" + vertex_id);
                poplar::Tensor v_tmp_4x4_2 = this->graph.addVariable(poplar::FLOAT, {4, 4}, "tmp_4x4_2" + vertex_id);
                poplar::Tensor v_tmp_4xN_0 = this->graph.addVariable(poplar::FLOAT, {4, nCols}, "tmp_4xN_0" + vertex_id);
                poplar::Tensor v_tmp_4xN_1 = this->graph.addVariable(poplar::FLOAT, {4, nCols}, "tmp_4xN_1" + vertex_id);
                this->graph.setTileMapping(v_p, i);
                this->graph.setTileMapping(v_m, i);
                this->graph.setTileMapping(v_p0, i);
                this->graph.setTileMapping(v_p_proj, i);
                this->graph.setTileMapping(v_p_projs, i);
                this->graph.setTileMapping(v_p_filt, i);
                this->graph.setTileMapping(v_p_filts, i);
                this->graph.setTileMapping(v_p_smooth, i);
                this->graph.setTileMapping(v_p_smooths, i);
                this->graph.setTileMapping(v_F, i);
                this->graph.setTileMapping(v_C, i);
                this->graph.setTileMapping(v_C0, i);
                this->graph.setTileMapping(v_C_proj, i);
                this->graph.setTileMapping(v_C_projs, i);
                this->graph.setTileMapping(v_C_filt, i);
                this->graph.setTileMapping(v_C_filts, i);
                this->graph.setTileMapping(v_FT, i);
                this->graph.setTileMapping(v_HTG, i);
                this->graph.setTileMapping(v_HTGH, i);
                this->graph.setTileMapping(v_tmp_4x4_0, i);
                this->graph.setTileMapping(v_tmp_4x4_1, i);
                this->graph.setTileMapping(v_tmp_4x4_2, i);
                this->graph.setTileMapping(v_tmp_4xN_0, i);
                this->graph.setTileMapping(v_tmp_4xN_1, i);

                // Create the graph constants.
                poplar::Tensor c_p0 = this->graph.addConstant<float>(poplar::FLOAT, {nRows, nCols}, p0.data());
                poplar::Tensor c_F = this->graph.addConstant<float>(poplar::FLOAT, {4, 4}, F.data());
                poplar::Tensor c_C0 = this->graph.addConstant<float>(poplar::FLOAT, {4, 4}, C0.data());
                poplar::Tensor c_FT = this->graph.addConstant<float>(poplar::FLOAT, {4, 4}, FT.data());
                poplar::Tensor c_HTG = this->graph.addConstant<float>(poplar::FLOAT, {4, 4}, HTG.eval().data());
                poplar::Tensor c_HTGH = this->graph.addConstant<float>(poplar::FLOAT, {4, 4}, HTGH.eval().data());
                this->graph.setTileMapping(c_p0, i);
                this->graph.setTileMapping(c_F, i);
                this->graph.setTileMapping(c_C0, i);
                this->graph.setTileMapping(c_FT, i);
                this->graph.setTileMapping(c_HTG, i);
                this->graph.setTileMapping(c_HTGH, i);

                // Copy the inital data.
                this->prog.add(poplar::program::Copy(c_p0, v_p0));
                this->prog.add(poplar::program::Copy(c_F, v_F));
                this->prog.add(poplar::program::Copy(c_C0, v_C0));
                this->prog.add(poplar::program::Copy(c_FT, v_FT));
                this->prog.add(poplar::program::Copy(c_HTG, v_HTG));
                this->prog.add(poplar::program::Copy(c_HTGH, v_HTGH));

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
                this->graph.setTileMapping(vtx, i);

                // Create a read to the smoothed tracks for testing.
                if (i == 0 && j == 0)
                    this->graph.createHostRead("p_smooths", v_p_smooths);
            }
        }
    }

    // Add the compute set to the program.
    this->prog.add(poplar::program::Execute(computeSet));
}
