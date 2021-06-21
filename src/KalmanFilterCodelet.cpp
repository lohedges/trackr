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

#include <ipudef.h>
#include <poplar/Vertex.hpp>

// Byte alignment of input vectors.
#define ALIGN 8

// Loop unroll factor.
#ifdef TEST
    #define UNROLL 2
#else
    #define UNROLL 1
#endif

// Handy typedefs for Poplar Input and InOut types. (These are rank-2 tensors.)
// We use the VectorLayout::SPAN so that we have access to the .size() member to
// work out the input size dynamically. For a production run, with a fixed
// number of tracks, we could instead use VectorLayout::ONE_PTR. See:
// https://docs.graphcore.ai/projects/assembly-programming/en/latest/vertex_vectors.html
using InputFloatTensor =
    poplar::Vector<poplar::Input<
    poplar::Vector<float, poplar::VectorLayout::SPAN, ALIGN>>,
    poplar::VectorLayout::SPAN, ALIGN>;
using InOutFloatTensor =
    poplar::Vector<poplar::InOut<
    poplar::Vector<float, poplar::VectorLayout::SPAN, ALIGN>>,
    poplar::VectorLayout::SPAN, ALIGN>;

// Templated helper functions prototypes.

// Copy 'in' into 'out'. The offsets allow us to specify the row
// index at which we begin reading or writing in the respective tensor.
template <typename T0, typename T1>
void copy(T0 &in, T1 &out, int offset_in, int offset_out);
template <typename T0, typename T1>
void row_copy(const T0 *__restrict in, T1 *__restrict out, int size);

// Sum 'in0' and 'in1', placing the result in 'out'.
template <typename T0, typename T1>
void sum(T0 &in0, T1 &in1, InOutFloatTensor &out);
template <typename T0, typename T1>
void row_sum(const T0 *__restrict in0, const T0 *__restrict in1, T1 *__restrict out, int size);

// Subtract 'in0' and 'in1', placing the result in 'out'.
template <typename T0, typename T1>
void sub(T0 &in0, T1 &in1, InOutFloatTensor &out);
template <typename T0, typename T1>
void row_sub(const T0 *__restrict in0, const T0 *__restrict in1, T1 *__restrict out, int size);

// Multiply 'in0' and 'in1', placing the result in 'out'.
template <typename T0, typename T1>
void mul(T0 &in0, T1 &in1, InOutFloatTensor &out);

// Compute the inverse of 'in', placing the result in 'out'.
// (Here 'in' is always a 4x4 matrix.)
template <typename T0, typename T1>
void inverse(T0 &in, T1 &out);

// Kalman Filter Codelet.

class KalmanFilter : public poplar::Vertex
{
public:
    // Input fields. (constant)
    InputFloatTensor p0;              // Intial state. (Planes every 4 rows, hits along cols.)
    InputFloatTensor F;               // Transfer matrix.
    InputFloatTensor FT;              // Transpose of transfer matrix.
    InputFloatTensor C0;              // Covariance matrix.
    InputFloatTensor HTG;             // (H^T)*G
    InputFloatTensor HTGH;            // (H^T)*G*H

    // InOut fields. (read/writeable)
    InOutFloatTensor p;               // The state tensor for hits on the current plane.
    InOutFloatTensor m;               // Hits for the current plane.
    InOutFloatTensor p_proj;          // The projected states for the current plane.
    InOutFloatTensor p_projs;         // The projected states for all planes.
    InOutFloatTensor p_filt;          // The filtered states for the current plane.
    InOutFloatTensor p_filts;         // The filtered states for all planes.
    InOutFloatTensor p_smooth;        // The smoothed states for the current plane.
    InOutFloatTensor p_smooths;       // The smoothed states for all planes.
    InOutFloatTensor C;               // The covariance matrix for the current plane.
    InOutFloatTensor C_proj;          // The projected covariance for the current plane.
    InOutFloatTensor C_projs;         // The projected covariance for all planes.
    InOutFloatTensor C_filt;          // The filtered covariance matrix for the current plane.
    InOutFloatTensor C_filts;         // The filtered covariance matrix for all planes.
    InOutFloatTensor tmp_4x4_0;       // A temporary 4x4 tensor.
    InOutFloatTensor tmp_4x4_1;       // A temporary 4x4 tensor.
    InOutFloatTensor tmp_4x4_2;       // A temporary 4x4 tensor.
    InOutFloatTensor tmp_4xN_0;       // A temporary 4xN tensor.
    InOutFloatTensor tmp_4xN_1;       // A temporary 4xN tensor.

    // Overloaded compute function.
    bool compute()
    {
        // Work out the number of planes from the size of the input
        // state vector. Hits are packed by column, planes by rows.
        int num_planes = p0.size() / 4;

        // Copy the state vector for the first plane into p.
        copy(p0, p, 0, 0);

        // Initialise the covariance matrix.
        copy(C0, C, 0, 0);

        // Step forward through the planes to linearly project the state
        // and correct (filter) using the measured hits and covariance.
        for (int i=0; i<num_planes; ++i)
        {
            // p_proj = F * p
            mul(F, p, p_proj);
            copy(p_proj, p_projs, 0, i*4);

            // C_proj = F * C * FT
            mul(F, C, tmp_4x4_0);
            mul(tmp_4x4_0, FT, C_proj);
            copy(C_proj, C_projs, 0, i*4);

            // C_filt = (C_proj.inverse() + HTGH).inverse()
            inverse(C_proj, tmp_4x4_0);
            sum(HTGH, tmp_4x4_0, C_proj);
            inverse(C_proj, C_filt);
            copy(C_filt, C_filts, 0, i*4);

            // tmp = (C_proj.inverse() * p_proj) + (HTG * m)
            copy(p0, m, i*4, 0);
            mul(tmp_4x4_0, p_proj, tmp_4xN_0);
            mul(HTG, m, tmp_4xN_1);
            sum(tmp_4xN_0, tmp_4xN_1, tmp_4xN_1);

            // p_filt = C_filt * tmp
            mul(C_filt, tmp_4xN_1, p_filt);
            copy(p_filt, p_filts, 0, i*4);

            // p = p_filt
            copy(p_filt, p, 0, 0);

            // C = C_filt
            copy(C_filt, C, 0, 0);
        }

        // Initialise the smmoothed hits for the final plane.
        copy(p_filt, p_smooth, 0, 0);
        copy(p_smooth, p_smooths, 0, (num_planes-1)*4);

        // Now propagate state and covariance updates backwards through the
        // states to find the globally optimum configuration.
        for (int i=num_planes-2; i>=0; --i)
        {
            // A = C_filts[i] * FT * C_projs[i+1].inverse()
            copy(C_filts, C_filt, i*4, 0);
            copy(C_projs, C_proj, (i+1)*4, 0);
            inverse(C_proj, tmp_4x4_1);
            mul(C_filt, FT, tmp_4x4_0);
            mul(tmp_4x4_0, tmp_4x4_1, tmp_4x4_2);

            // p_smooth = p_filts[i] + A*(p_smoothed[i+1] - p_projs[i+1]);
            copy(p_filts, p_filt, i*4, 0);
            copy(p_projs, p_proj, (i+1)*4, 0);
            sub(p_smooth, p_proj, p_smooth);
            mul(tmp_4x4_2, p_smooth, tmp_4xN_0);
            sum(p_filt, tmp_4xN_0, p_smooth);
            copy(p_smooth, p_smooths, 0, i*4);
        }

        return true;
    }
};

// Templated helper functions definitions.
// (Note that these are currently unoptimised.)

template <typename T0, typename T1>
void copy(T0 &in, T1 &out, int offset_in, int offset_out)
{
    // Work out the number of hits. (Each row has the same number of columns.)
    // This must be a multiple of 8, which is validated elsewhere. We divide by
    // two since we cast to float2, i.e. float2 contains two floats.
    int size = in[0].size() / 2;

    for (int i=0; i<4; ++i)
    {
        // Cast to float2 to to instruct compiler to emit 64-bit wide,
        // aligned instructions for 32-bit elements.
        float2 *p_in = const_cast<float2 *>(reinterpret_cast<const float2 *>(&in[i+offset_in][0]));
        float2 *p_out = reinterpret_cast<float2 *>(&out[i+offset_out][0]);

        // Process the row.
        row_copy(p_in, p_out, size);
    }
}

template <typename T0, typename T1>
void row_copy(const T0 *__restrict in, T1 *__restrict out, int size)
{
    for (int i=0; i<size; i+=UNROLL)
    {
        #pragma unroll UNROLL
        for (int j=0; j<UNROLL; ++j)
            out[i*UNROLL + j] = in[i*UNROLL + j];
    }
}

template <typename T0, typename T1>
void sum(T0 &in0, T1 &in1, InOutFloatTensor &out)
{
    // Work out the number of hits. (Each row has the same number of columns.)
    // This must be a multiple of 8, which is validated elsewhere. We divide by
    // two since we cast to float2, i.e. float2 contains two floats.
    int size = in0[0].size() / 2;

    for (int i=0; i<4; ++i)
    {
        // Cast to float2 to to instruct compiler to emit 64-bit wide,
        // aligned instructions for 32-bit elements.
        float2 *p_in0 = const_cast<float2 *>(reinterpret_cast<const float2 *>(&in0[i][0]));
        float2 *p_in1 = const_cast<float2 *>(reinterpret_cast<const float2 *>(&in1[i][0]));
        float2 *p_out = reinterpret_cast<float2 *>(&out[i][0]);

        // Process the row.
        row_sum(p_in0, p_in1, p_out, size);
    }
}

template <typename T0, typename T1>
void row_sum(const T0 *__restrict in0, const T0 *__restrict in1, T1 *__restrict out, int size)
{
    for (int i=0; i<size; i+=UNROLL)
    {
        #pragma unroll UNROLL
        for (int j=0; j<UNROLL; ++j)
            out[i*UNROLL + j] = in0[i*UNROLL + j] + in1[i*UNROLL + j];
    }
}

template <typename T0, typename T1>
void sub(T0 &in0, T1 &in1, InOutFloatTensor &out)
{
    // Work out the number of hits. (Each row has the same number of columns.)
    // This must be a multiple of 8, which is validated elsewhere. We divide by
    // two since we cast to float2, i.e. float2 contains two floats.
    int size = in0[0].size() / 2;

    for (int i=0; i<4; ++i)
    {
        // Cast to float2 to to instruct compiler to emit 64-bit wide,
        // aligned instructions for 32-bit elements.
        float2 *p_in0 = const_cast<float2 *>(reinterpret_cast<const float2 *>(&in0[i][0]));
        float2 *p_in1 = const_cast<float2 *>(reinterpret_cast<const float2 *>(&in1[i][0]));
        float2 *p_out = reinterpret_cast<float2 *>(&out[i][0]);

        // Process the row.
        row_sub(p_in0, p_in1, p_out, size);
    }
}

template <typename T0, typename T1>
void row_sub(const T0 *__restrict in0, const T0 *__restrict in1, T1 *__restrict out, int size)
{
    for (int i=0; i<size; i+=UNROLL)
    {
        #pragma unroll UNROLL
        for (int j=0; j<UNROLL; ++j)
            out[i*UNROLL + j] = in0[i*UNROLL + j] - in1[i*UNROLL + j];
    }
}

template <typename T0, typename T1>
void mul(T0 &in0, T1 &in1, InOutFloatTensor &out)
{
    int size = in1[0].size();

    for (int i=0; i<4; ++i)
    {
        for (int j=0; j<size; ++j)
        {
            out[i][j] = 0;
            for (int k=0; k<4; ++k)
            {
                out[i][j] += in0[i][k] * in1[k][j];
            }
        }
    }
}

template <typename T0, typename T1>
void inverse(T0 &in, T1 &out)
{
    // Adapted from: https://stackoverflow.com/a/60374938

    float A2323 = in[2][2] * in[3][3] - in[2][3] * in[3][2];
    float A1323 = in[2][1] * in[3][3] - in[2][3] * in[3][1];
    float A1223 = in[2][1] * in[3][2] - in[2][2] * in[3][1];
    float A0323 = in[2][0] * in[3][3] - in[2][3] * in[3][0];
    float A0223 = in[2][0] * in[3][2] - in[2][2] * in[3][0];
    float A0123 = in[2][0] * in[3][1] - in[2][1] * in[3][0];
    float A2313 = in[1][2] * in[3][3] - in[1][3] * in[3][2];
    float A1313 = in[1][1] * in[3][3] - in[1][3] * in[3][1];
    float A1213 = in[1][1] * in[3][2] - in[1][2] * in[3][1];
    float A2312 = in[1][2] * in[2][3] - in[1][3] * in[2][2];
    float A1312 = in[1][1] * in[2][3] - in[1][3] * in[2][1];
    float A1212 = in[1][1] * in[2][2] - in[1][2] * in[2][1];
    float A0313 = in[1][0] * in[3][3] - in[1][3] * in[3][0];
    float A0213 = in[1][0] * in[3][2] - in[1][2] * in[3][0];
    float A0312 = in[1][0] * in[2][3] - in[1][3] * in[2][0];
    float A0212 = in[1][0] * in[2][2] - in[1][2] * in[2][0];
    float A0113 = in[1][0] * in[3][1] - in[1][1] * in[3][0];
    float A0112 = in[1][0] * in[2][1] - in[1][1] * in[2][0];

    float det = in[0][0] * ( in[1][1] * A2323 - in[1][2] * A1323 + in[1][3] * A1223 )
              - in[0][1] * ( in[1][0] * A2323 - in[1][2] * A0323 + in[1][3] * A0223 )
              + in[0][2] * ( in[1][0] * A1323 - in[1][1] * A0323 + in[1][3] * A0123 )
              - in[0][3] * ( in[1][0] * A1223 - in[1][1] * A0223 + in[1][2] * A0123 );

    det = 1 / det;

    out[0][0] = det *   ( in[1][1] * A2323 - in[1][2] * A1323 + in[1][3] * A1223 );
    out[0][1] = det * - ( in[0][1] * A2323 - in[0][2] * A1323 + in[0][3] * A1223 );
    out[0][2] = det *   ( in[0][1] * A2313 - in[0][2] * A1313 + in[0][3] * A1213 );
    out[0][3] = det * - ( in[0][1] * A2312 - in[0][2] * A1312 + in[0][3] * A1212 );
    out[1][0] = det * - ( in[1][0] * A2323 - in[1][2] * A0323 + in[1][3] * A0223 );
    out[1][1] = det *   ( in[0][0] * A2323 - in[0][2] * A0323 + in[0][3] * A0223 );
    out[1][2] = det * - ( in[0][0] * A2313 - in[0][2] * A0313 + in[0][3] * A0213 );
    out[1][3] = det *   ( in[0][0] * A2312 - in[0][2] * A0312 + in[0][3] * A0212 );
    out[2][0] = det *   ( in[1][0] * A1323 - in[1][1] * A0323 + in[1][3] * A0123 );
    out[2][1] = det * - ( in[0][0] * A1323 - in[0][1] * A0323 + in[0][3] * A0123 );
    out[2][2] = det *   ( in[0][0] * A1313 - in[0][1] * A0313 + in[0][3] * A0113 );
    out[2][3] = det * - ( in[0][0] * A1312 - in[0][1] * A0312 + in[0][3] * A0112 );
    out[3][0] = det * - ( in[1][0] * A1223 - in[1][1] * A0223 + in[1][2] * A0123 );
    out[3][1] = det *   ( in[0][0] * A1223 - in[0][1] * A0223 + in[0][2] * A0123 );
    out[3][2] = det * - ( in[0][0] * A1213 - in[0][1] * A0213 + in[0][2] * A0113 );
    out[3][3] = det *   ( in[0][0] * A1212 - in[0][1] * A0212 + in[0][2] * A0112 );
}
