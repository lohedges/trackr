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

#include <poplar/Vertex.hpp>

// Handy typedefs for input and in-out types. (These are rank-2 tensors.)
using InputFloat = poplar::Vector<poplar::Input<poplar::Vector<float>>>;
using InOutFloat = poplar::Vector<poplar::InOut<poplar::Vector<float>>>;

// Unoptimised helper functions prototypes.

// Copy the tensor 'in' into 'out'. The offsets allow us to specify the row
// index at which we begin reading or writing in the respective tensor.
void copy(InputFloat &in, InOutFloat &out, int offset_in, int offset_out);
void copy(InOutFloat &in, InOutFloat &out, int offset_in, int offset_out);

// Sum tensors 'in0' and 'in1', placing the result in 'out'.
void sum(InputFloat &in0, InOutFloat &in1, InOutFloat &out);
void sum(InOutFloat &in0, InOutFloat &in1, InOutFloat &out);

// Subtract tensors 'in0' and 'in1', placing the result in 'out'.
void sub(InputFloat &in0, InOutFloat &in1, InOutFloat &out);
void sub(InOutFloat &in0, InOutFloat &in1, InOutFloat &out);

// Multiply matrices 'in0' and 'in1', placing the result in 'out'.
void mul(InputFloat &in0, InOutFloat &in1, InOutFloat &out);
void mul(InOutFloat &in0, InputFloat &in1, InOutFloat &out);
void mul(InOutFloat &in0, InOutFloat &in1, InOutFloat &out);

// Compute the inverse of matrix 'in', placing the result in 'out'.
void inverse(InputFloat &in, InOutFloat &out);
void inverse(InOutFloat &in, InOutFloat &out);

// Kalman Filter Codelet.

class KalmanFilter : public poplar::Vertex
{
public:
    // Input fields. (constant)
    InputFloat p0;
    InputFloat F;
    InputFloat FT;
    InputFloat C0;
    InputFloat HTG;
    InputFloat HTGH;

    // InOut fields. (read/writeable)
    InOutFloat p;
    InOutFloat m;
    InOutFloat p_proj;
    InOutFloat p_projs;
    InOutFloat p_filt;
    InOutFloat p_filts;
    InOutFloat p_smooth;
    InOutFloat p_smooths;
    InOutFloat C;
    InOutFloat C_proj;
    InOutFloat C_projs;
    InOutFloat C_filt;
    InOutFloat C_filts;
    InOutFloat tmp_4x4_0;
    InOutFloat tmp_4x4_1;
    InOutFloat tmp_4x4_2;
    InOutFloat tmp_4xN_0;
    InOutFloat tmp_4xN_1;

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


// Unoptimised helper functions definitions.

void copy(InputFloat &in, InOutFloat &out, int offset_in, int offset_out)
{
    for (int i=0; i<4; ++i)
    {
        for (int j=0; j<in[i].size(); ++j)
        {
            out[i+offset_out][j] = in[i+offset_in][j];
        }
    }
}

void copy(InOutFloat &in, InOutFloat &out, int offset_in, int offset_out)
{
    for (int i=0; i<4; ++i)
    {
        for (int j=0; j<in[i].size(); ++j)
        {
            out[i+offset_out][j] = in[i+offset_in][j];
        }
    }
}

void sum(InputFloat &in0, InOutFloat &in1, InOutFloat &out)
{
    for (int i=0; i<4; ++i)
    {
        for (int j=0; j<in0[i].size(); ++j)
        {
            out[i][j] = in0[i][j] + in1[i][j];
        }
    }
}

void sum(InOutFloat &in0, InOutFloat &in1, InOutFloat &out)
{
    for (int i=0; i<4; ++i)
    {
        for (int j=0; j<in0[i].size(); ++j)
        {
            out[i][j] = in0[i][j] + in1[i][j];
        }
    }
}

void sub(InputFloat &in0, InOutFloat &in1, InOutFloat &out)
{
    for (int i=0; i<4; ++i)
    {
        for (int j=0; j<in0[i].size(); ++j)
        {
            out[i][j] = in0[i][j] - in1[i][j];
        }
    }
}

void sub(InOutFloat &in0, InOutFloat &in1, InOutFloat &out)
{
    for (int i=0; i<4; ++i)
    {
        for (int j=0; j<in0[i].size(); ++j)
        {
            out[i][j] = in0[i][j] - in1[i][j];
        }
    }
}

void mul(InputFloat &in0, InOutFloat &in1, InOutFloat &out)
{
    for (int i=0; i<4; ++i)
    {
        for (int j=0; j<in1[0].size(); ++j)
        {
            out[i][j] = 0;
            for (int k=0; k<4; ++k)
            {
                out[i][j] += in0[i][k] * in1[k][j];
            }
        }
    }
}

void mul(InOutFloat &in0, InputFloat &in1, InOutFloat &out)
{
    for (int i=0; i<4; ++i)
    {
        for (int j=0; j<in1[0].size(); ++j)
        {
            out[i][j] = 0;
            for (int k=0; k<4; ++k)
            {
                out[i][j] += in0[i][k] * in1[k][j];
            }
        }
    }
}

void mul(InOutFloat &in0, InOutFloat &in1, InOutFloat &out)
{
    for (int i=0; i<4; ++i)
    {
        for (int j=0; j<in1[0].size(); ++j)
        {
            out[i][j] = 0;
            for (int k=0; k<4; ++k)
            {
                out[i][j] += in0[i][k] * in1[k][j];
            }
        }
    }
}

void inverse(InputFloat &in, InOutFloat &out)
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

void inverse(InOutFloat &in, InOutFloat &out)
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
