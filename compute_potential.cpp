// Last update: 2024/09/26
// Author: Ricard Arbat Carandell

// Master in Aerospace Engineering - Computational Engineering
// Universitat Politècnica de Catalunya (UPC)
// Overview: Stream computing functions

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "mesh.cpp"
using namespace std;

struct Coefficients
{
    double C_L, C_D;
};

/**
 * Compares two matrices and returns the maximum error between them.
 *
 * @param v1 Matrix (Vector of vectors) to compare
 * @param v2 Matrix (Vector of vectors) to compare
 */
double streamsError(vector<vector<double>> v1, vector<vector<double>> v2)
{
    double maxError = 0.0;
    for (int i = 0; i < v1.size(); i++)
    {
        for (int j = 0; j < v1[i].size(); j++)
        {
            double diff = abs(v1[i][j] - v2[i][j]);
            if (diff > maxError)
            {
                maxError = diff;
            }
        }
    }
    return maxError;
}

/**
 * Fills the initial stream values of the mesh nodes.
 *
 * @param initial_stream_value Matrix (Vector of vectors) containing the next stream values
 * @param mesh Mesh matrix (vector of vectors) to be filled with Node structs
 * @param p Parameters of the simulation
 */
void fillStream(vector<vector<double>> &initial_stream_value, vector<vector<node>> &mesh, Parameters p)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            if (mesh[i][j].is_solid == true)
            {
                initial_stream_value[i][j] = p.solid_stream;
            }
            else
            {
                initial_stream_value[i][j] = p.start_stream;
            }
        }
    }
}

/**
 * Computes the stream function of the mesh nodes. It uses the discretized stream function
 * equation to solve the stream values of the nodes. The function iterates until the error
 * is below a certain threshold.
 *
 * @param mesh Mesh matrix (vector of vectors) to be filled with Node structs
 * @param p Parameters of the simulation
 */
void computeStream(vector<vector<node>> &mesh, Parameters p)
{
    double a_n, a_e, a_s, a_w, a_p, b_p;
    double d_PE, d_Pe, d_Ee, d_PS, d_Ps, d_Ss, d_PW, d_Pw, d_Ww, d_PN, d_Pn, d_Nn;

    vector<vector<double>> next_stream_value(N, vector<double>(M));
    vector<vector<double>> last_stream_value(N, vector<double>(M));

    fillStream(next_stream_value, mesh, p);

    double error = p.initial_error;
    while (error > p.delta)
    {
        last_stream_value = next_stream_value;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {

                if (j == 0) // Bottom channel nodes
                {
                    next_stream_value[i][j] = 0;
                }
                else if (j == M - 1) // Top channel nodes
                {
                    next_stream_value[i][j] = p.v_in * p.H;
                }
                else if (i == 0) // Inlet nodes
                {
                    next_stream_value[i][j] = p.v_in * mesh[i][j].y;
                }
                else if (i == N - 1) // Outlet nodes
                {
                    next_stream_value[i][j] = last_stream_value[i - 1][j];
                }
                else // Internal nodes
                {
                    // East node
                    d_PE = p.dx;
                    d_Pe = p.dx / 2;
                    d_Ee = p.dx / 2;
                    a_e = (d_PE / ((d_Pe * mesh[i][j].rho + d_Ee * mesh[i + 1][j].rho) / p.rho_in)) * p.dy / d_PE;

                    // South node
                    d_PS = p.dy;
                    d_Ps = p.dy / 2;
                    d_Ss = p.dy / 2;
                    a_s = (d_PS / ((d_Ps * mesh[i][j].rho + d_Ss * mesh[i][j - 1].rho) / p.rho_in)) * p.dx / d_PS;

                    // West node
                    d_PW = p.dx;
                    d_Pw = p.dx / 2;
                    d_Ww = p.dx / 2;
                    a_w = (d_PW / ((d_Pw * mesh[i][j].rho + d_Ww * mesh[i - 1][j].rho) / p.rho_in)) * p.dy / d_PW;

                    // North node
                    d_PN = p.dy;
                    d_Pn = p.dy / 2;
                    d_Nn = p.dy / 2;
                    a_n = (d_PN / ((d_Pn * mesh[i][j].rho + d_Nn * mesh[i][j + 1].rho) / p.rho_in)) * p.dx / d_PN;

                    // Discretization of the stream function equation
                    a_p = a_n + a_e + a_s + a_w;
                    next_stream_value[i][j] = (last_stream_value[i + 1][j] * a_e + last_stream_value[i - 1][j] * a_w + last_stream_value[i][j + 1] * a_n + last_stream_value[i][j - 1] * a_s + b_p) / a_p;
                }
            }
        }
        error = streamsError(next_stream_value, last_stream_value);
        printf("Error = %f\n", error);
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            mesh[i][j].stream = next_stream_value[i][j];
        }
    }
}

/**
 * Calculates the velocity of the mesh nodes. It uses the stream function values
 * to calculate the velocity components of the nodes.
 *
 * @param mesh Mesh matrix (vector of vectors) to be filled with Node structs
 * @param p Parameters of the simulation
 */
void calculateVelocity(vector<vector<node>> &mesh, Parameters p)
{
    double vxn, vye, vxs, vyw;
    double d_PN, d_Pn, d_Nn, d_PE, d_Pe, d_Ee, d_PS, d_Ps, d_Ss, d_PW, d_Pw, d_Ww;
    for (int i = 1; i < N - 1; i++)
    {
        for (int j = 1; j < M - 1; j++)
        {

            d_PE = p.dy;
            d_Pe = p.dy / 2;
            d_Ee = p.dy / 2;
            vye = -(d_PE / ((d_Pe * mesh[i][j].rho + d_Ee * mesh[i + 1][j].rho) / p.rho_in)) * ((mesh[i + 1][j].stream - mesh[i][j].stream) / d_PE);

            d_PS = p.dx;
            d_Ps = p.dx / 2;
            d_Ss = p.dx / 2;
            vxs = -(d_PS / ((d_Ps * mesh[i][j].rho + d_Ss * mesh[i][j - 1].rho) / p.rho_in)) * ((mesh[i][j - 1].stream - mesh[i][j].stream) / d_PS);

            d_PW = p.dy;
            d_Pw = p.dy / 2;
            d_Ww = p.dy / 2;
            vyw = (d_PW / ((d_Pw * mesh[i][j].rho + d_Ww * mesh[i - 1][j].rho) / p.rho_in)) * ((mesh[i - 1][j].stream - mesh[i][j].stream) / d_PW);

            // North node
            d_PN = p.dx;
            d_Pn = p.dx / 2;
            d_Nn = p.dx / 2;
            vxn = (d_PN / ((d_Pn * mesh[i][j].rho + d_Nn * mesh[i][j + 1].rho) / p.rho_in)) * ((mesh[i][j + 1].stream - mesh[i][j].stream) / d_PN);

            mesh[i][j].u = (vxn + vxs) / 2;
            mesh[i][j].v = (vye + vyw) / 2;
        }
    }
}

/**
 * Calculates the circulation around the cylinder. It uses the velocity values
 * to calculate the circulation of each cell around the cylinder, adding all of them up.
 *
 * @param mesh Mesh matrix (vector of vectors) to be filled with Node structs
 * @param p Parameters of the simulation
 */
double cylinderCirculation(vector<vector<node>> &mesh, Parameters p)
{
    double vxn, vye, vxs, vyw, d_PN, d_Pn, d_Nn, d_PE, d_Pe, d_Ee, d_PS, d_Ps, d_Ss, d_PW, d_Pw, d_Ww;
    double circulation = 0;
    for (int i = 1; i < N - 1; i++)
    {
        for (int j = 1; j < M - 1; j++)
        {
            d_PN = p.dy;
            d_Pn = p.dy / 2;
            d_Nn = p.dy / 2;
            vxn = (d_PN / ((d_Pn * mesh[i][j].rho + d_Nn * mesh[i][j + 1].rho / p.rho_in))) * ((mesh[i][j + 1].stream - mesh[i][j].stream) / d_PN);

            d_PE = p.dx;
            d_Pe = p.dx / 2;
            d_Ee = p.dx / 2;
            vye = -(d_PE / ((d_Pe * mesh[i][j].rho + d_Ee * mesh[i + 1][j].rho / p.rho_in))) * ((mesh[i + 1][j].stream - mesh[i][j].stream) / d_PE);

            d_PS = p.dy;
            d_Ps = p.dy / 2;
            d_Ss = p.dy / 2;
            vxs = -(d_PS / ((d_Ps * mesh[i][j].rho + d_Ss * mesh[i][j - 1].rho / p.rho_in))) * ((mesh[i][j - 1].stream - mesh[i][j].stream) / d_PS);

            d_PW = p.dx;
            d_Pw = p.dx / 2;
            d_Ww = p.dx / 2;
            vyw = (d_PW / ((d_Pw * mesh[i][j].rho + d_Ww * mesh[i - 1][j].rho / p.rho_in))) * ((mesh[i - 1][j].stream - mesh[i][j].stream) / d_PW);

            if (mesh[i][j].is_solid == true)
            {
                if (mesh[i][j - 1].is_solid == false)
                {
                    circulation += vxs * p.dx;
                }
                if (mesh[i][j + 1].is_solid == false)
                {
                    circulation += -vxn * p.dx;
                }
                if (mesh[i - 1][j].is_solid == false)
                {
                    circulation += -vyw * p.dy;
                }
                if (mesh[i + 1][j].is_solid == false)
                {
                    circulation += vye * p.dy;
                }
            }
        }
    }
    return circulation;
}

/**
 * Calculates the pressure coefficient of the mesh nodes. It uses the velocity values
 * to calculate the pressure coefficient of the nodes.
 *
 * @param mesh Mesh matrix (vector of vectors) to be filled with Node structs
 * @param p Parameters of the simulation
 */
void calculateCp(vector<vector<node>> &mesh, Parameters p)
{
    double v = 0;
    for (int i = 1; i < N - 1; i++)
    {
        for (int j = 1; j < M - 1; j++)
        {
            v = sqrt(mesh[i][j].u * mesh[i][j].u + mesh[i][j].v * mesh[i][j].v);
            mesh[i][j].cp = 1 - (v / p.v_in) * (v / p.v_in);
        }
    }
}

/**
 * Calculates the forces around the cylinder. It uses the pressure coefficient values
 * to calculate the lift and drag coefficients of the cylinder.
 *
 * @param mesh Mesh matrix (vector of vectors) to be filled with Node structs
 * @param p Parameters of the simulation
 */
struct Coefficients cylinderForces(vector<vector<node>> &mesh, Parameters p)
{
    Coefficients c;

    double q = 0.5 * p.rho_in * p.p_in * p.v_in;
    for (int i = 1; i < N - 1; i++)
    {
        for (int j = 1; j < M - 1; j++)
        {
            if (mesh[i][j].is_solid == true)
            {
                if (mesh[i][j - 1].is_solid == false)
                {
                    c.C_L += (mesh[i][j - 1].cp * q + p.p_in) * p.dx;
                }
                if (mesh[i][j + 1].is_solid == false)
                {
                    c.C_L += -(mesh[i][j + 1].cp * q + p.p_in) * p.dx;
                }
                if (mesh[i - 1][j].is_solid == false)
                {
                    c.C_D += (mesh[i - 1][j].cp * q + p.p_in) * p.dy;
                }
                if (mesh[i + 1][j].is_solid == false)
                {
                    c.C_D += -(mesh[i + 1][j].cp * q + p.p_in) * p.dy;
                }
            }
        }
    }

    c.C_L = c.C_L / (q * 2 * p.cylinder_r);
    c.C_D = c.C_D / (q * 2 * p.cylinder_r);

    return c;
}

/**
 * Converts Cartesian coordinates to polar coordinates.
 *
 * @param x The x-coordinate in Cartesian coordinates.
 * @param y The y-coordinate in Cartesian coordinates.
 * @return A pair where the first element is the radius (r) and the second element is the angle (theta) in radians.
 */
pair<double, double> cartesianToPolar(double x, double y)
{
    double r = sqrt(x * x + y * y);
    double theta = atan2(y, x);
    return make_pair(r, theta);
}

/**
 * Computes the stream function of the mesh nodes. It uses the discretized stream function
 * equation to solve the stream values of the nodes. The function iterates until the error
 * is below a certain threshold.
 *
 * @param mesh Mesh matrix (vector of vectors) to be filled with Node structs
 * @param p Parameters of the simulation
 */
void computeAnalyticStream(vector<vector<node>> &mesh, Parameters p)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            pair<double, double> polar;
            pair<double, double> cartesian;

            polar = cartesianToPolar(mesh[i][j].x - 0.5, mesh[i][j].y - 0.5);

            if (j == 0) // Bottom channel nodes
            {
                mesh[i][j].stream = 0;
            }
            else if (j == M - 1) // Top channel nodes
            {
                mesh[i][j].stream = p.v_in * p.H;
            }
            else if (i == 0) // Inlet nodes
            {
                mesh[i][j].stream = p.v_in * mesh[i][j].y;
            }
            else if (i == N - 1) // Outlet nodes
            {
                mesh[i][j].stream = mesh[i - 1][j].stream;
            }
            else // Internal nodes
            {
                mesh[i][j].stream = abs(p.v_in * (polar.first - (p.cylinder_r * p.cylinder_r) / polar.first) * sin(polar.second));
            }
        }
    }
}

/**
 * Compares two matrices and returns the maximum error between them.
 *
 * @param v1 Matrix (Vector of vectors) to compare
 * @param v2 Matrix (Vector of vectors) to compare
 */
double analyticError(vector<vector<node>> numerical, vector<vector<node>> analytical)
{
    double maxError = 0.0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            double diff = abs(numerical[i][j].stream - analytical[i][j].stream);
            if (diff > maxError)
            {
                maxError = diff;
            }
        }
    }
    return maxError;
}