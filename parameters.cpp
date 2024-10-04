#include <iostream>
#include <cmath>
#include <fstream>
#include <array>

const int N = 100; // number of rows
const int M = N;   // number of columns

/**
 * Struct containing the parameters of the simulation. These values are used
 * by other functions to calculate the results of the simulation.
 *
 */

struct Parameters
{
    double L = 5;                                // length of the channel
    double H = 5;                                // height of the channel
    double cylinder_x = L / 2;                   // x position of the cylinder
    double cylinder_y = H / 2;                   // y position of the cylinder
    double cylinder_r = 0.1;                     // radius of the cylinder
    double p_in = 100000;                        // inlet pressure
    double t_in = 298;                           // inlet temperature
    double v_in = 10;                            // inlet velocity
    double rho_in = 1.225;                       // inlet density
    double initial_density = 1.225;              // initial density value
    double solid_density = pow(10, -10);         // solid density (very small value at solid nodes)
    double solid_stream = H * 0.5 * v_in;        // stream value at solid nodes
    double start_stream = 10;                    // initial stream value
    double delta = 0.00001;                      // maximum delta for the error
    double initial_error = 1.5;                  // initial error
    double dx = L / N;                           // x step
    double dy = L / M;                           // ystep
    double R = 287;                              // gas constant
    double gamma = 1.4;                          // specific heat ratio
    double specific_heat = 1005;                 // specific heat
    std::string folder = "output_final/";        // output file name
    std::string mesh_number = std::to_string(N); // mesh number
};