// Last update: 2024/09/26
// Author: Ricard Arbat Carandell

// Master in Aerospace Engineering - Computational Engineering
// Universitat Politècnica de Catalunya (UPC) - BarcelonaTech
// Overview: Potential flow solution in a channel around a cilinder using the stream function formulation.

// Libraries
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "compute_potential.cpp"
using namespace std;

/**
 * Main function of the program. It initializes the mesh, computes the stream
 * function and calculates the forces around the cilinder.
 */
int main(void)
{
    // SI Units
    Parameters p;
    vector<vector<node>> mesh(N, vector<node>(M));
    buildMesh(mesh, p); // creating the mesh
    setStream(mesh, p); // setting initial values
    setRho(mesh, p);

    computeStream(mesh, p);                            // stream solver
    calculateVelocity(mesh, p);                        // calculating velocity
    calculateCp(mesh, p);                              // calculating pressure coeficient distribution
    double circulation = cylinderCirculation(mesh, p); // calculate circulation around cilinder
    struct Coefficients c = cylinderForces(mesh, p);   // calculate forces around cilinder
    exportData(mesh);                                  // export data to file output.dat

    // Print final results
    std::cout << "### POTENTIAL FLOW RESULTS ###" << std::endl;
    std::cout << "C_L = " << c.C_L << std::endl;
    std::cout << "C_D = " << c.C_D << std::endl;
    std::cout << "Cylinder cicrulation = " << circulation << std::endl;

    return 0;
}
