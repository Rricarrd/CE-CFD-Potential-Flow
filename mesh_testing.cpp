// Last update: 2024/09/26
// Author: Ricard Arbat Carandell

// Master in Aerospace Engineering - Computational Engineering
// Universitat Politècnica de Catalunya (UPC) - BarcelonaTech
// Overview: Mesh testing

// Libraries
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "mesh.cpp"

using namespace std;

/**
 * Function used to test the mesh creation and the node structs.
 * Also tests the mesh export to a CSV file.
 */
int main(void)
{
    // SI Units
    Parameters p;
    vector<vector<node>> mesh(N, vector<node>(M));
    buildMesh(mesh, p); // creating the mesh
    setStream(mesh, p); // setting initial values
    setRho(mesh, p);
    exportData(mesh); // export data to file output.dat
}