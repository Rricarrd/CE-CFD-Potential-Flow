  // Compute analytic stream function
    vector<vector<node>> analytic_mesh(N, vector<node>(M));
    buildMesh(analytic_mesh, p);             // creating the mesh
    computeAnalyticStream(analytic_mesh, p); // stream solver