#include <iostream>
#include <vector>
#include <array>

struct Vertex
{
    double x, y; // Coordinates of the vertex
};

struct Face
{
    std::array<int, 2> vertices; // Indices of the vertices that form this face
    int leftCell, rightCell;     // Indices of neighboring cells (rightCell = -1 if boundary)

    Face(int v1, int v2, int lCell = -1, int rCell = -1)
        : vertices{v1, v2}, leftCell(lCell), rightCell(rCell) {}
};

struct Cell
{
    std::vector<int> vertices; // Indices of the vertices that form this cell
    std::vector<int> faces;    // Indices of the faces belonging to this cell
};

class Mesh
{
public:
    std::vector<Vertex> vertices; // List of all vertices
    std::vector<Face> faces;      // List of all faces
    std::vector<Cell> cells;      // List of all cells

    // Generate a simple structured mesh with nx x ny cells in a rectangular domain [0, width] x [0, height]
    void generateStructuredMesh(int nx, int ny, double width, double height);
    void printMeshInfo() const;
};

void Mesh::generateStructuredMesh(int nx, int ny, double width, double height)
{
    int numVerticesX = nx + 1;
    int numVerticesY = ny + 1;

    // Step 1: Create vertices
    for (int j = 0; j < numVerticesY; ++j)
    {
        for (int i = 0; i < numVerticesX; ++i)
        {
            double x = i * width / nx;
            double y = j * height / ny;
            vertices.push_back({x, y});
        }
    }

    // Step 2: Create cells and faces
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            int v1 = j * numVerticesX + i;
            int v2 = v1 + 1;
            int v3 = v1 + numVerticesX;
            int v4 = v3 + 1;

            // Create the cell (quadrilateral)
            Cell cell;
            cell.vertices = {v1, v2, v4, v3};
            int cellIndex = cells.size();
            cells.push_back(cell);

            // Create the faces for the cell and update neighbors
            // Left face (v1, v3)
            if (i == 0)
            {
                faces.push_back(Face(v1, v3, cellIndex, -1)); // Boundary face
                cell.faces.push_back(faces.size() - 1);
            }
            else
            {
                int leftCellIndex = j * nx + (i - 1);
                faces.push_back(Face(v1, v3, cellIndex, leftCellIndex));
                cells[leftCellIndex].faces.push_back(faces.size() - 1);
                cell.faces.push_back(faces.size() - 1);
            }

            // Bottom face (v1, v2)
            if (j == 0)
            {
                faces.push_back(Face(v1, v2, cellIndex, -1)); // Boundary face
                cell.faces.push_back(faces.size() - 1);
            }
            else
            {
                int bottomCellIndex = (j - 1) * nx + i;
                faces.push_back(Face(v1, v2, cellIndex, bottomCellIndex));
                cells[bottomCellIndex].faces.push_back(faces.size() - 1);
                cell.faces.push_back(faces.size() - 1);
            }

            // Right face (v2, v4)
            if (i == nx - 1)
            {
                faces.push_back(Face(v2, v4, cellIndex, -1)); // Boundary face
                cell.faces.push_back(faces.size() - 1);
            }

            // Top face (v3, v4)
            if (j == ny - 1)
            {
                faces.push_back(Face(v3, v4, cellIndex, -1)); // Boundary face
                cell.faces.push_back(faces.size() - 1);
            }
        }
    }
}

void Mesh::printMeshInfo() const
{
    std::cout << "Mesh Info:\n";
    std::cout << "Number of vertices: " << vertices.size() << "\n";
    std::cout << "Number of faces: " << faces.size() << "\n";
    std::cout << "Number of cells: " << cells.size() << "\n";

    std::cout << "\nVertices:\n";
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        std::cout << "Vertex " << i << ": (" << vertices[i].x << ", " << vertices[i].y << ")\n";
    }

    std::cout << "\nFaces:\n";
    for (size_t i = 0; i < faces.size(); ++i)
    {
        std::cout << "Face " << i << ": Vertices (" << faces[i].vertices[0] << ", " << faces[i].vertices[1] << "), ";
        std::cout << "Left Cell: " << faces[i].leftCell << ", Right Cell: " << faces[i].rightCell << "\n";
    }

    std::cout << "\nCells:\n";
    for (size_t i = 0; i < cells.size(); ++i)
    {
        std::cout << "Cell " << i << ": Vertices (";
        for (int v : cells[i].vertices)
        {
            std::cout << v << " ";
        }
        std::cout << "), Faces (";
        for (int f : cells[i].faces)
        {
            std::cout << f << " ";
        }
        std::cout << ")\n";
    }
}

int main()
{
    Mesh mesh;
    mesh.generateStructuredMesh(3, 3, 1.0, 1.0); // 3x3 cells in a 1x1 domain
    mesh.printMeshInfo();

    return 0;
}
