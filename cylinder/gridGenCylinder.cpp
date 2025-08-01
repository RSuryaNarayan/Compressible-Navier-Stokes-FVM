#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>

struct Node {
    long double x, y; // Position of the node
    int index; // index of the node
};

struct Cell {
    // Geometric information
    long double x_c, y_c; // Position of the cell centers
    long double volume; // cell volume
    std::vector<long double> n_x; // unit normal along x for each face
    std::vector<long double> n_y; // unit normal along y for each face
    std::vector<long double> S; // face areas
    std::vector<int> faces_index; // indices of the faces bounding the cell
    long double faceVectorSum; // sum_faces ((nxS)^2+ (nyS)^2))
    std::vector<Node*> nodes; // Nodes of the cell in clockwise direction
    int index; // index of the current cell
    std::vector<int> neighbors; // vector storing indices of the neighboring cells
};

struct Face {
    int index; // index of the face
    std::vector<Node*> nodes; // Nodes of the face in clockwise direction
    long double S_f; // face area
    long double n_x; // unit normal along x for each face
    long double n_y; // unit normal along y for each face
    int i_L; // index of the cell to the left of the face
    int i_R; // index of the cell to the right of the face
    int bc_index;   // bc info
                    // 0: internal cell
                    // 1: left end (inlet)
                    // 2: right end (wall)
                    // 3: bottom end (outlet)
                    // 4: top end (outlet)
};

class CircularGridAroundCircle {
public:
    CircularGridAroundCircle(int nx, int ny, double r_circle, double r_x, double r_y, double stretchFactor = 1.2)
        : nx(nx), ny(ny), r_circle(r_circle), r_x(r_y), r_y(r_x), stretchFactor(stretchFactor) {

        ncells = ((nx - 1) * (ny - 1));
        nodes.resize(nx * ny);
        generateGrid();
    }

    void generateGrid() {

        // we want to generate points so that the ij-ordering in the computational space is retained
        // as the ramp case. Following this, we note that the outermost boundary needs to be filled first 
        // traversing all the way to the circle. similarly theta runs clockwise from bottom to top 

        // essentially this means dx is stretched to get closer resolution near the cylinder
        // dy is uniform and is decided by theta

        //without stretching
        long double dTheta = M_PI / (ny-1) / 2.0L;
        for (int i = 0; i < nx; ++i) { 
            for (int j = 0; j < ny; ++j) {
                
                //create node
                Node& node = nodes[i*ny + j];
                node.index = i*ny +j;
                //compute params for getting x and y
                long double theta =  j*dTheta + 0.5L * M_PI;

                // long double R_xi = r_circle + (nx-1-i)*(r_x-r_circle)/(nx-1); // no stretching 
                // long double R_yi = r_circle + (ny-1-i)*(r_y-r_circle)/(nx-1); // no stretching

                long double R_xi = r_circle + (r_x-r_circle) * (1.0L - pow(stretchFactor, (long double) (nx-1-i) ))/ (1.0L - pow(stretchFactor, (long double) (nx-1)));
                long double R_yi = r_circle + (r_y-r_circle) * (1.0L - pow(stretchFactor, (long double) (nx-1-i) ))/ (1.0L - pow(stretchFactor, (long double) (nx-1)));
                
                long double R_i = R_xi * R_yi / std::sqrt( R_xi*R_xi*sin(theta)*sin(theta) + R_yi*R_yi*cos(theta)*cos(theta) );

                //assign
                node.x = -1.0L * R_i * sin(M_PI-theta);
                node.y =  1.0L * R_i * cos(M_PI-theta);
            }
        }

        //convert to cells
        #include "convertToCells.H"

        //fill face info
        #include "nodesToFaces.H"

        //fill connectivity
        #include "connectivity.H"
    }

    void outputToCSV(const std::string& filename) const {
        std::ofstream file(filename);
        
        if (file.is_open()) {
            file << "x,y,index\n"; // CSV header

            for (const auto& node : nodes) {
                file << node.x << "," << node.y << "," << node.index << "\n";
            }

            file.close();
            std::cout << "Grid data written to " << filename << std::endl;
        } else {
            std::cerr << "Unable to open file for writing: " << filename << std::endl;
        }
    }

private:
    // int nr, nTheta;
    int nx, ny; // Radial and angular resolution
    double r_circle; // Radius of the inner circle (starting boundary)
    double r_x; // Outer radius of the grid along x
    double r_y; // Outer radius of the grid along y
    double stretchFactor; // Factor to control grid stretching

    //unstructured representation
    int nfaces, ncells;
    std::vector<Node> nodes;
    std::vector<Cell> cells;
    std::vector<Face> faces;
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1; // Exit if no input file is provided
    }

    // Coarse grid parameters
    int nr, nTheta;
    double r_circle, r_x, r_y;
    std::string csvFilename;

    // Read parameters from input file specified in command line
    std::ifstream inputFile(argv[1]);
    if (inputFile.is_open()) {
        std::string line;
        while (std::getline(inputFile, line)) {
            std::istringstream iss(line);
            std::string key, value;
            if (std::getline(iss, key, '=') && std::getline(iss, value)) {
                // Remove leading/trailing whitespace
                key.erase(0, key.find_first_not_of(" \t"));
                key.erase(key.find_last_not_of(" \t") + 1);
                value.erase(0, value.find_first_not_of(" \t"));
                value.erase(value.find_last_not_of(" \t") + 1);

                // Assign values based on keys
                if (key == "nx") nr = std::stoi(value);
                else if (key == "ny") nTheta = std::stoi(value);
                else if (key == "r_circle") r_circle = std::stold(value);
                else if (key == "r_x") r_x = std::stold(value);
                else if (key == "r_y") r_y = std::stold(value); // Read as degrees
                else if (key == "gridFile") csvFilename = value;
            }
        }
        inputFile.close();
    } else {
        std::cerr << "Unable to open input file: " << argv[1] << std::endl;
        return 1; // Exit if file cannot be opened
    }

    CircularGridAroundCircle circularGrid(nr, nTheta, r_circle, r_x, r_y);
    circularGrid.outputToCSV(csvFilename);  // Outputs grid data to CSV file

    return 0;
}
