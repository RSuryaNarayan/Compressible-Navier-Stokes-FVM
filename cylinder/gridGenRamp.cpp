#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <fstream>
#include <sstream> // Include for string stream
#include <string>  // Include for std::string

struct Node {
    long double x, y; // Position of the node
};

class GridGenerator {
public:
    GridGenerator(int nx, int ny, long double length, long double inletHeight, long double rampAngle, long double rampStart, long double rampEnd)
      : nx(nx), ny(ny), length(length), inletHeight(inletHeight), rampAngle(rampAngle), rampStart(rampStart), rampEnd(rampEnd) {
        nodes.resize(nx * ny);
        generateGrid();
    }

    void generateGrid() {
        long double dx = length / (nx-1);

        //set lower boundary, j=0
        for (int i=0; i<nx; i++)
        {
            Node& node = nodes[i * ny];
            node.x = i*dx; 
            node.y = node.x<=rampStart? 0 : (node.x>=rampEnd? (rampEnd-rampStart) * tan(rampAngle) : (node.x-rampStart) * tan(rampAngle));
        }

        //assign grid nodes
        for (int i = 0; i < nx; ++i) {
            for (int j = 1; j < ny; ++j) {
                Node& node = nodes[i * ny + j];
                node.x = i * dx;
                node.y = nodes[i*ny].y + ((inletHeight-nodes[i*ny].y)/(ny-1)) * j; // Adjust for ramp slope
            }
        }
    }

    void outputToCSV(const std::string& filename) const {
        std::ofstream file(filename);

        if (file.is_open()) {
            file << "x,y\n"; // Updated CSV header for cell info

            for (const auto& node : nodes) {
                file << node.x << "," << node.y << "\n"; // Output node info
            }

            file.close();
            std::cout << "Grid data written to " << filename << std::endl;
        } else {
            std::cerr << "Unable to open file for writing: " << filename << std::endl;
        }
    }

private:
    int nx, ny;
    long double length, inletHeight, rampAngle, rampStart, rampEnd;
    std::vector<Node> nodes;
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1; // Exit if no input file is provided
    }

    // Coarse grid parameters
    int nx, ny;
    long double length, inletHeight, rampAngleDegrees, rampStart, rampEnd;
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
                if (key == "nx") nx = std::stoi(value);
                else if (key == "ny") ny = std::stoi(value);
                else if (key == "length") length = std::stold(value);
                else if (key == "inletHeight") inletHeight = std::stold(value);
                else if (key == "rampAngle") rampAngleDegrees = std::stold(value); // Read as degrees
                else if (key == "rampStart") rampStart = std::stold(value);
                else if (key == "rampEnd") rampEnd = std::stold(value);
                else if (key == "gridFile") csvFilename = value;
            }
        }
        inputFile.close();
    } else {
        std::cerr << "Unable to open input file: " << argv[1] << std::endl;
        return 1; // Exit if file cannot be opened
    }

    // Convert rampAngle from degrees to radians
    long double rampAngle = rampAngleDegrees * (M_PI / 180.0);

    GridGenerator Grid(nx, ny, length, inletHeight, rampAngle, rampStart, rampEnd);
    Grid.outputToCSV(csvFilename);  // Outputs grid data to CSV file
    
    return 0;
}
