#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>

struct Cell {
    double x, y;      // Position of the cell in Cartesian coordinates
    double area;      // Area of the cell
};

class CircularGridAroundCircle {
public:
    CircularGridAroundCircle(int nr, int nTheta, double r_circle, double r_x, double r_y, double stretchFactor = 2.0)
        : nr(nr), nTheta(nTheta), r_circle(r_circle), r_x(r_y), r_y(r_x), stretchFactor(stretchFactor) {
        cells.resize(nr * nTheta);
        generateGrid();
    }

    void generateGrid() {
        double dTheta = M_PI / (nTheta-1);           // Angular step size for a full circle

        for (int j = 0; j < nTheta; ++j) {
            double theta = j * dTheta;               // Angle for each point on the circle
            double r_outer = r_x * r_y / std::sqrt(r_x*r_x*sin(theta)*sin(theta)+r_y*r_y*cos(theta)*cos(theta));
            double dr = (r_outer - r_circle)/nr;
            for (int i = 0; i < nr; ++i) {
                double r = r_circle + (i + 0.5) * dr * std::pow((double)i / (nr - 1), stretchFactor); // Radial position with stretching
                Cell& cell = cells[i * nTheta + j];
            
                // Convert polar coordinates to Cartesian
                cell.y = r * cos(theta);
                cell.x = -1.0*r * sin(theta);
                
                // Calculate area of each cell (approximate as sector area)
                cell.area = r * dr * dTheta;
            }
        }
    }

    void outputToCSV(const std::string& filename) const {
        std::ofstream file(filename);
        
        if (file.is_open()) {
            file << "x,y,area\n"; // CSV header

            for (const auto& cell : cells) {
                file << cell.x << "," << cell.y << "," << cell.area << "\n";
            }

            file.close();
            std::cout << "Grid data written to " << filename << std::endl;
        } else {
            std::cerr << "Unable to open file for writing: " << filename << std::endl;
        }
    }

private:
    int nr, nTheta;          // Radial and angular resolution
    double r_circle;         // Radius of the inner circle (starting boundary)
    double r_x;              // Outer radius of the grid along x
    double r_y;
    double stretchFactor; // Factor to control grid stretching
    std::vector<Cell> cells;
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
