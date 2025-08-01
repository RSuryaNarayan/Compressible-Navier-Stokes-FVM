#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <fstream>
#include <sstream>
#include <string>  
#include <cstdlib> 
#include <Eigen/Dense> 
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <chrono> 

#define URHO 0
#define UMX 1
#define UMY 2
#define UEINT 3

#define QRHO 0
#define QU 1
#define QV 2
#define QP 3
#define QT 4

#define GAMMA 1.4
#define RU 8.314
#define R 287.05 

struct Node {
    long double x, y; // Position of the node
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

    //data that each face stores
    Eigen::Matrix<long double, 4, 1> F; // flux at the face
    Eigen::Matrix<long double, 4, 4> Ap; // Jacobian plus
    Eigen::Matrix<long double, 4, 4> Am; // Jacobian minus
    
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

    //data that each cell stores
    Eigen::Matrix<long double, 4, 1> U; // conserved state vector
    Eigen::Matrix<long double, 4, 1> dU_old; // conserved state vector
    Eigen::Matrix<long double, 4, 1> dU_new; // conserved state vector
    Eigen::Matrix<long double, 5, 1> V; // primitive state vector
    Eigen::Matrix<long double, 4, 1> Ri = Eigen::Matrix<long double, 4, 1>::Zero(); // residual vector containing fluxes
    Eigen::Matrix<long double, 4, 4> Di = Eigen::Matrix<long double, 4, 4>::Zero(); // diagonal vector 
    Eigen::Matrix<long double, 4, 4> Di_guess = Eigen::Matrix<long double, 4, 4>::Zero(); // diagonal vector 
    
    //neighboring jacobians for the cell
    std::vector<Eigen::Matrix<long double, 4, 4>>N;
    
};

class GridReader {
public:
    GridReader(int nx, int ny, std::string gridFile)
      : nx(nx), ny(ny), gridFile(gridFile) {

        ncells = ((nx - 1) * (ny - 1));
        nodes.resize(nx * ny);

        generateCells();
    }

    void generateCells() {

        //read the grid into nodes
        #include "readGrid.H"

        //convert to nodes to finite volume representation
        #include "convertToCells.H"

        // Now populate face information
        #include "nodesToFaces.H"

        //now populate the connectivity information into the Cell structure by assigning faces to cells
        #include "connectivity.H"

        // Can begin actual solution. First set initial condition
        #include "initialCondition.H"

        // fill in the other cell matrices
        dx = nodes[nx*ny-1].x-nodes[nx*ny-2].x; //rough estimate
        dt = CFL*dx/M_inlet/(GAMMA*R*T_inlet); // rough dt estimate
        std::cout<<"\nnodes[nx*ny-1].x = "<< nodes[nx*ny-1].x;
        std::cout<<"\nnodes[nx*ny-2].x = "<< nodes[nx*ny-2].x;
        std::cout<<"dx = "<<dx<<"\n";
        std::cout<<"dt = "<<dt<<"\n";
        long double time = 0.0;
        for (int k=1; k<=nsteps; k++)
        {
            time = time + dt;
            std::cout<<"TIMESTEP="<<k<<" TIME="<<time<<std::endl;
            //loop over cells to fill D with the volume term. 
            //All other terms need flux info
            for (auto& cell : cells)
            {
                for (int i=0; i<4; i++)
                {
                    cell.Di(i,i) = cell.volume/dt;
                    cell.Di_guess(i,i) = cell.volume/dt;
                    cell.Ri(i) = 0.0;
                }
            }

            //fill fluxes by looping over faces
            for (int iif=0; iif<nfaces; iif++)
            {
                Face& face = faces[iif];

                //interior face flux
                if (face.bc_index==0)
                {
                    #include "interiorFluxAtFace.H"
                }

                //left end of the domain (inlet flux)
                if (face.bc_index==1)
                {
                    #include "leftEnd.H"
                }

                //right end of the domain (wall flux)
                if (face.bc_index==2)
                {
                    #include "rightEnd.H"
                }

                // bottom end of the domain (symmetry flux)
                if (face.bc_index==3)
                {
                    #include "bottomEnd.H"
                }

                // top end of the domain (outlet flux)
                if (face.bc_index==4)
                {
                    #include "topEnd.H"
                }
            }

            // for (int iic=0; iic<ncells; iic++)
            // {
            //     Cell cell = cells[iic];
            //     std::cout<<"Cell "<<iic<<" has neighbors ";
            //     for (size_t inbr=0; inbr<cell.N.size(); inbr++)
            //     {
            //         std::cout<<cell.neighbors[inbr]<<" Jacobian="<<cell.N[inbr]<<"\n";
            //     }
            //     std::cout<<"\n";
            // }

            // for (auto& face : faces)
            // {
            //     std::cout<<"\n Face "<<face.index<<" has flux=\n "<<face.F;
            // }

            // for (auto& cell: cells)
            // {
            //    std::cout<<"\nCell "<<cell.index<<" Ri = \n"<<cell.Ri;
            // }

            //Implicit update step
            // #include "solveSystemSparse.H"

            //point implicit relaxation
            //explicit guess for initial dU
            for  (auto& cell : cells)
            {
                cell.dU_old = cell.Di_guess.inverse()*(cell.Ri);
            }
            long double error=10.0L;
            for (int itr=1; itr<=max_iters; itr++)
            {
                error=0;
                for (auto& cell : cells)
                {
                    Eigen::Matrix<long double, 4, 1> RHS = cell.Ri;
                    for (size_t inbr=0; inbr<cell.neighbors.size(); inbr++)
                    {
                       int nbr_idx = cell.neighbors[inbr]; 
                       RHS  = RHS - cell.N[inbr]*cells[nbr_idx].dU_old;
                    }
                    cell.dU_new = cell.Di.inverse()*RHS;
                    error = error + std::abs( cell.dU_new(0) - cell.dU_old(0)) 
                    + std::abs( cell.dU_new(1) - cell.dU_old(1)) + std::abs( cell.dU_new(2) - cell.dU_old(2)) +std::abs( cell.dU_new(3) - cell.dU_old(3));
                    cell.dU_old = cell.dU_new;  
                }
                // std::cout<<"STEP = "<<k<<" iter = "<<itr<<" error="<<error<<std::endl;
                if (error<=tol) {std::cout<<"converged in "<<itr<<" iters!"<<std::endl; break;}
            }

            long double res=0;
            for (auto& cell : cells)
            {
                
                for (int i=0; i<4; i++)
                {
                    cell.U(i) = cell.U(i) + cell.dU_new(i);
                }
                Eigen::Matrix<long double, 5, 1> V_new;
                U2V(cell.U, V_new);
                cell.V = V_new;

                //compute residuals
                res += cell.Ri(0)*cell.Ri(0)/cell.volume/cell.volume;

                //clear out the neighbors and their Jacobian matrices 
                cell.neighbors.clear();
                cell.neighbors.shrink_to_fit();
                cell.N.clear();
                cell.N.shrink_to_fit();
                
                // Reset Di and Di_guess to zero
                cell.Di.setZero(); // Reset Di to zero
                cell.Di_guess.setZero(); // Reset Di_guess to zero
            }
            res = std::sqrt(res);
            residuals.push_back(res);
            
            std::cout<<"residual = "<<res<<std::endl;
            //output the initial condition to check
            if (k % 10000 == 0) {
                outputPrimitiveStateToCSV("solution" + std::to_string(k) + ".csv");
            }
            std::cout<<"Finished STEP = "<<k<<"\n\n";
        }

        // After computing residuals
        outputResidualsToCSV("residuals.csv"); // Output residuals to CSV

    }

    void computeMinMaxFaceVectorSum() {
        long double minFaceVector = std::numeric_limits<long double>::max();
        long double maxFaceVector = std::numeric_limits<long double>::lowest();

        for (const auto& cell : cells) {
            if (cell.faceVectorSum < minFaceVector) minFaceVector = cell.faceVectorSum;
            if (cell.faceVectorSum > maxFaceVector) maxFaceVector = cell.faceVectorSum;
        }

        std::cout << "Min Face Vector Sum: " << minFaceVector << std::endl;
        std::cout << "Max Face Vector Sum: " << maxFaceVector << std::endl;
    }

    void outputToCSV(const std::string& filename) const {
        std::ofstream file(filename);

        if (file.is_open()) {  
            file << "x_c,y_c,volume,skewness\n"; // Updated CSV header for cell info

            for (const auto& cell : cells) {
                file << cell.x_c << "," << cell.y_c << "," << cell.volume << "," << cell.faceVectorSum << "\n"; // Output cell info
            }

            file.close();
            std::cout << "Grid data written to " << filename << std::endl;
        } else {
            std::cerr << "Unable to open file for writing: " << filename << std::endl;
        }
    }

    void outputPrimitiveStateToCSV(const std::string& filename) const {
        std::ofstream file(filename);

        if (file.is_open()) {
            file << "x_c,y_c,Ri,rho,u,v,p,T\n"; // CSV header for primitive state

            for (const auto& cell : cells) {
                file << cell.x_c << "," << cell.y_c << "," // Output cell centers for visualization
                     << cell.Ri(0)/cell.volume << "," << cell.V(0) << "," << cell.V(1) << "," << cell.V(2) << "," 
                     << cell.V(3) << "," << cell.V(4) << "\n"; // Output primitive state
            }

            file.close();
            // std::cout << "\nPrimitive state data written to " << filename << std::endl;
        } else {
            std::cerr << "Unable to open file for writing: " << filename << std::endl;
        }
    }

    // Function to convert conserved state to primitive state
    void U2V(const Eigen::Matrix<long double, 4, 1>& U, Eigen::Matrix<long double, 5, 1>& V) const {
        V(QRHO) = U(URHO); // rho
        V(QU) = U(UMX) / U(URHO); // u
        V(QV) = U(UMY) / U(URHO); // v
        V(QP) = (U(UEINT) - 0.5 * (U(UMX) * U(UMX) + U(UMY) * U(UMY)) / U(URHO)) * (GAMMA - 1); // P
        V(QT) = V(QP) / R / V(QRHO); // T
    }

    // Function to convert primitive state to conserved state
    void V2U(const Eigen::Matrix<long double, 5, 1>& V, Eigen::Matrix<long double, 4, 1>& U) const {
        U(URHO) = V(0); // rho
        U(UMX) = V(0) * V(1); // rhoU
        U(UMY) = V(0) * V(2); // rhoV
        U(UEINT) = V(3) / (GAMMA - 1) + 0.5 * V(0) * (V(1) * V(1) + V(2) * V(2)); // E
    }

    // Function to output residuals to a CSV file
    void outputResidualsToCSV(const std::string& filename) const {
        std::ofstream file(filename);

        if (file.is_open()) {
            file << "Index,Residual\n"; // CSV header for residuals

            for (size_t i = 0; i < residuals.size(); ++i) {
                file << i << "," << residuals[i] << "\n"; // Output index and residual
            }

            file.close();
            std::cout << "\nResiduals data written to " << filename << std::endl;
        } else {
            std::cerr << "Unable to open file for writing: " << filename << std::endl;
        }
    }

private:
    int nx, ny, nfaces, ncells;

    std::string gridFile;
    std::vector<Node> nodes;
    std::vector<Cell> cells;
    std::vector<Face> faces;

    //initial condition stuff
    long double M_inlet = 5.0;
    long double T_inlet = 300; //K
    long double P_inlet = 10*1e3; //Pa
    
    //simulation time
    int nsteps = 550000;
    long double dx;
    long double CFL=7.0;
    long double dt;
    long double tol = 1e-8;
    int max_iters = 100;
    //residuals
    std::vector<long double>residuals;

};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1; // Exit if no input file is provided
    }

    // Coarse grid parameters
    int nx, ny;
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
                else if (key == "gridFile") csvFilename = value; // Read CSV filename
            }
        }
        inputFile.close();
    } else {
        std::cerr << "Unable to open input file: " << argv[1] << std::endl;
        return 1; // Exit if file cannot be opened
    }

    GridReader grid(nx, ny, csvFilename);
    grid.computeMinMaxFaceVectorSum();
    grid.outputToCSV(csvFilename.erase(csvFilename.find_last_of('.')) + "_cells.csv");  // Outputs grid data to CSV file

    return 0;
}
