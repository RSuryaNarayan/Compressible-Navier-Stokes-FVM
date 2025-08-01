#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include<chrono>

#include"grid.hpp"
#include"solver.hpp"
#include"indexDefines.hpp"

using namespace std;


int main(int argc, char *argv[])
{
    //measure time of execution
    auto start = std::chrono::high_resolution_clock::now();
    
    //declare an object with dummy inputs first
    solver s(0,0,"dummy.txt");
    int choice;
    std::string fname;

    //we will now read input from the command line as a text file
    std::string line;
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " input_file" << std::endl;
        return 1;
    }

    std::ifstream inputFile(argv[1]);
    if (!inputFile.is_open()) {
        std::cerr << "Unable to open input file: " << argv[1] << std::endl;
        return 1;
    }

     while (std::getline(inputFile, line)) {
        // Ignore lines starting with '#'
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        std::string key, value;
        
        if (std::getline(iss, key, '=') && std::getline(iss, value, '#')) {
            // Trim leading and trailing spaces from key and value
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);

            // Assign values to the Solver object based on key
            if (key == "choice") {
                choice = std::stod(value);
                //re-declare the solver object with the right grids based on choice
                if (choice==0)
                {
                    std::cout<<"Solving on a coarse grid (33x21)\n";
                    fname = "../grids/coarse.txt";
                    // std::string fname_test = "../grids/test_grid.txt";
                    s = solver(33,21, fname); 
                }
                else
                {
                    std::cout<<"Solving on a fine grid (71x48)\n";
                    fname = "../grids/fine.txt";
                    s = solver(71,48, fname); 
                }
            } else if (key == "cfl") {
                s.cfl = std::stod(value);
            } else if (key == "nsteps") {
                s.nsteps = std::stoi(value);
            } else if (key == "plt_int") {
                s.plt_int = std::stoi(value);
            } else if (key == "eps") {
                s.eps = std::stod(value);
            } else if (key == "beta") {
                s.beta = std::stod(value);    
            } else if (key == "k") {
                s.k = std::stod(value);
            } else if (key == "use_flux_lim") {
                s.use_flux_lim = std::stoi(value);
            } else if (key == "tol") {
                s.tol = std::stod(value);
            }
        }
    }

    // Close the file
    inputFile.close();

    std::cout<<"\nBeginning solve with the following parameters:\n";
    std::cout<<"CFL= "<<s.cfl<<"\n";
    std::cout<<"Max steps= "<<s.nsteps<<"\n";
    std::cout<<"Plotting interval= "<<s.plt_int<<"\n";
    std::cout<<"Epsilon (0-->first order accurate | 1-->Second order accurate | no flux limiter)= "<<s.eps<<"\n";
    std::cout<<"Using flux limiter? (0-->no | 1-->yes)= "<<s.use_flux_lim<<"\n";
    std::cout<<"k value used in flux limiter= "<<s.k<<"\n";
    std::cout<<"Tolerance criterion for deciding steady-state/convergence = "<<s.tol<<"\n";

    //solve!
    s.solve();
    
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(stop - start);
    std::cout<<"\nExecution time = "<<duration.count()/1000000<<" seconds\n";
    //end
    return 0;
}