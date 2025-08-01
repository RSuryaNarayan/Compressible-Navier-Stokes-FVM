//Solve the equation A(dU) = B
BigA.reserve(4 * ncells * 4 * (1 + 4)); 
for (auto& cell : cells)
{

    //fill diagonal terms
    for (int i=0; i<4; i++)
    {
        BigA.insert(4*cell.index+i, 4*cell.index+i) = cell.Di(i, i);
    }

    //fill off-diagonal terms
    for (size_t inbr=0; inbr<cell.neighbors.size(); inbr++)
    {
        int nbr_idx = cell.neighbors[inbr];
        // Use a single insert for off-diagonal terms
        for (int i=0; i<4; i++)
        {
            for (int j=0; j<4; j++)
            {
                BigA.insert(4*cell.index+i, 4*nbr_idx+j) = cell.N[inbr](i, j);
            }
        }
    }
    // Fill residual vector
    for (int i = 0; i < 4; ++i) {
        BigB(4 * cell.index + i) = cell.Ri(i);
    }
}

// Solve the linear system BigA * X = BigB
BigA.makeCompressed();  // Convert to compressed format

// Solve the system using BiCGSTAB solver
Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
solver.compute(BigA);

if (solver.info() != Eigen::Success) {
    std::cerr << "Solver setup failed!" << std::endl;
    // return -1;
}
std::cout<<"\nsolving...";
dU = solver.solve(BigB);

if (solver.info() != Eigen::Success) {
    std::cerr << "Solving failed!" << std::endl;
    exit(1);
}

//clear out the matrices so they can be reused next time-step
BigA.setZero();
BigB.setZero();

