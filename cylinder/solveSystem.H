// Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> BigA;
// Eigen::Matrix<long double, Eigen::Dynamic, 1> BigB;

// Initialize the sizes at runtime
// BigA = Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>::Zero(4 * ncells, 4 * ncells);
// BigB = Eigen::Matrix<long double, Eigen::Dynamic, 1>::Zero(4 * ncells, 1);

//fill up the matrix entries
for (int iic=0; iic<ncells; iic++)
{
    //fill diagonal terms
    Cell cell = cells[iic];
    BigA.block(4*iic,4*iic,4,4) = cell.Di;

    //fill off-diagonal terms
    for (size_t inbr=0; inbr<cell.neighbors.size(); inbr++)
    {
        BigA.block(4*iic,4*cell.neighbors[inbr],4,4) = cell.N[inbr];
    }

    //fill residuals
    BigB.block(4*iic,0,4,1) = cell.Ri;
}

// Solve the linear system BigA * X = BigB
// Eigen::Matrix<long double, Eigen::Dynamic, 1> dU;
dU = BigA.colPivHouseholderQr().solve(BigB);