for (auto& cell : cells)
{
    // Cell& cell = cells[iic]; // Use reference to modify the cell directly
    long double rho =  P_inlet/R/T_inlet;
    long double u = M_inlet * sqrt(GAMMA*R*T_inlet);
    cell.U(URHO) = rho; 
    cell.U(UMX) = rho*u;
    cell.U(UMY) = 0.0L;
    cell.U(UEINT) = P_inlet/(GAMMA-1) + 0.5 * rho * (u*u);

    //delU is zero initially
    cell.dU_old(URHO) = 0.0; 
    cell.dU_old(UMX) = 0.0;
    cell.dU_old(UMY) = 0.0;
    cell.dU_old(UEINT) = 0.0;

    //delU is zero initially
    cell.dU_new(URHO) = 0.0; 
    cell.dU_new(UMX) = 0.0;
    cell.dU_new(UMY) = 0.0;
    cell.dU_new(UEINT) = 0.0;
    
    // Populate primitive state using U2V function
    Eigen::Matrix<long double, 5, 1> V; // Primitive state vector
    U2V(cell.U, V); // Convert conserved state to primitive state
    cell.V = V; // Store the primitive state in the cell
}
//output the initial condition to check
outputPrimitiveStateToCSV("solution0.csv"); // Output initial conditions to CSV