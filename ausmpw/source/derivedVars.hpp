//stagnation enthalpy
long double solver::h0(std::vector<long double> Q)
{
    return (Q[QPRES]/Q[QRHO])*(GAMMA/(GAMMA-1)) + 0.5 * (Q[QU]*Q[QU] + Q[QV]*Q[QV]);
}

//speed of sound
long double solver::c(std::vector<long double> Q)
{
    long double c2 = GAMMA * Q[QPRES] /Q[QRHO];
    if (c2<0)
    {
        std::cout<<"p = "<<Q[QPRES]<<", rho="<<Q[QRHO]<<", c="<<std::sqrt(c2)<<"\n";
        std::cout<<"Cannot compute speed of sound!!! ABORTING!!"<<"\n";
        //exit(0);
    }

    return std::sqrt(c2); 
}

//copy Un1 to Un
void solver::copy(std::vector<std::vector<std::vector<long double>>> Un, std::vector<std::vector<std::vector<long double>>> &Un1)
{
    for (int i=i_min; i<=i_max; i++)
    {
        for (int j=j_min; j<=j_max; j++)
        {
            for (int k=0; k<NUVAR; k++)
            {
                Un1[i][j][k] = Un[i][j][k];
            }
        }
    }
}

//helper functions for computing the flux
long double solver::sign(long double a)
{
    return a==0? 0 : (a>0)? 1 : -1;
}

long double solver::findMin(long double a, long double b, long double c, long double d)
{
    long double min1 = (a < b) ? a : b; // Find minimum between a and b
    long double min2 = (c < d) ? c : d; // Find minimum between c and d

    return (min1 < min2) ? min1 : min2; // Find minimum between min1 and min2
}

//minmod function
long double solver::minmod(long double x, long double y)
{
    long double min1 = std::min(x * sign(y), y * sign(x));
    long double max1  = 0>=min1? 0 : min1;
    return sign(x) * max1;
}

//fill derived variables into the primitive state
void solver::fillDerivedVars(std::vector<std::vector<std::vector<long double>>> &prim)
{

    for (int i=i_min; i<=i_max; i++)
    {
        for (int j=j_min; j<=j_max; j++)
        {
            prim[i][j][QC] = 1.0*c(prim[i][j]);
            prim[i][j][QHO] = h0(prim[i][j]);
            prim[i][j][QUT] = xi_x[i][j] * prim[i][j][QU] + xi_y[i][j] * prim[i][j][QV];
            prim[i][j][QVT] = eta_x[i][j] * prim[i][j][QU] + eta_y[i][j] * prim[i][j][QV];
            prim[i][j][QM] = std::sqrt(prim[i][j][QU]*prim[i][j][QU] + prim[i][j][QV]*prim[i][j][QV]) /prim[i][j][QC];
        }
    }
}
