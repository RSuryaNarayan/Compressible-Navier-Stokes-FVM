//conserved state to primitive conversion
//The primitive state is what is primarily used for plotting 

void solver::cons_to_prim (std::vector<long double> cons,
                           std::vector<long double> &prim)
{
    prim[QRHO] = cons[URHO];
    prim[QU] = cons[UMX]/cons[URHO];
    prim[QV] = cons[UMY]/cons[URHO];
    prim[QPRES] = (GAMMA-1)*(cons[UE]-0.5*cons[URHO]*(prim[QU]*prim[QU] + prim[QV]*prim[QV]));
}