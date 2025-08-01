//primitive state to conservative conversion

void solver::prim_to_cons(std::vector<long double> prim,
                          std::vector<long double> &cons)
{
   
    cons[URHO] = prim[QRHO];
    cons[UMX]= prim[QU] * prim[QRHO];
    cons[UMY]= prim[QV] * prim[QRHO];
    cons[UE]= prim[QPRES]/(GAMMA-1) + 0.5 * prim[QRHO] * (prim[QU]*prim[QU] + prim[QV]*prim[QV]);
}