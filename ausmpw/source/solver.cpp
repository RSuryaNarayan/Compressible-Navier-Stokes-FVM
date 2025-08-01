/* Constructor definition for solver class*/
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include<algorithm>
#include"indexDefines.hpp"

#include"solver.hpp"
#include"setInitialCondition.hpp"
#include"consToPrim.hpp"
#include"primToCons.hpp"
#include"derivedVars.hpp"
#include"solve.hpp"
#include"computeFlux.hpp"
#include"advance.hpp"
#include"writeState.hpp"
#include"setBoundaryCondition.hpp"
#include"computeDtMin.hpp"
#include"computeResidual.hpp"
#include"writeWallPressure.hpp"

solver::solver(int i, int j, std::string fname) : grid(i, j, fname)
{
    Un.resize(i+1, std::vector<std::vector<long double> >(j+1, std::vector<long double>(NUVAR)));
    Un1.resize(i+1, std::vector<std::vector<long double> >(j+1, std::vector<long double>(NUVAR)));
    UL.resize(i+1, std::vector<std::vector<long double> >(j+1, std::vector<long double>(NUVAR)));
    UR.resize(i+1, std::vector<std::vector<long double> >(j+1, std::vector<long double>(NUVAR)));
    Qn.resize(i+1, std::vector<std::vector<long double> >(j+1, std::vector<long double>(NQVAR)));
    QL.resize(i+1, std::vector<std::vector<long double> >(j+1, std::vector<long double>(NQVAR)));
    QR.resize(i+1, std::vector<std::vector<long double> >(j+1, std::vector<long double>(NQVAR)));
    UB.resize(i+1, std::vector<std::vector<long double> >(j+1, std::vector<long double>(NUVAR)));
    UT.resize(i+1, std::vector<std::vector<long double> >(j+1, std::vector<long double>(NUVAR)));
    QB.resize(i+1, std::vector<std::vector<long double> >(j+1, std::vector<long double>(NQVAR)));    
    QT.resize(i+1, std::vector<std::vector<long double> >(j+1, std::vector<long double>(NQVAR)));
    F.resize(i+1, std::vector<std::vector<long double> >(j+1, std::vector<long double>(NQVAR)));    
    G.resize(i+1, std::vector<std::vector<long double> >(j+1, std::vector<long double>(NQVAR)));
    //Qn1.resize(i+1, std::vector<std::vector<long double> >(j+1, std::vector<long double>(NQVAR)));
};