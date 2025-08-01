#ifndef SOVLER_HPP
#define SOVLER_HPP
#include<iostream>
#include<string>
#include<vector>
#include"grid.hpp"

class solver : public grid
{
    private:
        void set_inital_condition();
        void set_boundary_condition();
        void cons_to_prim(std::vector<long double> c, std::vector<long double> &p);
        void prim_to_cons(std::vector<long double> p, std::vector<long double> &c);
        void fillDerivedVars(std::vector<std::vector<std::vector<long double>>> &p);
        void copy(std::vector<std::vector<std::vector<long double>>> p, std::vector<std::vector<std::vector<long double>>> &c);
        void advance();
        void compute_flux();
        void compute_dt_min();
        std::vector<long double> computeResidual(std::vector<std::vector<std::vector<long double>>> Un, std::vector<std::vector<std::vector<long double>>> Un1);
        long double h0(std::vector<long double> Q);
        long double c(std::vector<long double> Q);
        long double sign(long double a);
        long double findMin(long double a, long double b, long double c, long double d);
        long double minmod(long double a, long double b);
    protected:
        // to keep memory used low, we only carry the conserved states at n, n+1 separately. A single primitive state is used.
        std::vector<std::vector<std::vector<long double>>> Un, Un1, Qn, F, G, QL, QR, UL, UR, QB, QT, UB, UT;
    public:
        solver(int i,int j, std::string filename);
        void solve();
        void writeState(std::vector<std::vector<std::vector<long double>>> s, std::string fname);
        void writeWallPressure();
        int nsteps=1;
        int plt_int=10;
        int use_flux_lim=0;
        long double eps = 0; 
        long double cfl=0.8;
        long double dt=100;
        long double k=-1;
        long double beta=1;
        long double tol=1e-6;
};
#endif