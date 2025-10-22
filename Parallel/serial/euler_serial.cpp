#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>   // for timing

using namespace std;
using namespace std::chrono;

const double gamma_ = 1.4;
enum { RHO=0, MOM=1, EN=2 };
const int NVARS = 3;

// ------------------------------------------------------------
// Convert conserved → primitive
// ------------------------------------------------------------
void cons_to_prim(const vector<double>& U, double& rho, double& u, double& p) {
    rho = U[RHO];
    u = U[MOM] / rho;
    double E = U[EN];
    p = (gamma_ - 1.0) * (E - 0.5 * rho * u * u);
}

// ------------------------------------------------------------
// Compute flux F(U)
// ------------------------------------------------------------
vector<double> flux(const vector<double>& U) {
    double rho, u, p;
    cons_to_prim(U, rho, u, p);
    double E = U[EN];
    vector<double> F(NVARS);
    F[RHO] = rho * u;
    F[MOM] = rho * u * u + p;
    F[EN]  = u * (E + p);
    return F;
}

// ------------------------------------------------------------
// HLL numerical flux
// ------------------------------------------------------------
vector<double> HLL_flux(const vector<double>& UL, const vector<double>& UR) {
    double rhoL, uL, pL;
    double rhoR, uR, pR;
    cons_to_prim(UL, rhoL, uL, pL);
    cons_to_prim(UR, rhoR, uR, pR);

    double cL = sqrt(gamma_ * pL / rhoL);
    double cR = sqrt(gamma_ * pR / rhoR);

    double SL = min(uL - cL, uR - cR);
    double SR = max(uL + cL, uR + cR);

    vector<double> FL(NVARS), FR(NVARS), FH(NVARS);
    double EL = UL[EN], ER = UR[EN];
    FL[RHO] = rhoL * uL;
    FL[MOM] = rhoL * uL * uL + pL;
    FL[EN]  = uL * (EL + pL);

    FR[RHO] = rhoR * uR;
    FR[MOM] = rhoR * uR * uR + pR;
    FR[EN]  = uR * (ER + pR);

    if (SL >= 0.0) return FL;
    if (SR <= 0.0) return FR;

    for (int k=0; k<NVARS; ++k)
        FH[k] = (SR*FL[k] - SL*FR[k] + SL*SR*(UR[k]-UL[k])) / (SR - SL);
    return FH;
}

// ------------------------------------------------------------
// Write solution to file
// ------------------------------------------------------------
void write_solution(const vector<vector<double>>& U, double x0, double dx, double t, int N) {
    stringstream ss;
    ss << "solution_t" << fixed << setprecision(2) << t << ".dat";
    ofstream out(ss.str());
    out << "# x rho u p\n";

    for (int i=0; i<N; ++i) {
        double rho, u, p;
        cons_to_prim(U[i], rho, u, p);
        double x = x0 + (i + 0.5) * dx;
        out << x << " " << rho << " " << u << " " << p << "\n";
    }
    out.close();
    cout << "  -> Wrote " << ss.str() << endl;
}

// ------------------------------------------------------------
// Main solver
// ------------------------------------------------------------
int main() {
    auto start_total = high_resolution_clock::now();

    const double x0 = 0.0, x1 = 1.0;
    const int N = 200;
    const double CFL = 0.5;
    const double t_final = 0.2;
    const double dx = (x1 - x0) / N;
    const double t_output_interval = 0.02;

    vector<vector<double>> U(N, vector<double>(NVARS, 0.0));
    vector<vector<double>> U_new(N, vector<double>(NVARS, 0.0));

    // Initial condition: Sod shock tube
    for (int i=0; i<N; ++i) {
        double x = x0 + (i + 0.5) * dx;
        double rho, u, p;
        if (x < 0.5) { rho = 1.0; u = 0.0; p = 1.0; }
        else         { rho = 0.125; u = 0.0; p = 0.1; }
        U[i][RHO] = rho;
        U[i][MOM] = rho * u;
        U[i][EN]  = p / (gamma_ - 1.0) + 0.5 * rho * u * u;
    }

    double t = 0.0, next_output = 0.0;
    int step = 0;

    write_solution(U, x0, dx, t, N);

    auto start_loop = high_resolution_clock::now();

    while (t < t_final) {
        double max_speed = 0.0;
        for (int i=0; i<N; ++i) {
            double rho, u, p;
            cons_to_prim(U[i], rho, u, p);
            double c = sqrt(gamma_ * p / rho);
            max_speed = max(max_speed, fabs(u) + c);
        }
        double dt = CFL * dx / max_speed;
        if (t + dt > t_final) dt = t_final - t;

        vector<vector<double>> fluxes(N+1, vector<double>(NVARS, 0.0));
        vector<double> UL = U[0];
        vector<double> UR = U[0];
        fluxes[0] = HLL_flux(UL, UR);

        for (int i=0; i<N-1; ++i)
            fluxes[i+1] = HLL_flux(U[i], U[i+1]);

        UL = U[N-1]; UR = U[N-1];
        fluxes[N] = HLL_flux(UL, UR);

        for (int i=0; i<N; ++i)
            for (int k=0; k<NVARS; ++k)
                U_new[i][k] = U[i][k] - (dt/dx) * (fluxes[i+1][k] - fluxes[i][k]);

        U.swap(U_new);
        t += dt; step++;

        if (t >= next_output - 1e-10) {
            write_solution(U, x0, dx, t, N);
            next_output += t_output_interval;
        }
    }

    auto end_loop = high_resolution_clock::now();
    auto end_total = high_resolution_clock::now();

    double loop_time = duration<double>(end_loop - start_loop).count();
    double total_time = duration<double>(end_total - start_total).count();

    cout << "\n================ TIMING REPORT ================\n";
    cout << "  Cells:           " << N << "\n";
    cout << "  Steps:           " << step << "\n";
    cout << "  Simulation time: " << t_final << " s\n";
    cout << "  Compute time:    " << loop_time << " s\n";
    cout << "  Total time:      " << total_time << " s\n";
    cout << "==============================================\n";

    return 0;
}
