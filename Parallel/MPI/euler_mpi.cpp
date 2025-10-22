#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <map>
#include <stdexcept> // For runtime_error
#include <mpi.h>

using namespace std;

const double gamma_ = 1.4;
enum { RHO=0, MOM=1, EN=2 };
const int NVARS = 3;

// --- A struct to hold all simulation parameters ---
struct SimParams {
    double x0, x1, CFL, t_final, t_output_interval;
    int N_total;
};

// ==========================================================
// NEW FUNCTION: Rank 0 reads the input file
// ==========================================================
SimParams read_input_file(const string& filename) {
    SimParams params;
    map<string, bool> found_params = {
        {"x0", false}, {"x1", false}, {"N_total", false},
        {"CFL", false}, {"t_final", false}, {"t_output_interval", false}
    };

    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Could not open input file: " + filename);
    }

    string line;
    while (getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') {
            continue;
        }

        stringstream ss(line);
        string key;
        ss >> key;

        if (key == "x0")                  { ss >> params.x0; found_params["x0"] = true; }
        else if (key == "x1")             { ss >> params.x1; found_params["x1"] = true; }
        else if (key == "N_total")        { ss >> params.N_total; found_params["N_total"] = true; }
        else if (key == "CFL")            { ss >> params.CFL; found_params["CFL"] = true; }
        else if (key == "t_final")        { ss >> params.t_final; found_params["t_final"] = true; }
        else if (key == "t_output_interval") { ss >> params.t_output_interval; found_params["t_output_interval"] = true; }
        else {
            cerr << "Warning: Unknown key in input file: " << key << endl;
        }
    }
    file.close();

    // Check if all parameters were found
    for (auto const& [key, val] : found_params) {
        if (!val) {
            throw runtime_error("Missing parameter in input file: " + key);
        }
    }
    return params;
}


void cons_to_prim(const vector<double>& U, double& rho, double& u, double& p) {
    rho = U[RHO]; u = U[MOM]/rho; double E=U[EN]; p=(gamma_-1.0)*(E-0.5*rho*u*u);
}

vector<double> HLL_flux(const vector<double>& UL, const vector<double>& UR) {
    double rhoL,uL,pL,rhoR,uR,pR; cons_to_prim(UL,rhoL,uL,pL); cons_to_prim(UR,rhoR,uR,pR);
    double cL=sqrt(gamma_*pL/rhoL), cR=sqrt(gamma_*pR/rhoR);
    double SL=min(uL-cL,uR-cR), SR=max(uL+cL,uR+cR);
    vector<double> FL(NVARS),FR(NVARS),FH(NVARS);
    double EL=UL[EN], ER=UR[EN];
    FL[RHO]=rhoL*uL; FL[MOM]=rhoL*uL*uL+pL; FL[EN]=uL*(EL+pL);
    FR[RHO]=rhoR*uR; FR[MOM]=rhoR*uR*uR+pR; FR[EN]=uR*(ER+pR);
    if(SL>=0) return FL; if(SR<=0) return FR;
    for(int k=0;k<NVARS;k++) FH[k]=(SR*FL[k]-SL*FR[k]+SL*SR*(UR[k]-UL[k]))/(SR-SL);
    return FH;
}

// ------------------------------------------------------------
void write_global_solution(const vector<double>& global_buffer, double x0, double dx,
                           double t, int N_total) {
    ofstream out;
    stringstream ss;
    ss << "solution_t" << fixed << setprecision(2) << t << ".dat";
    out.open(ss.str());
    out << "# x rho u p\n";
    for(int i=0;i<N_total;i++){
        double rho = global_buffer[i*NVARS+RHO];
        double mom = global_buffer[i*NVARS+MOM];
        double E   = global_buffer[i*NVARS+EN];
        double u = mom/rho;
        double p = (gamma_-1.0)*(E-0.5*rho*u*u);
        double x = x0 + (i + 0.5) * dx;
        out << x << " " << rho << " " << u << " " << p << "\n";
    }
    out.close();
    cout << "  -> wrote " << ss.str() << endl;
}

// ------------------------------------------------------------
int main(int argc, char** argv){
    MPI_Init(&argc,&argv);
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    // ==========================================================
    // READ AND BROADCAST PARAMETERS
    // ==========================================================
    if (argc != 2) {
        if (rank == 0) {
            cerr << "Usage: " << argv[0] << " <input_file.txt>" << endl;
        }
        MPI_Finalize();
        return 1;
    }
    string input_filename = argv[1];
    SimParams params;

    if (rank == 0) {
        try {
            params = read_input_file(input_filename);
            cout << "--- Simulation Parameters ---" << endl;
            cout << "  Domain:    [" << params.x0 << ", " << params.x1 << "]" << endl;
            cout << "  Cells:     " << params.N_total << endl;
            cout << "  T_final:   " << params.t_final << endl;
            cout << "  CFL:       " << params.CFL << endl;
            cout << "  Output Int: " << params.t_output_interval << endl;
            cout << "-----------------------------" << endl;
        } catch (const std::exception& e) {
            cerr << "Error: " << e.what() << endl;
            MPI_Abort(MPI_COMM_WORLD, 1); // Abort all processes
        }
    }

    // Broadcast the entire params struct from rank 0 to all other ranks
    // We send it as a block of bytes
    MPI_Bcast(&params, sizeof(SimParams), MPI_BYTE, 0, MPI_COMM_WORLD);

    // --- Replace hard-coded constants with values from params ---
    const double x0 = params.x0;
    const double x1 = params.x1;
    const int N_total = params.N_total;
    const double CFL = params.CFL;
    const double t_final = params.t_final;
    const double t_output_interval = params.t_output_interval;
    // -----------------------------------------------------------

    const double dx=(x1-x0)/N_total;

    // Compute local sizes (non-divisible N_total)
    int base = N_total / size;
    int remainder = N_total % size;
    int N_local = base + (rank < remainder ? 1 : 0);
    int offset = rank*base + min(rank,remainder);

    vector<vector<double>> U(N_local+2, vector<double>(NVARS,0.0));
    vector<vector<double>> U_new(N_local+2, vector<double>(NVARS,0.0));

    // ------------------------------
    // Initialize local domain
    // ------------------------------
    for(int i=1;i<=N_local;i++){
        double x = x0 + (offset + i-1 + 0.5)*dx;
        double rho,u,p;
        if(x<0.5 * (x0 + x1)){ // Use midpoint relative to domain
             rho=1.0; u=0.0; p=1.0; 
        } else {
             rho=0.125; u=0.0; p=0.1; 
        }
        U[i][RHO]=rho; U[i][MOM]=rho*u; U[i][EN]=p/(gamma_-1.0)+0.5*rho*u*u;
    }

    // ------------------------------
    // Initialize ghost cells (all ranks)
    // ------------------------------
    if(rank == 0)        U[0] = U[1];            
    if(rank == size-1)   U[N_local+1] = U[N_local]; 
    if(rank > 0)         U[0] = U[1];
    if(rank < size-1)    U[N_local+1] = U[N_local];

    double t=0.0, next_output=0.0;
    int step=0;
    double t_total_start = MPI_Wtime();
    double t_compute = 0.0, t_output = 0.0;

    // ------------------------------
    // Write initial condition (t=0) before timestep loop
    // ------------------------------
    if (t_output_interval > 0) { // Only write if output is requested
        vector<double> local_buffer(N_local*NVARS);
        for(int i=0;i<N_local;i++)
            for(int k=0;k<NVARS;k++)
                local_buffer[i*NVARS+k]=U[i+1][k];

        vector<int> recvcounts(size), displs(size);
        if(rank==0){
            for(int r=0;r<size;r++){
                int nl = base + (r<remainder?1:0);
                recvcounts[r]=nl*NVARS;
                displs[r]=(r*base+min(r,remainder))*NVARS;
            }
        }

        vector<double> global_buffer;
        if(rank==0) global_buffer.resize(N_total*NVARS);

        MPI_Gatherv(local_buffer.data(), N_local*NVARS, MPI_DOUBLE,
                    global_buffer.data(), recvcounts.data(), displs.data(), MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        if(rank==0) write_global_solution(global_buffer,x0,dx,0.0,N_total);
        next_output += t_output_interval;
    }

    // ------------------------------
    // Main timestep loop
    // ------------------------------
    while(t<t_final){

        // ------------------------------
        // Exchange ghost cells
        // ------------------------------
        int left=(rank>0)?rank-1:MPI_PROC_NULL;
        int right=(rank<size-1)?rank+1:MPI_PROC_NULL;

        MPI_Sendrecv(&U[N_local][0], NVARS, MPI_DOUBLE, right, 0,
                     &U[0][0], NVARS, MPI_DOUBLE, left, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&U[1][0], NVARS, MPI_DOUBLE, left, 1,
                     &U[N_local+1][0], NVARS, MPI_DOUBLE, right, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // ------------------------------
        // Compute fluxes / update
        // ------------------------------
        auto t_comp_start = MPI_Wtime();

        double local_max=0.0;
        for(int i=1;i<=N_local;i++){
            double rho,u,p; cons_to_prim(U[i],rho,u,p);
            local_max=max(local_max,fabs(u)+sqrt(gamma_*p/rho));
        }
        double global_max;
        MPI_Allreduce(&local_max,&global_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        double dt = CFL*dx/global_max;
        if(t+dt>t_final) dt=t_final-t;
        if(t+dt <= t) {
            if(rank==0) cerr << "Warning: dt is too small, exiting." << endl;
            break; 
        }

        vector<vector<double>> fluxes(N_local+1, vector<double>(NVARS,0.0));
        for(int i=0;i<=N_local;i++) fluxes[i]=HLL_flux(U[i],U[i+1]);

        for(int i=1;i<=N_local;i++)
            for(int k=0;k<NVARS;k++)
                U_new[i][k]=U[i][k]-(dt/dx)*(fluxes[i][k]-fluxes[i-1][k]);
        
        U.swap(U_new);
        
        // FIX: Re-apply physical boundary conditions
        if(rank == 0)        U[0] = U[1];
        if(rank == size-1)   U[N_local+1] = U[N_local];

        auto t_comp_end = MPI_Wtime();
        t_compute += t_comp_end - t_comp_start;

        t += dt; step++;

        // ------------------------------
        // Gather & write output AFTER ghost exchange and update
        // ------------------------------
        
        // Only do output if interval is positive and time is reached
        if (t_output_interval > 0 && t >= next_output - 1e-10) {
            auto t_out_start = MPI_Wtime();

            vector<double> local_buffer(N_local*NVARS);
            for(int i=0;i<N_local;i++)
                for(int k=0;k<NVARS;k++)
                    local_buffer[i*NVARS+k]=U[i+1][k];

            vector<int> recvcounts(size), displs(size);
            if(rank==0){
                for(int r=0;r<size;r++){
                    int nl = base + (r<remainder?1:0);
                    recvcounts[r]=nl*NVARS;
                    displs[r]=(r*base+min(r,remainder))*NVARS;
                }
            }

            vector<double> global_buffer;
            if(rank==0) global_buffer.resize(N_total*NVARS);

            MPI_Gatherv(local_buffer.data(), N_local*NVARS, MPI_DOUBLE,
                        global_buffer.data(), recvcounts.data(), displs.data(), MPI_DOUBLE,
                        0, MPI_COMM_WORLD);

            if(rank==0) write_global_solution(global_buffer,x0,dx,t,N_total);

            auto t_out_end = MPI_Wtime();
            t_output += t_out_end - t_out_start;

            next_output += t_output_interval;
        }
    }
    
    // Add a final output step if t_final was reached and not just written
    if (t_output_interval > 0 && abs(t - t_final) < 1e-10 && abs(t - (next_output - t_output_interval)) > 1e-10) {
        auto t_out_start = MPI_Wtime();
        vector<double> local_buffer(N_local*NVARS);
        for(int i=0;i<N_local;i++)
            for(int k=0;k<NVARS;k++)
                local_buffer[i*NVARS+k]=U[i+1][k];
        vector<int> recvcounts(size), displs(size);
        if(rank==0){
            for(int r=0;r<size;r++){
                int nl = base + (r<remainder?1:0);
                recvcounts[r]=nl*NVARS;
                displs[r]=(r*base+min(r,remainder))*NVARS;
            }
        }
        vector<double> global_buffer;
        if(rank==0) global_buffer.resize(N_total*NVARS);
        MPI_Gatherv(local_buffer.data(), N_local*NVARS, MPI_DOUBLE,
                    global_buffer.data(), recvcounts.data(), displs.data(), MPI_DOUBLE,
                    0, MPI_COMM_WORLD);
        if(rank==0) write_global_solution(global_buffer,x0,dx,t,N_total);
        auto t_out_end = MPI_Wtime();
        t_output += t_out_end - t_out_start;
    }


    double t_total_end = MPI_Wtime();
    double total_time = t_total_end - t_total_start;

    double max_compute, max_output, max_total;
    MPI_Reduce(&t_compute, &max_compute, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_output, &max_output, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_time, &max_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(rank==0){
        cout << "\n================ TIMING REPORT ================\n";
        cout << "  Cells:           " << N_total << "\n";
        cout << "  Ranks:           " << size << "\n";
        cout << "  Steps:           " << step << "\n";
        cout << "  Final time:      " << t << " s\n";
        cout << "  Compute time:    " << max_compute << " s\n";
        cout << "  Output time:     " << max_output << " s\n";
        cout << "  Total time:      " << max_total << " s\n";
        cout << "==============================================\n";
    }

    MPI_Finalize();
    return 0;
}
