//imposes inflow BC
void solver::set_boundary_condition()
{
    //inflow face i.e. i=i_min
    for (int j=j_min; j<=j_max; j++)
    {
        Qn[i_min][j][QRHO]  = 4.0;
        Qn[i_min][j][QU]    = 1.0;
        Qn[i_min][j][QV]    = 0.0;
        Qn[i_min][j][QPRES] = 1/GAMMA;

        //transfer info to conserved state
        prim_to_cons(Qn[i_min][j], Un[i_min][j]);
    }
    
    // outflow face i.e. i=i_max
    for (int j=j_min; j<=j_max; j++)
    {
        for (int k=0; k<NUVAR; k++)
        {
            // fill in conserved state
            Un[i_max][j][k]=Un[i_max-1][j][k];

            // fill in primitive state
            cons_to_prim(Un[i_max][j], Qn[i_max][j]);
        }
    }

    //bottom wall i.e. j=j_min
    for (int i=i_min; i<=i_max; i++)
    {

        long double s1 = std::sqrt( xi_x[i][j_min]   * xi_x[i][j_min]   + xi_y[i][j_min]   * xi_y[i][j_min]   );
        long double s2 = std::sqrt( xi_x[i][j_min+1] * xi_x[i][j_min+1] + xi_y[i][j_min+1] * xi_y[i][j_min+1] );
        long double s3 = std::sqrt( xi_x[i][j_min+2] * xi_x[i][j_min+2] + xi_y[i][j_min+2] * xi_y[i][j_min+2] );

        long double delta2 = s1/s2;
        long double delta3 = s1/s3;

        long double u_t2 = xi_x[i][j_min+1]*Qn[i][j_min+1][QU] + xi_y[i][j_min+1]*Qn[i][j_min+1][QV];
        long double u_t3 = xi_x[i][j_min+2]*Qn[i][j_min+2][QU] + xi_y[i][j_min+2]*Qn[i][j_min+2][QV];

        long double vtemp = (2*u_t2*delta2 - 1*u_t3*delta3);

        long double u1   = (eta_y[i][j_min]/J[i][j_min])*vtemp;
        long double v1   = -1*(eta_x[i][j_min]/J[i][j_min])*vtemp;
        long double p1   = Qn[i][j_min+1][QPRES];
        long double h02  = h0(Qn[i][j_min+1]);
        long double rho1 = p1*(GAMMA/(GAMMA-1))/(h02-0.5*(u1*u1+v1*v1));

        // std::cout<<"u1 = "<<u1<<"\n";
        // std::cout<<"v1 = "<<v1<<"\n";
        // std::cout<<"p1 = "<<p1<<"\n";
        // std::cout<<"h02 = "<<h02<<"\n";
        // std::cout<<"rho1 = "<<rho1<<"\n";
        // std::cout<<"h02-0.5*(u1*u1+v1*v1) = "<<h02-0.5*(u1*u1+v1*v1)<<"\n";
        // std::cout<<"---------"<<"\n";

        //fill in primitive
        Qn[i][j_min][QRHO]  = rho1;
        Qn[i][j_min][QU]    = u1;
        Qn[i][j_min][QV]    = v1;
        Qn[i][j_min][QPRES] = p1;

        //zero-gradient at the bottom wall
        // Qn[i][j_min][QRHO]  = Qn[i][j_min+1][QRHO];
        // Qn[i][j_min][QU]    = Qn[i][j_min+1][QU];
        // Qn[i][j_min][QV]    = Qn[i][j_min+1][QV];
        // Qn[i][j_min][QPRES] = Qn[i][j_min+1][QPRES];
        
        //fill in cons
        prim_to_cons(Qn[i][j_min],Un[i][j_min]);
    }

    //top wall i.e. j=j_max
    for (int i=i_min; i<=i_max; i++)
    {

        long double s1 = std::sqrt(xi_x[i][j_max]   * xi_x[i][j_max]   + xi_y[i][j_max]   * xi_y[i][j_max]  );
        long double s2 = std::sqrt(xi_x[i][j_max-1] * xi_x[i][j_max-1] + xi_y[i][j_max-1] * xi_y[i][j_max-1]);
        long double s3 = std::sqrt(xi_x[i][j_max-2] * xi_x[i][j_max-2] + xi_y[i][j_max-2] * xi_y[i][j_max-2]);

        long double delta2 = s1/s2;
        long double delta3 = s1/s3;

        long double u_t1 = xi_x[i][j_max-1]*Qn[i][j_max-1][QU] + xi_y[i][j_max-1]*Qn[i][j_max-1][QV];
        long double u_t2 = xi_x[i][j_max-2]*Qn[i][j_max-2][QU] + xi_y[i][j_max-2]*Qn[i][j_max-2][QV];

        long double vtemp = (2*u_t1*delta2 - 1*u_t2*delta3);

        long double u1   = (eta_y[i][j_max]/J[i][j_max])*vtemp;
        long double v1   = -1*(eta_x[i][j_max]/J[i][j_max])*vtemp;
        long double p1   = Qn[i][j_max-1][QPRES];
        long double h02  = h0(Qn[i][j_max-1]);
        long double rho1 = p1*(GAMMA/(GAMMA-1))/(h02-0.5*(u1*u1+v1*v1));

        //fill in primitive
        Qn[i][j_max][QRHO]  = rho1;
        Qn[i][j_max][QU]    = u1;
        Qn[i][j_max][QV]    = v1;
        Qn[i][j_max][QPRES] = p1;

        //zero-gradient at the top wall
        // Qn[i][j_max][QRHO]  = Qn[i][j_max-1][QRHO];
        // Qn[i][j_min][QU]    = Qn[i][j_max-1][QU];
        // Qn[i][j_min][QV]    = Qn[i][j_max-1][QV];
        // Qn[i][j_min][QPRES] = Qn[i][j_max-1][QPRES];

        //fill in cons
        prim_to_cons(Qn[i][j_max],Un[i][j_max]);
    }
}