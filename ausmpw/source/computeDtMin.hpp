//computes the minimum dt required given a primitive state
void solver::compute_dt_min()
{
    long double dt_min = 100; //set large value initially

    for (int i=i_min; i<=i_max; i++)
    {
        for (int j=j_min; j<=j_max; j++)
        {
            long double u_tilde = xi_x[i][j]*Qn[i][j][QU]+xi_y[i][j]*Qn[i][j][QV];
            long double v_tilde = eta_x[i][j]*Qn[i][j][QU]+eta_y[i][j]*Qn[i][j][QV];
            long double uc = c(Qn[i][j]);
            long double t1 = std::abs(u_tilde);
            long double t2 = std::abs(v_tilde);
            long double t3 = uc
                        *std::sqrt(
                        xi_x[i][j]*xi_x[i][j] + xi_y[i][j]*xi_y[i][j] 
                        + eta_x[i][j]*eta_x[i][j] + eta_y[i][j]*eta_y[i][j]
                        + 2*std::abs(xi_x[i][j]*eta_x[i][j] + xi_y[i][j]*eta_y[i][j])
                        );
            long double dt_new = cfl/(t1 + t2 + t3);

            if (dt_new<dt_min)
            {
                dt_min = dt_new;
            }
        }
    }

    if (dt_min<0)
    {
        std::cout<<"dt_min<0!!! ABORTING!!!\n";
        exit(0);
    }

    std::cout<<"Computed dt_min = "<<dt_min<<"\n";
    dt=dt_min;
}