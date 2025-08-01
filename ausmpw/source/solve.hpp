//the main driver of the solver class
void solver::solve()
{
    //set initial conditions
    set_inital_condition();
    writeState(Qn, "plt_initial.txt");

    std::ofstream resout("residuals.txt");

    //outer time loop
    for (int n=1; n<=nsteps; n++)
    {
        std::cout<<"\nComputing step "<<n<<"\n";

        //copy the full state at n to a new container
        copy(Un, Un1);

        //compute fluxes
        compute_flux(); 

        //compute dt_min
        compute_dt_min();

        //advance by one step the interor points
        advance();

        // update boundaries
        set_boundary_condition();

        // write output
        if (n%plt_int==0)
        {
            writeState(Qn, "plt_"+std::to_string(n)+".txt");
        }
        
        //get residuals for the step
        auto res = computeResidual(Un, Un1);
        resout<<n<<","<<res[0]<<","<<res[1]<<","<<res[2]<<","<<res[3]<<"\n";

        if (res[0]<tol && res[1]<tol && res[2]<tol && res[3]<tol)
        {
            std::cout<<"\nConverged after "<<n<<" iterations!\n";
            writeState(Qn, "plt_converged.txt");
            break;
        }
    }

    resout.close();

    //write out wall pressure
    writeWallPressure();
}