//subroutine that outputs the state at any time instant to a text file for visualization
//Note: takes in only the *primitive state*
void solver::writeState(std::vector<std::vector<std::vector<long double>>> state, std::string fname)
{
    //get derived variables
    fillDerivedVars(state);

    //open file
    std::ofstream fout(fname);

    //write header
    fout<<"X,"<<"Y,"<<"RHO,"<<"U,"<<"V,"<<"P,"<<"C,"<<"H0,"<<"U_tilda,"<<"V_tilda,"<<"M"<<"\n";

    //write data
    for (int i=i_min; i<=i_max; i++)
    {
        for (int j=j_min; j<=j_max; j++)
        {
            fout<<x[i][j]<<","<<y[i][j]<<","<<state[i][j][QRHO]
                                         <<","<<state[i][j][QU]
                                         <<","<<state[i][j][QV]
                                         <<","<<state[i][j][QPRES]
                                         <<","<<state[i][j][QC]
                                         <<","<<state[i][j][QHO]
                                         <<","<<state[i][j][QUT]
                                         <<","<<state[i][j][QVT]
                                         <<","<<state[i][j][QM]<<",\n";
        }
    }

    //close
    fout.close();
}