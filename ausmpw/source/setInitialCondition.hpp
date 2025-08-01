//sets the initial condition
void solver::set_inital_condition()
{
    for(int i=i_min; i<=i_max; i++)
    {
        for (int j=j_min; j<=j_max; j++)
        {
            Qn[i][j][QRHO] = 4.0;
            Qn[i][j][QU] = 1.0;
            Qn[i][j][QV] = 0.0;
            Qn[i][j][QPRES] = 1/GAMMA;

            //fill in the conserved state
            prim_to_cons(Qn[i][j], Un[i][j]);
        }
    }
}