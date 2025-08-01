void solver::advance()
{
    for (int i=i_min+1; i<=i_max-1; i++)
    {
        for (int j=j_min+1; j<=j_max-1; j++)
        {
            for (int k=0; k<NUVAR; k++)
            {
                Un[i][j][k] = Un[i][j][k] 
                            - dt * J[i][j] * 
                            ( (F[i][j][k]-F[i-1][j][k])
                         +    (G[i][j][k]-G[i][j-1][k]) );
            }

            //convert to primitive immediately
            cons_to_prim(Un[i][j], Qn[i][j]);
        }
    }
}