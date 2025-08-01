std::vector<long double> solver::computeResidual(std::vector<std::vector<std::vector<long double>>> Un,
                                                 std::vector<std::vector<std::vector<long double>>> Un1)
{
    std::vector<long double> res(4, 0.0);

    for (int i=i_min; i<=i_max; i++)
    {
        for (int j=j_min; j<=j_max; j++)
        {
            for (int k=0; k<NUVAR; k++)
            {
                res[k] = res[k] + (Un1[i][j][k] - Un[i][j][k])*(Un1[i][j][k] - Un[i][j][k]);
            }
        }
    }

    for (int k=0; k<NUVAR; k++)
    {
        res[k] = std::sqrt(res[k]/i_max/j_max);
    }

    return res;
}