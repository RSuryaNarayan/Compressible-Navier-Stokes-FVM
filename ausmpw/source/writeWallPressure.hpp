void solver::writeWallPressure()
{
    std::ofstream wallOut("wallPressure.txt");

    for (int i=i_min;i<=i_max; i++)
    {
        wallOut<<x[i][j_min]<<","<<Qn[i][j_min][QPRES]<<"\n";
    }
    wallOut.close();
}