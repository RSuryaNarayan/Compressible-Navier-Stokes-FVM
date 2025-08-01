#include<fstream>
void grid::writeFields()
{
    //open file
    std::ofstream outfile("plot_metrics.txt");

    //header
    outfile<<"X,"<<"Y,"<<"xi_x,"<<"xi_y,"<<"eta_x,"<<"eta_y,"<<"J\n";

    //data
    for (int i=i_min; i<=i_max;i++)
    {
        for (int j=j_min;j<=j_max;j++)
        {
            outfile<<x[i][j]<<",\t"<<y[i][j]<<",\t"<<xi_x[i][j]<<",\t"<<xi_y[i][j]<<",\t"<<eta_x[i][j]<<",\t"<<eta_y[i][j]<<",\t"<<J[i][j]<<"\n";
        }
    }

    //close
    outfile.close();
}