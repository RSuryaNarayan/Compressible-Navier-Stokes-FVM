#include<iostream>
#include<fstream>

int main()
{
    std::ofstream fout("test_grid.txt");

    for (int i=1; i<=33; i++)
    {
        for (int j=1;j<=21; j++)
        {
            long double x = -0.5 + long double(i-1)*(2.0)/32.0;
            long double y = long double(j-1)*2.0/20.0;
            fout<<x<<"\t"<<y<<"\n";
        }
    }
    fout.close();
    return 0;
}