#include<fstream>
#include <sstream>

void grid::initGrid(std::string fname)
{
    std::ifstream infile(fname);
    std::string line;

    if (!infile.is_open()) {
        std::cerr << "Error opening the file!" << std::endl;
        return;
    }

    //read values into the file
    long double xcoord, ycoord;
    while ( infile>>xcoord && infile>>ycoord) {
        x1D.push_back(xcoord);
        y1D.push_back(ycoord);
    }
    infile.close();

    //re-shape them into two-dimensional arrays of size i_max x j_max;
    int k=0;
    for (int i=i_min; i<=i_max; i++)
    {
        for (int j=j_min; j<=j_max; j++)
        {
            x[i][j] = x1D[k];
            y[i][j] = y1D[k++];
        }
    }

    //get grid metics
    computeMetrics();
}