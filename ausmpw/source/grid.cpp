/* Constructor definition for grid class*/
#include<string>
#include"grid.hpp"
#include"initGrid.hpp"
#include"computeMetrics.hpp"
#include"writeField.hpp"


grid::grid(int i, int j, std::string fname)
{
    i_min = 1;
    j_min = 1;
    i_max = i;
    j_max = j;

    //add 1 everywhere to enable indexing from 1 to i_max, j_max
    x.resize(i_max+1,std::vector<long double>(j_max+1,0));
    y.resize(i_max+1,std::vector<long double>(j_max+1,0));
    x_xi.resize(i_max+1,std::vector<long double>(j_max+1,0));
    y_xi.resize(i_max+1,std::vector<long double>(j_max+1,0));
    x_eta.resize(i_max+1,std::vector<long double>(j_max+1,0));
    y_eta.resize(i_max+1,std::vector<long double>(j_max+1,0));
    xi_x.resize(i_max+1,std::vector<long double>(j_max+1,0));
    xi_y.resize(i_max+1,std::vector<long double>(j_max+1,0));
    eta_x.resize(i_max+1,std::vector<long double>(j_max+1,0));
    eta_y.resize(i_max+1,std::vector<long double>(j_max+1,0));
    J.resize(i_max+1,std::vector<long double>(j_max+1,0));

    initGrid(fname);
};