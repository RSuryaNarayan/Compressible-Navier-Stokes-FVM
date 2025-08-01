/* Contains the main solver class definition */
#ifndef GRID_HPP
#define GRID_HPP

#include<iostream>
#include<string>
#include<vector>

class grid{
    protected:
        int i_min, i_max, j_min, j_max;
        std::vector<std::vector<long double>> x, y, x_xi, y_xi, x_eta, y_eta, J, xi_x, xi_y, eta_x, eta_y;
        std::vector<long double> x1D, y1D;
    public:
        grid(int,int, std::string);
        void initGrid(std::string);
        void writeFields();
        void computeMetrics();
};

#endif