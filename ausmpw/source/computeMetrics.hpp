void grid::computeMetrics()
{
    /*---------------------------------------------*\
                    xi derivatives
    \*---------------------------------------------*/
    for (int i=i_min; i<=i_max; i++)
    {
        for (int j=j_min; j<=j_max; j++)
        {
            if(i==i_min)
            {
                x_xi[i][j] = (4*x[i+1][j] - 3*x[i][j] - x[i+2][j])/2.0;
                y_xi[i][j] = (4*y[i+1][j] - 3*y[i][j] - y[i+2][j])/2.0;
            }
            else if(i==i_max)
            {
                x_xi[i][j] = (3*x[i][j] - 4*x[i-1][j] + x[i-2][j])/2.0;
                y_xi[i][j] = (3*y[i][j] - 4*y[i-1][j] + y[i-2][j])/2.0;
            }
            else
            {
                x_xi[i][j] = (x[i+1][j] - x[i-1][j])/2.0;
                y_xi[i][j] = (y[i+1][j] - y[i-1][j])/2.0;
            }
        }
    }
    
    /*---------------------------------------------*\
                    eta derivatives
    \*---------------------------------------------*/
    for (int i=i_min; i<=i_max; i++)
    {
        for (int j=j_min; j<=j_max; j++)
        {
            if(j==j_min)
            {
                x_eta[i][j] = (4*x[i][j+1] - 3*x[i][j] - x[i][j+2])/2.0;
                y_eta[i][j] = (4*y[i][j+1] - 3*y[i][j] - y[i][j+2])/2.0;
            }
            else if(j==j_max)
            {
                x_eta[i][j] = (3*x[i][j] - 4*x[i][j-1] + x[i][j-2])/2.0;
                y_eta[i][j] = (3*y[i][j] - 4*y[i][j-1] + y[i][j-2])/2.0;
            }
            else
            {
                x_eta[i][j] = (x[i][j+1] - x[i][j-1])/2.0;
                y_eta[i][j] = (y[i][j+1] - y[i][j-1])/2.0;
            }
        }
    }

    //inverse grid metrics and the jacobian 

    for (int i=i_min; i<=i_max; i++)
    {
        for (int j=j_min; j<=j_max; j++)
        {
            J[i][j]     =  1/(x_xi[i][j]*y_eta[i][j] - x_eta[i][j]*y_xi[i][j]);
            xi_x[i][j]  =  y_eta[i][j]*J[i][j]; 
            xi_y[i][j]  = -x_eta[i][j]*J[i][j];
            eta_x[i][j] = - y_xi[i][j]*J[i][j];
            eta_y[i][j] =   x_xi[i][j]*J[i][j];
        }
    }
    writeFields();
}