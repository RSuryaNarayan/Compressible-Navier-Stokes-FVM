void solver::compute_flux()
{
    //calculate F -- the x-direction flux
    for (int i=i_min; i<=i_max-1; i++)
    {
        for (int j=j_min+1; j<=j_max-1; j++)
        {
            // get the extrapolated Left and Right states at a face
            for (int k=0; k<NUVAR; k++)
            {
                if (i==i_min || i==i_max-1)
                {
                    //At boundaries, we will drop to an upwind formulation
                    UL[i][j][k] = Un[i][j][k];
                    UR[i][j][k] = Un[i+1][j][k];
                }
                else
                {
                    
                    if (use_flux_lim==1)
                    {
                        //left state
                        double d1 = minmod( (Un[i][j][k] - Un[i-1][j][k]), beta * (Un[i+1][j][k] - Un[i][j][k]) );
                        double d2 = minmod( (Un[i+1][j][k] - Un[i][j][k]), beta * (Un[i][j][k] - Un[i-1][j][k]) );

                        UL[i][j][k] = Un[i][j][k]   + (1/4)*( (1-k)*d1 + (1+k)*d2 );
                                                    

                        //right state
                        d1 = minmod( (Un[i+1][j][k] - Un[i][j][k]), beta * (Un[i+2][j][k] - Un[i+1][j][k]) );
                        d2 = minmod( (Un[i+2][j][k] - Un[i+1][j][k]), beta * (Un[i+1][j][k] - Un[i][j][k]) );

                        UR[i][j][k] = Un[i+1][j][k] - (1/4)*( (1+k) * d1 + (1-k) * d2 );
                    }
                    else 
                    {
                        UL[i][j][k] = Un[i][j][k]   + 0.5 * eps * ( Un[i][j][k]   - Un[i-1][j][k] );
                        UR[i][j][k] = Un[i+1][j][k] + 0.5 * eps * ( Un[i+1][j][k] - Un[i+2][j][k] );
                    }

                }
            }

            //convert to corresponding primitive forms
            cons_to_prim(UL[i][j], QL[i][j]); 
            cons_to_prim(UR[i][j], QR[i][j]);

            //h0 at left and right
            long double h0L = h0(QL[i][j]);
            long double h0R = h0(QR[i][j]);

            //calculate metrics at the faces
            long double xi_xavg  = 0.5 * (xi_x[i][j]  + xi_x[i+1][j]  );
            long double xi_yavg  = 0.5 * (xi_y[i][j]  + xi_y[i+1][j]  );
            long double eta_xavg = 0.5 * (eta_x[i][j] + eta_y[i+1][j] );
            long double eta_yavg = 0.5 * (eta_y[i][j] + eta_y[i+1][j] );
            long double Javg     = 0.5 * (    J[i][j] +     J[i+1][j] );

            //normalization factors
            long double A1 = std::sqrt( xi_xavg  * xi_xavg   + xi_yavg  * xi_yavg  );
            long double A2 = std::sqrt( eta_xavg * eta_xavg  + eta_yavg * eta_yavg );
            
            //normalized contravarient velocities
            long double utL = (xi_xavg  * QL[i][j][QU]  + xi_yavg  * QL[i][j][QV]) /A1;
            long double vtL = (eta_xavg * QL[i][j][QU]  + eta_yavg * QL[i][j][QV]) /A2;
            long double utR = (xi_xavg  * QR[i][j][QU]  + xi_yavg  * QR[i][j][QV]) /A1;
            long double vtR = (eta_xavg * QR[i][j][QU]  + eta_yavg * QR[i][j][QV]) /A2;

            //calculate h0norm and cs
            long double h0norm = 0.5 * (h0L - 0.5*vtL*vtL + h0R - 0.5*vtR*vtR);
            long double cs     = std::sqrt(2*h0norm*(GAMMA-1)/(GAMMA+1));
            long double cavg   = 0.5*(utL+utR)>=0? (cs*cs)/std::max(std::abs(utL),cs) : (cs*cs)/std::max(std::abs(utR),cs);

            //cell face mach numbers 
            long double MtL = utL/cavg;
            long double MtR = utR/cavg;

            //split Mach numbers
            long double MtLp = std::abs(MtL)<=1?    0.25*std::pow((MtL+1),2) : 0.5*(MtL+std::abs(MtL));
            long double MtRm = std::abs(MtR)<=1? -1*0.25*std::pow((MtR-1),2) : 0.5*(MtR-std::abs(MtR));

            //split pressures
            long double Pp = std::abs(MtL)<=1?  0.25*pow((MtL+1),2)*(2-MtL) : 0.5*(1 + sign(MtL));
            long double Pm = std::abs(MtR)<=1?  0.25*pow((MtR-1),2)*(2+MtR) : 0.5*(1 - sign(MtR));

            //calculate pressure weighting terms
            long double pL = QL[i][j][QPRES];
            long double pR = QR[i][j][QPRES];
            long double p_min = findMin(Qn[i][j-1][QPRES], Qn[i][j+1][QPRES], Qn[i+1][j-1][QPRES], Qn[i+1][j+1][QPRES]);
            long double ps = Pp * pL + Pm * pR;
            long double min1 = pL<=pR? pL : pR;
            long double min2 = 1 <= p_min/min1 ? 1 : p_min/min1;
            long double fL = ps<=0? 0 : (pL/ps-1)*std::pow(min2,2);
            long double fR = ps<=0? 0 : (pR/ps-1)*std::pow(min2,2);
            long double min3 = pL/pR <= pR/pL ?  pL/pR : pR/pL;
            long double w = 1 - std::pow(min3, 3);

            //compute averged mach number
            long double Mtavg = MtLp + MtRm;

            //compute weighting mach numbers
            long double MtLpB = Mtavg>=0? MtLp + MtRm*( (1-w)*(1+fR)-fL) : MtLp*w*(1+fL);
            long double MtRmB = Mtavg>=0? MtRm*w*(1+fR) : MtRm + MtLp*( (1-w)*(1+fL)-fR);

            //compute fluxes
            F[i][j][URHO] = (MtLpB * cavg * A1/Javg) * UL[i][j][URHO]
                        +   (MtRmB * cavg * A1/Javg) * UR[i][j][URHO];
            
            F[i][j][UMX] = (MtLpB * cavg * A1/Javg) * UL[i][j][UMX]
                        +  (MtRmB * cavg * A1/Javg) * UR[i][j][UMX]
                        +  (Pp/Javg) * (xi_xavg * QL[i][j][QPRES])
                        +  (Pm/Javg) * (xi_xavg * QR[i][j][QPRES]);

            F[i][j][UMY] = (MtLpB * cavg * A1/Javg) * UL[i][j][UMY]
                        +  (MtRmB * cavg * A1/Javg) * UR[i][j][UMY]
                        +  (Pp/Javg) * (xi_yavg * QL[i][j][QPRES])
                        +  (Pm/Javg) * (xi_yavg * QR[i][j][QPRES]);

            F[i][j][UE] = (MtLpB * cavg * A1/Javg) * (UL[i][j][UE] + QL[i][j][QPRES])
                        + (MtRmB * cavg * A1/Javg) * (UR[i][j][UE] + QR[i][j][QPRES]);
        }
    }

    //calculate G -- the y-direction flux
    for (int i=i_min+1; i<=i_max-1; i++)
    {
        for (int j=j_min; j<=j_max-1; j++)
        {
            // get the extrapolated Bottom and Top states at a face
            for (int k=0; k<NUVAR; k++)
            {
                if (j==j_min || j==j_max-1)
                {
                    //At boundaries, we will drop to an upwind formulation
                    UB[i][j][k] = Un[i][j][k];
                    UT[i][j][k] = Un[i][j+1][k];
                }
                else
                {
                    if (use_flux_lim==1)
                    {
                        //left state
                        double d1 = minmod( (Un[i][j][k] - Un[i][j-1][k]), beta * (Un[i][j+1][k] - Un[i][j][k]) );
                        double d2 = minmod( (Un[i][j+1][k] - Un[i][j][k]), beta * (Un[i][j][k] - Un[i][j-1][k]) );

                        UB[i][j][k] = Un[i][j][k]   + (1/4)*( (1-k)*d1 + (1+k)*d2 );
                                                    

                        //right state
                        d1 = minmod( (Un[i][j+1][k] - Un[i][j][k]), beta * (Un[i][j+2][k] - Un[i][j+1][k]) );
                        d2 = minmod( (Un[i][j+2][k] - Un[i][j+1][k]), beta * (Un[i][j+1][k] - Un[i][j][k]) );

                        UT[i][j][k] = Un[i][j+1][k] - (1/4)*( (1+k) * d1 + (1-k) * d2 );
                    }
                    else 
                    {
                        UB[i][j][k] = Un[i][j][k]   + 0.5 * eps * ( Un[i][j][k]   - Un[i][j-1][k] );
                        UT[i][j][k] = Un[i][j+1][k] + 0.5 * eps * ( Un[i][j+1][k] - Un[i][j+2][k] );
                    }
                }
            }

            //convert to corresponding primitive forms
            cons_to_prim(UB[i][j], QB[i][j]); 
            cons_to_prim(UT[i][j], QT[i][j]);

            //h0 at left and right
            long double h0B = h0(QB[i][j]);
            long double h0T = h0(QT[i][j]);

            //calculate metrics at the faces
            long double xi_xavg  =  0.5 * ( xi_x [i][j]  + xi_x [i][j+1]  );
            long double xi_yavg  =  0.5 * ( xi_y [i][j]  + xi_y [i][j+1]  );
            long double eta_xavg =  0.5 * ( eta_x[i][j]  + eta_x[i][j+1]  );
            long double eta_yavg =  0.5 * ( eta_y[i][j]  + eta_y[i][j+1]  );
            long double Javg     =  0.5 * (     J[i][j]  +     J[i][j+1]  );

            //normalization factors
            long double A1 = std::sqrt( xi_xavg  * xi_xavg  + xi_yavg  * xi_yavg  );
            long double A2 = std::sqrt( eta_xavg * eta_xavg + eta_yavg * eta_yavg );
            
            //normalized contravarient velocities
            long double utB = (xi_xavg   * QB[i][j][QU] + xi_yavg  * QB[i][j][QV])/A1;
            long double vtB = (eta_xavg  * QB[i][j][QU] + eta_yavg * QB[i][j][QV])/A2;
            long double utT = (xi_xavg   * QT[i][j][QU] + xi_yavg  * QT[i][j][QV])/A1;
            long double vtT = (eta_xavg  * QT[i][j][QU] + eta_yavg * QT[i][j][QV])/A2;

            //calculate h0norm and cs
            long double h0norm = 0.5*(h0B - 0.5 * utB * utB + h0T - 0.5*utT*utT);
            long double cs     = std::sqrt(2 * h0norm * (GAMMA-1)/(GAMMA+1));
            long double cavg   = 0.5*(vtB+vtT)>=0? (cs*cs)/std::max(std::abs(vtB),cs) : (cs*cs)/std::max(std::abs(vtT),cs);

            //cell face mach numbers 
            long double MtB = vtB/cavg;
            long double MtT = vtT/cavg;

            //split Mach numbers
            long double MtBp = std::abs(MtB)<=1?    0.25*std::pow((MtB+1),2) : 0.5*(MtB+std::abs(MtB));
            long double MtTm = std::abs(MtT)<=1? -1*0.25*std::pow((MtT-1),2) : 0.5*(MtT-std::abs(MtT));

            //split pressures
            long double Pp = std::abs(MtB)<=1?  0.25*pow((MtB+1),2)*(2-MtB) : 0.5*(1 + sign(MtB));
            long double Pm = std::abs(MtT)<=1?  0.25*pow((MtT-1),2)*(2+MtT) : 0.5*(1 - sign(MtT));

            //calculate pressure weighting terms
            long double pB = QB[i][j][QPRES];
            long double pT = QT[i][j][QPRES];
            long double p_min = findMin(Qn[i-1][j][QPRES], Qn[i+1][j][QPRES], Qn[i-1][j+1][QPRES], Qn[i+1][j+1][QPRES]);
            long double ps = Pp * pB + Pm * pT;
            long double min1 = pB <= pT ? pB : pT;
            long double min2 = 1 <= p_min/min1 ? 1 : p_min/min1;
            long double fB = ps<=0? 0 : (pB/ps-1)*std::pow(min2,2);
            long double fT = ps<=0? 0 : (pT/ps-1)*std::pow(min2,2);
            long double min3 = pB/pT <= pT/pB ? pB/pT : pT/pB;
            long double w = 1 - std::pow(min3, 3);

            //compute averged mach number
            long double Mtavg = MtBp + MtTm;

            //compute weighting mach numbers
            long double MtBpB = Mtavg>=0? MtBp + MtTm * ( (1-w)*(1+fT)-fB ) : MtBp*w*(1+fB);
            long double MtTmB = Mtavg>=0? MtTm*w*(1+fT) : MtTm + MtBp * ( (1-w)*(1+fB)-fT);

            //compute fluxes
            G[i][j][URHO] = (MtBpB * cavg * A2/Javg) * UB[i][j][URHO]
                        +   (MtTmB * cavg * A2/Javg) * UT[i][j][URHO];
            
            G[i][j][UMX] = (MtBpB * cavg * A2/Javg) * UB[i][j][UMX]
                        +  (MtTmB * cavg * A2/Javg) * UT[i][j][UMX]
                        +  (Pp/Javg) * (eta_xavg * QB[i][j][QPRES])
                        +  (Pm/Javg) * (eta_xavg * QT[i][j][QPRES]);

            G[i][j][UMY] = (MtBpB * cavg * A2/Javg) * UB[i][j][UMY]
                        +  (MtTmB * cavg * A2/Javg) * UT[i][j][UMY]
                        +  (Pp/Javg) * (eta_yavg * QB[i][j][QPRES])
                        +  (Pm/Javg) * (eta_yavg * QT[i][j][QPRES]);

            G[i][j][UE] = (MtBpB * cavg * A2/Javg) * (UB[i][j][UE] + QB[i][j][QPRES])
                        + (MtTmB * cavg * A2/Javg) * (UT[i][j][UE] + QT[i][j][QPRES]);
        }
    }
    
    //optional code block to write out the fluxes to debug
    /*std::ofstream fout("F.txt");
    std::ofstream gout("G.txt");
    for (int i=i_min; i<=i_max; i++)
    {
        for (int j=j_min; j<=j_max; j++)
        {
            fout<<x[i][j]<<"\t"<<y[i][j]<<"\t"<<F[i][j][URHO]<<"\t"<<F[i][j][UMX]<<"\t"<<F[i][j][UMY]<<"\t"<<F[i][j][UE]<<"\n";
            gout<<x[i][j]<<"\t"<<y[i][j]<<"\t"<<G[i][j][URHO]<<"\t"<<G[i][j][UMX]<<"\t"<<G[i][j][UMY]<<"\t"<<G[i][j][UE]<<"\n";
        }
    }*/
}