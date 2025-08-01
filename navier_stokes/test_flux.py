import numpy as np

def InteriorFaceFlux():

    gamma = 1.4
    #left cell state
    rho_L = 0.266
    u_L = 0.927
    v_L = 0.0
    P_L = 0.303
    a_L = np.sqrt(gamma*P_L/rho_L)
    E_L = P_L/(gamma-1)+0.5*rho_L*(u_L*u_L)

    #right cell state
    rho_R = 0.125
    u_R = 0.0
    v_R = 0.0
    P_R = 0.1
    a_R = np.sqrt(gamma*P_R/rho_R)
    E_R = P_R/(gamma-1)+0.5*rho_R*(u_R*u_R)

    #modified steger warming
    #average on conserved state
    rho_f = 0.5*(rho_L+rho_R)
    rhou_f = 0.5*(rho_L*u_L+rho_R*u_R)
    rhov_f = 0.5*(rho_L*v_L+rho_R*v_R)
    E_f = 0.5*(E_L+E_R)

    #convert to prim
    u_f = rhou_f/rho_f
    v_f = rhov_f/rho_f
    p_f = (E_f - 0.5 * (rhou_f**2 + rhov_f**2)/rho_f) * (gamma - 1)
    a_f = np.sqrt(gamma*p_f/rho_f)

    #both left and right will be evaluated at this face value now
    rho_L = rho_f
    u_L = u_f
    v_L = v_f
    a_L = a_f
    P_L = p_f
    E_L = E_f

    rho_R = rho_f
    u_R = u_f
    v_R = v_f
    a_R = a_f
    P_R = p_f
    E_R = E_f

    nx = 1.0
    ny = 0.0

    u_dL = u_L*nx + v_L*ny
    u_dR = u_R*nx + v_R*ny

    #left eigen values 
    if (u_dL > 0):
        lambda_L = u_dL
    else:
        lambda_L = 0
    
    if ( (u_dL + a_L) > 0):
        lambdap_L = u_dL + a_L
    else:
        lambdap_L = 0

    if ( (u_dL - a_L) > 0):
        lambdam_L = u_dL - a_L
    else:
        lambdam_L = 0
    
    #right eigen values
    if (u_dR < 0):
        lambda_R = u_dR
    else:
        lambda_R = 0
    
    if ( (u_dR + a_R) < 0):
        lambdap_R = u_dR + a_R
    else:
        lambdap_R = 0

    if ( (u_dR - a_R) < 0):
        lambdam_R = u_dR - a_R
    else:
        lambdam_R = 0

    #set l_t, l_c
    #left cell
    l_tL = 0.5*(lambdap_L-lambdam_L)
    l_cL = 0.5*(lambdap_L+lambdam_L) - lambda_L

    #right cell
    l_tR = 0.5*(lambdap_R-lambdam_R)
    l_cR = 0.5*(lambdap_R+lambdam_R) - lambda_R

    # print("int flux")
    # print(l_tL)
    # print(l_cL)
    # print(l_tR)
    # print(l_cR)

    #Prho, PE, h0
    Prho_L = 0.5*(gamma - 1)*(u_L**2 + v_L**2)
    PE_L = (gamma - 1)
    h0_L = (E_L + P_L)/rho_L

    Prho_R = 0.5*(gamma - 1)*(u_R**2 + v_R**2)
    PE_R = (gamma - 1)
    h0_R = (E_R + P_R)/rho_R

    # print("int flux")
    # print(Prho_L)
    # print(PE_L)
    # print(h0_L)
    # print(Prho_R)
    # print(PE_R)
    # print(h0_R)
    
    #forming matrices for jacobian estimation
    #A+
    column1 = np.array( [ [l_cL/(a_L**2) , ( u_L*l_cL + a_L*nx*l_tL )/(a_L**2), ( v_L*l_cL + a_L*ny*l_tL )/(a_L**2), ( h0_L*l_cL + a_L*u_dL*l_tL )/(a_L**2)]]).T
    
    row1 = np.array( [ [  Prho_L  ,   -u_L*PE_L ,  -v_L*PE_L ,  PE_L ] ])

    column2 = np.array( [ [l_tL/a_L],
                          [u_L*l_tL/a_L + nx*l_cL], 
                          [v_L*l_tL/a_L + ny*l_cL],  
                          [h0_L*l_tL/a_L + u_dL*l_cL] 
                        ])
    
    row2 = np.array(   [ [   -u_dL,  nx, ny, 0   ] ]  )

    mat = lambda_L*np.identity(4)

    # print(column1.shape)
    # print("\n")
    # print(row1.shape)
    # print("\n")
    # print(column2.shape)
    # print("\n")
    # print(row2.shape)
    # print("\n")

    Ap = np.matmul(column1,row1) + np.matmul(column2,row2) + mat

    #A-
    column11 = np.array( [ [l_cR/(a_R**2)],
                            [ ( u_R*l_cR + a_R*nx*l_tR )/(a_R**2)],
                            [ ( v_R*l_cR + a_R*ny*l_tR )/(a_R**2)],
                            [ ( h0_R*l_cR + a_R*u_dR*l_tR )/(a_R**2)]
                            ])
    
    row11 = np.array( [ [  Prho_R  ,   -u_R*PE_R ,  -v_R*PE_R ,  PE_R ] ])

    column22 = np.array( [ [l_tR/a_R],
                          [u_R*l_tR/a_R + nx*l_cR], 
                          [v_R*l_tR/a_R + ny*l_cR],  
                          [h0_R*l_tR/a_R + u_dR*l_cR] 
                        ])
    
    row22 = np.array(  [  [   -u_dR,  nx, ny, 0   ] ] )

    mat1 = lambda_R*np.identity(4)

    Am = np.matmul(column11,row11) + np.matmul(column22,row22) + mat1

    print(Am)
    #left cell state
    rho_L = 0.266
    u_L = 0.927
    v_L = 0.0
    P_L = 0.303
    a_L = np.sqrt(gamma*P_L/rho_L)
    E_L = P_L/(gamma-1)+0.5*rho_L*(u_L*u_L)

    #right cell state
    rho_R = 0.125
    u_R = 0.0
    v_R = 0.0
    P_R = 0.1
    a_R = np.sqrt(gamma*P_R/rho_R)
    E_R = P_R/(gamma-1)+0.5*rho_R*(u_R*u_R)

    #estimate flux
    U_L = np.array([    rho_L,
                        rho_L*u_L,
                        rho_L*v_L,
                        P_L/(gamma-1)+0.5*rho_L*(u_L*u_L)
                    ])
    
    #estimate flux
    U_R = np.array([    rho_R,
                        rho_R*u_R,
                        rho_R*v_R,
                        P_R/(gamma-1)+0.5*rho_R*(u_R*u_R)
                    ])
    
    print("U_L=",U_L)
    print("U_R=",U_R)
    # print(Ap.shape)
    # print("\n")
    # print(U_L.shape)
    # print("\n")
    # print(Am.shape)
    # print("\n")
    # print(U_R.shape)
    # print("\n")
    
    flux = ( np.matmul(Ap,U_L) + np.matmul(Am,U_R) ) * 1.0
    print(flux)


InteriorFaceFlux()