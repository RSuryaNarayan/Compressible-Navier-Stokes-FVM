// get left and right states
Eigen::Matrix<long double, 4, 1> U_L  = cells[face.i_L].U;
Eigen::Matrix<long double, 4, 1> U_R  = cells[face.i_R].U;

Eigen::Matrix<long double, 5, 1> V_L;
Eigen::Matrix<long double, 5, 1> V_R;

//convert to primitive
U2V(U_L, V_L);
U2V(U_R, V_R);

//average onto the face for modified steger warming
// Eigen::Matrix<long double, 4, 1> U_f  = (U_L + U_R)/2.0L;
// Eigen::Matrix<long double, 5, 1> V_f = (V_L + V_R)/2.0L;
// U2V(U_f, V_f);

//do steger warming instead
long double rho_f = V_L(QRHO);
long double u_f = V_L(QU);
long double v_f = V_L(QV);
long double p_f = V_L(QP);
long double a_f = sqrt(GAMMA*p_f/rho_f);

//get unit normals of the face
long double n_x = face.n_x;
long double n_y = face.n_y;

//calculate u_dash and eigen values
long double ud = u_f * n_x + v_f * n_y;
long double lam = ud;
long double lam_p= ud + a_f;
long double lam_m= ud - a_f;

// eigen value matrix
Eigen::Matrix<long double, 4, 4> Lambda = Eigen::Matrix<long double, 4, 4>::Zero();
Lambda(0, 0) = lam;
Lambda(1, 1) = lam_p;
Lambda(2, 2) = lam;
Lambda(3, 3) = lam_m;

//split the eigen value matrix into plus and minus components
Eigen::Matrix<long double, 4, 4> Lambda_plus  = (Lambda + Lambda.cwiseAbs()) / 2.0L;
Eigen::Matrix<long double, 4, 4> Lambda_minus = (Lambda - Lambda.cwiseAbs()) / 2.0L;

// compute Aplus
long double lambda = Lambda_plus(0,0);
long double lambda_plus = Lambda_plus(1,1);
long double lambda_minus = Lambda_plus(3,3);
long double lambda_tilda = 0.5L*(lambda_plus-lambda_minus);
long double lambda_hat = 0.5L*(lambda_plus+lambda_minus)-lambda;

long double Prho = 0.5L*(GAMMA-1)*(u_f*u_f+v_f*v_f);
long double PE = (GAMMA-1);
long double h0 = (U_L(UEINT)+p_f)/rho_f;

Eigen::Matrix<long double, 4, 1> colvec1 = Eigen::Matrix<long double, 4, 1>::Zero();
colvec1 <<
lambda_hat/a_f/a_f,
(u_f*lambda_hat+a_f*n_x*lambda_tilda)/a_f/a_f,
(v_f*lambda_hat+a_f*n_y*lambda_tilda)/a_f/a_f,
(h0*lambda_hat +a_f*ud*lambda_tilda)/a_f/a_f;

Eigen::Matrix<long double, 1, 4> rowvec1 = Eigen::Matrix<long double, 1, 4>::Zero();
rowvec1 << Prho, -u_f*PE, -v_f*PE, PE;

Eigen::Matrix<long double, 4, 1> colvec2 = Eigen::Matrix<long double, 4, 1>::Zero();
colvec2 <<
lambda_tilda/a_f,
u_f*lambda_tilda/a_f+n_x*lambda_hat,
v_f*lambda_tilda/a_f+n_y*lambda_hat,
h0 *lambda_tilda/a_f+ud*lambda_hat;

Eigen::Matrix<long double, 1, 4> rowvec2 = Eigen::Matrix<long double, 1, 4>::Zero();
rowvec2 << -ud, n_x, n_y, 0;

Eigen::Matrix<long double, 4, 4> lam_mat =  Eigen::Matrix<long double, 4, 4>::Zero();
lam_mat(0,0) = lambda; lam_mat(1,1) = lambda; lam_mat(2,2) = lambda; lam_mat(3,3) = lambda;

//fill Jacobians in the face
face.Ap = colvec1*rowvec1 + colvec2*rowvec2 + lam_mat;

// compute Aminus
//do steger warming instead
rho_f = V_R(QRHO);
u_f = V_R(QU);
v_f = V_R(QV);
p_f = V_R(QP);
a_f = sqrt(GAMMA*p_f/rho_f);

//calculate u_dash and eigen values
ud = u_f * n_x + v_f * n_y;
lam = ud;
lam_p= ud + a_f;
lam_m= ud - a_f;

// eigen value matrix
Lambda = Eigen::Matrix<long double, 4, 4>::Zero();
Lambda(0, 0) = lam;
Lambda(1, 1) = lam_p;
Lambda(2, 2) = lam;
Lambda(3, 3) = lam_m;

//split the eigen value matrix into plus and minus components
Lambda_minus  = (Lambda - Lambda.cwiseAbs()) / 2.0L;

lambda = Lambda_minus(0,0);
lambda_plus = Lambda_minus(1,1);
lambda_minus = Lambda_minus(3,3);
lambda_tilda = 0.5L*(lambda_plus-lambda_minus);
lambda_hat = 0.5L*(lambda_plus+lambda_minus)-lambda;

Prho = 0.5L*(GAMMA-1)*(u_f*u_f+v_f*v_f);
PE = (GAMMA-1);
h0 = (U_R(UEINT)+p_f)/rho_f;

colvec1 = Eigen::Matrix<long double, 4, 1>::Zero();
colvec1 <<
lambda_hat/a_f/a_f,
(u_f*lambda_hat+a_f*n_x*lambda_tilda)/a_f/a_f,
(v_f*lambda_hat+a_f*n_y*lambda_tilda)/a_f/a_f,
(h0*lambda_hat +a_f*ud*lambda_tilda)/a_f/a_f;

rowvec1 = Eigen::Matrix<long double, 1, 4>::Zero();
rowvec1 << Prho, -u_f*PE, -v_f*PE, PE;

colvec2 = Eigen::Matrix<long double, 4, 1>::Zero();
colvec2 <<
lambda_tilda/a_f,
u_f*lambda_tilda/a_f+n_x*lambda_hat,
v_f*lambda_tilda/a_f+n_y*lambda_hat,
h0 *lambda_tilda/a_f+ud*lambda_hat;

rowvec2 = Eigen::Matrix<long double, 1, 4>::Zero();
rowvec2 << -ud, n_x, n_y, 0;

lam_mat =  Eigen::Matrix<long double, 4, 4>::Zero();
lam_mat(0,0) = lambda; lam_mat(1,1) = lambda; lam_mat(2,2) = lambda; lam_mat(3,3) = lambda;

//fill Jacobians in the face
face.Am = colvec1*rowvec1 + colvec2*rowvec2 + lam_mat;

//fill flux
face.F  = face.Ap * U_L + face.Am * U_R;

// we can now assign the fluxes and jacobian terms to the cells. 
// Note that by default our fluxes are in the outward normal direction for the L cell and 
// in the inward normal direction for the R cells. So flip the flux sign for the R cell

// fill the diagonal matrix Di = I*V_i/delta t + sum(A_plus*S_f)
cells[face.i_L].Di  = cells[face.i_L].Di + face.Ap * face.S_f;
cells[face.i_R].Di  = cells[face.i_R].Di - face.Am * face.S_f;

cells[face.i_L].Di_guess  = cells[face.i_L].Di_guess + face.Ap * face.S_f;
cells[face.i_R].Di_guess  = cells[face.i_R].Di_guess - face.Am * face.S_f;

//fill the off-diagonal jacobians and their locations.
cells[face.i_L].neighbors.push_back(face.i_R);
cells[face.i_L].N.push_back(face.Am*face.S_f);
cells[face.i_R].neighbors.push_back(face.i_L);
cells[face.i_R].N.push_back(-1*face.Ap*face.S_f);

//fill the residual Ri = -sum_faces(F*S_f).
cells[face.i_R].Ri  = cells[face.i_R].Ri + face.F * face.S_f;
cells[face.i_L].Ri  = cells[face.i_L].Ri - face.F * face.S_f;