// get left and right states
Eigen::Matrix<long double, 4, 1> U_L; 
Eigen::Matrix<long double, 5, 1> V_L; 
Eigen::Matrix<long double, 4, 1> U_R  = cells[face.i_R].U;
Eigen::Matrix<long double, 5, 1> V_R  = cells[face.i_R].V;


//derive ghost cell velocities
long double u_i = cells[face.i_R].U(UMX)/cells[face.i_R].U(URHO);
long double v_i = cells[face.i_R].U(UMY)/cells[face.i_R].U(URHO);
long double p_i = cells[face.i_R].V(QP);
long double ud = u_i * face.n_x + v_i * face.n_y;
long double u_ghost = -u_i;//u_i - 2*ud*face.n_x;
long double v_ghost = -v_i;//v_i - 2*ud*face.n_y;
long double p_ghost = cells[face.i_R].U(URHO) * R * (2*T_wall - V_R(QT));
long double E_ghost = p_ghost/(GAMMA-1) + 0.5 * cells[face.i_R].U(URHO) * (u_ghost *u_ghost + v_ghost * v_ghost);

//set left state
U_L(URHO) = U_R(URHO);
U_L(UMX) = U_R(URHO)*u_ghost;
U_L(UMY) = U_R(URHO)*v_ghost;
U_L(UEINT) = E_ghost; 

// V_L(QRHO) = V_R(QRHO);
// V_L(QU) = 0.0;
// V_L(QV) = 0.0;
// V_L(QP) = V_R(QRHO) * R * (2*T_wall -V_R(QT));
// V_L(QT) = 2*T_wall - V_R(QT);
// V2U(V_L, U_L);

//average onto the face for modified steger warming
Eigen::Matrix<long double, 4, 1> U_f = (U_L+U_R)/2.0;
Eigen::Matrix<long double, 5, 1> V_f;
U2V(U_f, V_f);

long double rho_f = V_f(QRHO);
long double u_f = V_f(QU);
long double v_f = V_f(QV);
long double T_f = V_f(QT);
long double p_f = rho_f * R * T_f;
long double a_f = sqrt(GAMMA*p_f/rho_f);
U_f(3) = p_f/(GAMMA-1) + 0.5 * rho_f * (u_f * u_f + v_f * v_f);

//get unit normals of the face
long double n_x = face.n_x;
long double n_y = face.n_y;

//calculate u_dash and eigen values
ud = u_f * n_x + v_f * n_y;
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
long double h0 = (U_f(UEINT)+p_f)/rho_f;

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
lambda = Lambda_minus(0,0);
lambda_plus = Lambda_minus(1,1);
lambda_minus = Lambda_minus(3,3);
lambda_tilda = 0.5L*(lambda_plus-lambda_minus);
lambda_hat = 0.5L*(lambda_plus+lambda_minus)-lambda;

Prho = 0.5L*(GAMMA-1)*(u_f*u_f+v_f*v_f);
PE = (GAMMA-1);
h0 = (U_f(UEINT)+p_f)/rho_f;

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

//fill the diagonal matrix Di = I*V_i/delta t + sum(A_plus*S_f)
// Eigen::Matrix<long double, 4, 4> Ew = Eigen::Matrix<long double, 4, 4>::Zero();

//for inviscid
// Ew(0,0)=1.0L; Ew(1,1)=1.0L-2.0L*face.n_x*face.n_x; Ew(1,2) = -2.0L*face.n_x*face.n_y; 
// Ew(2,1) = -2.0L*face.n_x*face.n_y; Ew(2,2)=1.0L-2*face.n_y*face.n_y; Ew(3,3)=1.0L;

//for viscous
// Ew(0,0) =1.0L; Ew(1,1)=-1.0L; Ew(2,2)=-1.0L; Ew(3,3)=-1.0;

//compute ghost cell transformation Jacobian N_L
//get Jacobian N_L
// long double rho_L = V_L(QRHO);
// long double u_L = V_L(QU);
// long double v_L = V_L(QV);
// long double T_L = V_L(QT);

// Eigen::Matrix<long double, 4, 4> N_L = Eigen::Matrix<long double, 4, 4>::Zero();
// N_L << 
// 1.0L, 0.0L, 0.0L, 0.0L,
// -u_L/rho_L, 1/rho_L, 0.0L, 0.0L,
// -v_L/rho_L, 0.0L, 1/rho_L, 0.0L,
// (u_L*u_L+v_L*v_L)/2.0L/CV/rho_L - T_L/rho_L, -u_L/CV/rho_L, -v_L/CV/rho_L, 1/CV/rho_L;

// //compute interior cell transformation Jacobian N_R
// long double rho_R = V_R(QRHO);
// long double u_R = V_R(QU);
// long double v_R = V_R(QV);
// long double T_R = V_R(QT);

// Eigen::Matrix<long double, 4, 4> N_R = Eigen::Matrix<long double, 4, 4>::Zero();
// N_R << 
// 1.0L, 0.0L, 0.0L, 0.0L,
// -u_R/rho_R, 1/rho_R, 0.0L, 0.0L,
// -v_R/rho_R, 0.0L, 1/rho_R, 0.0L,
// (u_R*u_R+v_R*v_R)/2.0L/CV/rho_R - T_R/rho_R, -u_R/CV/rho_R, -v_R/CV/rho_R, 1/CV/rho_R;

// cells[face.i_R].Di  = cells[face.i_R].Di + face.Ap * face.S_f + face.Am * N_L.inverse() * Ew * N_R * face.S_f;
// cells[face.i_R].Di_guess  = cells[face.i_R].Di_guess + face.Ap * face.S_f;

//fill the residual Ri = -sum_faces(F*S_f).
cells[face.i_R].Ri  = cells[face.i_R].Ri + face.F * face.S_f;

// Viscous flux computation

//calculate lambda, mu and k 
long double mu = 1.458 * 1e-6 * std::pow(T_f, 1.5) / (T_f + 110.3);
lambda = -(2.0L/3.0L) * mu;
long double kappa = (CP)*mu/0.7;

// // need M_xx, M_xy, M_yx, M_yy and gradient
// Eigen::Matrix<long double, 4, 4> M_xx = Eigen::Matrix<long double, 4, 4>::Zero();
// M_xx << 
// 0.0L, 0.0L, 0.0L, 0.0L, 
// 0.0L, lambda+2*mu, 0.0L, 0.0L,
// 0.0L, 0.0L, mu, 0.0L,
// 0.0, u_f*(lambda+2*mu),  v_f*mu, kappa;

// Eigen::Matrix<long double, 4, 4> M_xy = Eigen::Matrix<long double, 4, 4>::Zero();
// M_xy << 
// 0.0L, 0.0L, 0.0L, 0.0L, 
// 0.0L, 0.0L, lambda, 0.0L,
// 0.0L, mu, 0.0L, 0.0L,
// 0.0, v_f*mu,  u_f*lambda, 0.0L;

// Eigen::Matrix<long double, 4, 4> M_yx = Eigen::Matrix<long double, 4, 4>::Zero();
// M_yx << 
// 0.0L, 0.0L, 0.0L, 0.0L, 
// 0.0L, 0.0L, mu, 0.0L,
// 0.0L, lambda, 0.0L, 0.0L,
// 0.0, v_f*lambda,  u_f*mu, 0.0L;

// Eigen::Matrix<long double, 4, 4> M_yy = Eigen::Matrix<long double, 4, 4>::Zero();
// M_yy << 
// 0.0L, 0.0L, 0.0L, 0.0L, 
// 0.0L, mu, 0.0, 0.0L,
// 0.0L, 0.0, lambda+2*mu, 0.0L,
// 0.0, u_f*mu,  v_f*(lambda+2*mu), kappa;

// //get M_dash
// Eigen::Matrix<long double, 4, 4> M_dash = Eigen::Matrix<long double, 4, 4>::Zero();
// M_dash = M_xx * face.n_x * face.n_x + (M_xy + M_yx)*face.n_x * face.n_y + M_yy * face.n_y * face.n_y;

// //get Jacobian N 
// Eigen::Matrix<long double, 4, 4> N = Eigen::Matrix<long double, 4, 4>::Zero();
// N << 
// 1.0L, 0.0L, 0.0L, 0.0L,
// -u_f/rho_f, 1/rho_f, 0.0L, 0.0L,
// -v_f/rho_f, 0.0L, 1/rho_f, 0.0L,
// (u_f*u_f+v_f*v_f)/2.0L/CV/rho_f - T_f/rho_f, -u_f/CV/rho_f, -v_f/CV/rho_f, 1/CV/rho_f;

// //get gradient
// Eigen::Matrix<long double, 4, 1> grad_U = Eigen::Matrix<long double, 4, 1>::Zero();
// grad_U  = (U_R - U_L)/sqrt( (cells[face.i_R+1].x_c-cells[face.i_R].x_c)*(cells[face.i_R+1].x_c-cells[face.i_R].x_c) 
//                           + (cells[face.i_R+1].y_c-cells[face.i_R].y_c)*(cells[face.i_R+1].y_c-cells[face.i_R].y_c) );

// Eigen::Matrix<long double, 4, 1> F_v= Eigen::Matrix<long double, 4, 1>::Zero();
// F_v = -1.0L*M_dash * N * grad_U;

//viscous flux using green gauss gradients

Eigen::Matrix<long double, 2, 1> gradU = (cells[face.i_R].gradU);
Eigen::Matrix<long double, 2, 1> gradV = (cells[face.i_R].gradV);
Eigen::Matrix<long double, 2, 1> gradT = (cells[face.i_R].gradT);

long double sigma11 = 2*mu*(gradU(0)) - (2/3)*mu*(gradU(0)+gradV(1));
long double sigma12 = mu*(gradU(1)+gradV(0));
long double sigma21 = mu*(gradU(1)+gradV(0));
long double sigma22 = 2*mu*(gradV(0)) - (2/3)*mu*(gradU(0)+gradV(1));
long double q1 = kappa*gradT(0);
long double q2 = kappa*gradT(1);

Eigen::Matrix<long double, 4, 1> vec1 = Eigen::Matrix<long double, 4, 1>::Zero();
vec1<<0,-sigma11,-sigma12,-(sigma11*u_f+sigma21*v_f)-q1;

Eigen::Matrix<long double, 4, 1> vec2 = Eigen::Matrix<long double, 4, 1>::Zero();
vec2<<0,-sigma12,-sigma22,-(sigma12*u_f+sigma22*v_f)-q2;

Eigen::Matrix<long double, 4, 1> F_v= Eigen::Matrix<long double, 4, 1>::Zero();
F_v = (vec1 * face.n_x + vec2 * face.n_y) * face.S_f;

cells[face.i_R].Ri  = cells[face.i_R].Ri + F_v * face.S_f;