// #include <Eigen/Sparse>
// #include <iostream>

// int main() {
//     // Define a sparse matrix A
//     Eigen::SparseMatrix<double> A(3, 3);
//     A.insert(0, 0) = 4;
//     A.insert(0, 1) = -1;
//     A.insert(1, 0) = -1;
//     A.insert(1, 1) = 6;
//     A.insert(1, 2) = 1;
//     A.insert(2, 1) = 1;
//     A.insert(2, 2) = 5;

//     // Define the right-hand side vector B
//     Eigen::VectorXd B(3);
//     B << 7, 4, 3;

//     // Use SparseLU to solve the system
//     Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
//     solver.compute(A);

//     if (solver.info() != Eigen::Success) {
//         std::cerr << "Decomposition failed!" << std::endl;
//         return -1;
//     }

//     Eigen::VectorXd x = solver.solve(B);

//     if (solver.info() != Eigen::Success) {
//         std::cerr << "Solving failed!" << std::endl;
//         return -1;
//     }

//     std::cout << "Solution x:\n" << x << std::endl;

//     return 0;
// }

#include <iostream>
#include <Eigen/Dense>

int main() {
    // Define matrices with long double type and dynamic sizes
    Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> mat1(2, 3); // 2x3 matrix
    Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> mat2(3, 2); // 3x2 matrix

    // Initialize matrices
    mat1 << 1.1L, 2.2L, 3.3L,
            4.4L, 5.5L, 6.6L;

    mat2 << 7.7L, 8.8L,
            9.9L, 10.10L,
            11.11L, 12.12L;

    // Multiply the matrices
    Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> result = mat1 * mat2;

    // Print the result
    std::cout << "Result of multiplication:\n" << result << std::endl;

    return 0;
}
