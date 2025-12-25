#include <iostream>
#include "matrix.hpp"

int main(int argc, char** argv)
{
    if (argc < 4) {
        std::cout << "Usage:\n";
        std::cout << "  " << argv[0] << " <double|float> <method> <filename>\n";
        std::cout << "Methods:\n";
        std::cout << "  print    - вывести матрицу\n";
        std::cout << "  gauss    - метод Гаусса\n";             // работоспособен (2)
        std::cout << "  jordan   - метод Жордана\n";            // работоспособен (2)
        std::cout << "  lu       - LU-разложение\n";            // работоспособен (3)
        std::cout << "  tridiag  - метод прогонки (трёхдиагональные СЛАУ)\n";   // работоспособен (4)
        std::cout << "  sqrt     - метод квадратных корней\n";  // работоспособен
        std::cout << "  iter     - метод простых итераций\n";   // работоспособен
        std::cout << "  seidel   - метод Зейделя\n";            // работоспособен
        std::cout << "  inverse  - обращение матрицы (метод присоединённой матрицы)\n";     // работоспособен
        std::cout << "  inverse_lu - обращение матрицы (метод LU-разложения)\n";    // работоспособен
        std::cout << "  inverse_sqrt - обращение матрицы (метод квадратных корней)\n";  // работоспособен
        std::cout << "  qr_decomp - QR-разложение матрицы\n";
        std::cout << "  qr_alg    - QR-алгоритм (вывод матрицы с собственными значениями на диагонали)\n";
        std::cout << "  qr        - QR-алгоритм (вывод только собственных значений)\n";
        return 0;
    }

    std::string type   = argv[1];
    std::string method = argv[2];
    std::string file   = argv[3];

    try {
        if (type == "double") {
            math::matrix<double> M;
            M.load_from_file(file);

            if (method == "print") M.print();
            else if (method == "gauss")   { M.gauss();         M.print(); }
            else if (method == "jordan")  { M.jordan_gauss();  M.print(); }
            else if (method == "lu")      { M.lu_solve();      M.print(); }
            else if (method == "tridiag") { M.tridiagonal_solve(); M.print(); }
            else if (method == "sqrt")    { M.sqrt_method();   M.print(); }
            else if (method == "iter")    { M.iter_method();   M.print(); }
            else if (method == "seidel")  { M.seidel_method(); M.print(); }
            else if (method == "inverse") { auto inv = M.inverse_by_adjoint(); inv.print(); }
            else if (method == "inverse_lu") { auto inv = M.inverse_by_lu(); inv.print(); }
            else if (method == "inverse_sqrt") { auto inv = M.inverse_by_sqrt(); inv.print(); }
            else if (method == "qr") {
                auto eigvals = M.qr_eigenvalues();
                std::cout << "Eigenvalues:\n";
                for (size_t i = 0; i < eigvals.size(); i++) {
                    std::cout << "  λ" << i+1 << " = " << eigvals[i] << "\n";
                }
            }
            else if (method == "qr_decomp") {
                math::matrix<double> Q(M.size_x(), M.size_y());
                math::matrix<double> R(M.size_x(), M.size_y());
                M.qr_decomposition(Q, R);
                std::cout << "Q matrix:\n";
                Q.print();
                std::cout << "R matrix:\n";
                R.print();
            }
            else if (method == "qr_alg") {
                M.qr_algorithm();
                std::cout << "Matrix after QR algorithm (eigenvalues on diagonal):\n";
                M.print();
            }
            else std::cerr << "Unknown method\n";
        } else if (type == "float") {
            math::matrix<float> M;
            M.load_from_file(file);

            if (method == "print") M.print();
            else if (method == "gauss")   { M.gauss();         M.print(); }
            else if (method == "jordan")  { M.jordan_gauss();  M.print(); }
            else if (method == "lu")      { M.lu_solve();      M.print(); }
            else if (method == "tridiag") { M.tridiagonal_solve(); M.print(); }
            else if (method == "sqrt")    { M.sqrt_method();   M.print(); }
            else if (method == "iter")    { M.iter_method();   M.print(); }
            else if (method == "seidel")  { M.seidel_method(); M.print(); }
            else if (method == "inverse") { auto inv = M.inverse_by_adjoint(); inv.print(); }
            else if (method == "inverse_lu") { auto inv = M.inverse_by_lu(); inv.print(); }
            else if (method == "inverse_sqrt") { auto inv = M.inverse_by_sqrt(); inv.print(); }
            else if (method == "qr") {
                auto eigvals = M.qr_eigenvalues();
                std::cout << "Eigenvalues:\n";
                for (size_t i = 0; i < eigvals.size(); i++) {
                    std::cout << "  λ" << i+1 << " = " << eigvals[i] << "\n";
                }
            }
            else if (method == "qr_decomp") {
                math::matrix<float> Q(M.size_x(), M.size_y());
                math::matrix<float> R(M.size_x(), M.size_y());
                M.qr_decomposition(Q, R);
                std::cout << "Q matrix:\n";
                Q.print();
                std::cout << "R matrix:\n";
                R.print();
            }
            else if (method == "qr_alg") {
                M.qr_algorithm();
                std::cout << "Matrix after QR algorithm (eigenvalues on diagonal):\n";
                M.print();
            }
            else std::cerr << "Unknown method\n";
        } else {
            std::cerr << "Unsupported type (must be float or double)\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
    }

    return 0;
}
