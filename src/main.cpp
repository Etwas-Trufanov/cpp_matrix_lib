#include <iostream>
#include "matrix.hpp"

using namespace math;

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
        std::cout << "  tridiag  - метод прогонки (трёхдиагональные СЛАУ)\n";
        std::cout << "  sqrt     - метод квадратных корней\n";
        std::cout << "  iter     - метод простых итераций\n";
        std::cout << "  seidel   - метод Зейделя\n";
        return 0;
    }

    std::string type   = argv[1];
    std::string method = argv[2];
    std::string file   = argv[3];

    try
    {
        if (type == "double")
        {
            matrix<double> M;
            M.load_from_file(file);

            if (method == "print") M.print();
            else if (method == "gauss")   { M.gauss();         M.print(); }
            else if (method == "jordan")  { M.jordan_gauss();  M.print(); }
            else if (method == "lu")      { M.lu_solve();      M.print(); }
            else if (method == "tridiag") { M.tridiagonal_solve(); M.print(); }
            else if (method == "sqrt")    { M.sqrt_method();   M.print(); }
            else if (method == "iter")    { M.iter_method();   M.print(); }
            else if (method == "seidel")  { M.seidel_method(); M.print(); }
            else std::cerr << "Unknown method\n";
        }
        else if (type == "float")
        {
            matrix<float> M;
            M.load_from_file(file);

            if (method == "print") M.print();
            else if (method == "gauss")   { M.gauss();         M.print(); }
            else if (method == "jordan")  { M.jordan_gauss();  M.print(); }
            else if (method == "lu")      { M.lu_solve();      M.print(); }
            else if (method == "tridiag") { M.tridiagonal_solve(); M.print(); }
            else if (method == "sqrt")    { M.sqrt_method();   M.print(); }
            else if (method == "iter")    { M.iter_method();   M.print(); }
            else if (method == "seidel")  { M.seidel_method(); M.print(); }
            else std::cerr << "Unknown method\n";
        }
        else {
            std::cerr << "Unsupported type (must be float or double)\n";
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "ERROR: " << e.what() << "\n";
    }

    return 0;
}
