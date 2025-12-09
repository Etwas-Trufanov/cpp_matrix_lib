#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <type_traits>
#include <cmath>

namespace math {

    template<typename T>
    requires std::is_arithmetic_v<T>
    class matrix {
    private:
        std::vector<std::vector<T>> m;
        std::size_t sx{}, sy{};

    public:

        // ------------------------
        //  Конструктор
        // ------------------------
        matrix(std::size_t size_x = 0, std::size_t size_y = 0)
        : sx(size_x), sy(size_y) {
            m.assign(sx, std::vector<T>(sy, T{}));
        }

        // Доступ
        T& at(std::size_t x, std::size_t y) { return m.at(x).at(y); }
        const T& at(std::size_t x, std::size_t y) const { return m.at(x).at(y); }

        std::size_t size_x() const { return sx; }
        std::size_t size_y() const { return sy; }

        // ------------------------
        //  Загрузка из файла
        // ------------------------
        void load_from_file(const std::string& filename) {
            std::ifstream f(filename);
            if (!f.is_open())
                throw std::runtime_error("Cannot open file: " + filename);

            f >> sx >> sy;
            m.assign(sx, std::vector<T>(sy));

            for (std::size_t y = 0; y < sy; y++)
                for (std::size_t x = 0; x < sx; x++)
                    f >> m[x][y];
        }

        // ------------------------
        //  Прямая загрузка
        // ------------------------
        void load_direct() {
            for (std::size_t y = 0; y < sy; y++)
                for (std::size_t x = 0; x < sx; x++)
                    std::cin >> m[x][y];
        }

        // ------------------------
        // Печать
        // ------------------------
        void print() const {
            std::cout << "----------------------------------------\n";
            for (std::size_t y = 0; y < sy; y++) {
                for (std::size_t x = 0; x < sx; x++)
                    std::cout << std::setw(10) << m[x][y] << " ";
                std::cout << "\n";
            }
        }

        // ======================================================
        // Базовые операции
        // ======================================================

        matrix<T> operator+(const matrix<T>& b) const {
            if (sx != b.sx || sy != b.sy)
                throw std::runtime_error("Matrix sizes incompatible for +");

            matrix<T> r(sx, sy);

            for (std::size_t y = 0; y < sy; y++)
                for (std::size_t x = 0; x < sx; x++)
                    r.m[x][y] = m[x][y] + b.m[x][y];

            return r;
        }


        matrix<T> operator*(const matrix<T>& b) const {
            if (sx != b.sy)
                throw std::runtime_error("Matrix sizes incompatible for *");

            matrix<T> r(b.sx, sy);

            for (std::size_t y = 0; y < sy; y++) {
                for (std::size_t x = 0; x < b.sx; x++) {
                    T acc = 0;
                    for (std::size_t k = 0; k < sx; k++)
                        acc += m[k][y] * b.m[x][k];
                    r.m[x][y] = acc;
                }
            }
            return r;
        }

        void add_value(T v) {
            for (auto& col : m)
                for (auto& e : col)
                    e += v;
        }

        void scale(T v) {
            for (auto& col : m)
                for (auto& e : col)
                    e *= v;
        }

        void regularize(T v) {
            if (sx != sy)
                throw std::runtime_error("regularize: matrix must be square");
            for (std::size_t i = 0; i < sx; i++)
                m[i][i] += v;
        }

        // =======================================
        // Операции со строками (как в C-версии)
        // =======================================

        void swap_rows(std::size_t a, std::size_t b) {
            for (std::size_t x = 0; x < sx; x++)
                std::swap(m[x][a], m[x][b]);
        }

        void multiply_row(std::size_t r, T v) {
            for (std::size_t x = 0; x < sx; x++)
                m[x][r] *= v;
        }

        void row_add_mult(std::size_t fr, std::size_t to, T v) {
            for (std::size_t x = 0; x < sx; x++)
                m[x][to] += m[x][fr] * v;
        }

        // ------------------------
        // Подматрица
        // ------------------------
        matrix<T> submatrix(std::size_t fx, std::size_t fy, std::size_t tx, std::size_t ty) const {
            if (fx > tx || fy > ty || tx >= sx || ty >= sy)
                throw std::runtime_error("Invalid submatrix");

            std::size_t nx = tx - fx + 1;
            std::size_t ny = ty - fy + 1;

            matrix<T> r(nx, ny);
            for (std::size_t y = 0; y < ny; y++)
                for (std::size_t x = 0; x < nx; x++)
                    r.m[x][y] = m[x + fx][y + fy];

            return r;
        }

        // ======================================================
        //  Гаусс — приведение к верхней треугольной
        // ======================================================
        void gauss() {
            std::size_t n = sy;
            for (std::size_t k = 0; k < n; k++)
            {
                if (std::abs(m[k][k]) < 1e-12) {
                    for (std::size_t r = k + 1; r < n; r++)
                        if (std::abs(m[k][r]) > 1e-12)
                            swap_rows(k, r);
                }

                for (std::size_t y = k + 1; y < n; y++) {
                    T c = m[k][y] / m[k][k];
                    for (std::size_t x = k; x < sx; x++)
                        m[x][y] -= m[x][k] * c;
                }
            }
        }

        // ======================================================
        // Жордан — к диагональному виду
        // ======================================================
        void jordan_gauss() {
            std::size_t n = sy;

            for (std::size_t k = 0; k < n; k++) {
                if (std::abs(m[k][k]) < 1e-12) {
                    for (std::size_t r = k + 1; r < n; r++)
                        if (std::abs(m[k][r]) > 1e-12)
                            swap_rows(k, r);
                }

                T d = m[k][k];
                for (std::size_t x = k; x < sx; x++)
                    m[x][k] /= d;

                for (std::size_t y = 0; y < n; y++) {
                    if (y == k) continue;
                    T c = m[k][y];
                    for (std::size_t x = k; x < sx; x++)
                        m[x][y] -= m[x][k] * c;
                }
            }
        }

        // ======================================================
        // LU-разложение + решение Ax=b
        // Матрица должна быть n x (n+1)
        // ======================================================
        void lu_solve()
        {
            std::size_t n = sy;
            if (sx != n + 1)
                throw std::runtime_error("LU: matrix must be n×(n+1)");

            std::vector<std::vector<T>> L(n, std::vector<T>(n, 0));
            std::vector<std::vector<T>> U(n, std::vector<T>(n, 0));

            for (std::size_t i = 0; i < n; i++)
                U[i][i] = 1;

            for (std::size_t k = 0; k < n; k++) {
                for (std::size_t i = k; i < n; i++) {
                    T sum = 0;
                    for (std::size_t t = 0; t < k; t++)
                        sum += L[t][i] * U[k][t];
                    L[k][i] = m[k][i] - sum;
                }

                for (std::size_t j = k + 1; j < n; j++) {
                    T sum = 0;
                    for (std::size_t t = 0; t < k; t++)
                        sum += L[t][k] * U[j][t];
                    U[j][k] = (m[j][k] - sum) / L[k][k];
                }
            }

            // Решения
            std::vector<T> Y(n), X(n);

            for (std::size_t i = 0; i < n; i++) {
                T sum = 0;
                for (std::size_t t = 0; t < i; t++)
                    sum += L[t][i] * Y[t];
                Y[i] = (m[n][i] - sum) / L[i][i];
            }

            for (int i = n - 1; i >= 0; i--) {
                T sum = 0;
                for (std::size_t t = i + 1; t < n; t++)
                    sum += U[t][i] * X[t];
                X[i] = Y[i] - sum;
            }

            for (std::size_t i = 0; i < n; i++)
                m[n][i] = X[i];
        }
        // ======================================================
        // Метод квадратных корней (Холецкого)
        // Решает Ax = b, где матрица — n×(n+1)
        // ======================================================
        void sqrt_method() {
            std::size_t n = sy;
            if (sx != n + 1)
                throw std::runtime_error("sqrt_method: matrix must be n×(n+1)");

            // A — первые n столбцов
            // b — последний столбец

            // Матрица L
            std::vector<std::vector<T>> L(n, std::vector<T>(n, 0));

            // ===== Построение L =====
            for (std::size_t k = 0; k < n; k++) {
                // Диагональ: l[k][k] = sqrt( a[k][k] - sum(l[k][m]^2) )
                T sum_diag = 0;
                for (std::size_t m2 = 0; m2 < k; m2++)
                    sum_diag += L[k][m2] * L[k][m2];

                T diag = m[k][k] - sum_diag;
                if (diag <= 0)
                    throw std::runtime_error("sqrt_method: matrix is not positive definite");

                L[k][k] = std::sqrt(diag);

                // Вне диагонали: l[i][k] = (a[i][k] - sum(l[i][m]*l[k][m])) / l[k][k]
                for (std::size_t i = k + 1; i < n; i++) {
                    T sum = 0;
                    for (std::size_t m2 = 0; m2 < k; m2++)
                        sum += L[i][m2] * L[k][m2];

                    L[i][k] = (m[k][i] - sum) / L[k][k];
                }
            }

            // ============ Решаем LY = B ============
            std::vector<T> Y(n);
            for (std::size_t i = 0; i < n; i++) {
                T sum = 0;
                for (std::size_t m2 = 0; m2 < i; m2++)
                    sum += L[i][m2] * Y[m2];

                Y[i] = (m[n][i] - sum) / L[i][i];
            }

            // ============ Решаем Lᵀ X = Y ============
            std::vector<T> X(n);

            for (int i = n - 1; i >= 0; i--) {
                T sum = 0;

                for (std::size_t m2 = i + 1; m2 < n; m2++)
                    sum += L[m2][i] * X[m2];

                X[i] = (Y[i] - sum) / L[i][i];
            }

            // Записываем результат в последний столбец
            for (std::size_t i = 0; i < n; i++)
                m[n][i] = X[i];

            // Обнулим A и поставим единицы на диагонали
            for (std::size_t y = 0; y < n; y++)
                for (std::size_t x = 0; x < n; x++)
                    m[x][y] = (x == y ? 1 : 0);
        }
        // ======================================================
        // Метод простых итераций (Jacobi)
        // ======================================================
        void iter_method(T eps = 1e-6, int max_iter = 10000) {
            size_t n = sy;
            if (sx != n + 1)
                throw std::runtime_error("iter_method: matrix must be n×(n+1)");

            std::vector<T> X(n, 0);
            std::vector<T> XX(n, 0);

            for (int step = 0; step < max_iter; step++) {
                int good = 0;

                for (size_t i = 0; i < n; i++) {
                    XX[i] = m[n][i]; // B[i]
                    for (size_t j = 0; j < n; j++)
                        if (j != i)
                            XX[i] -= m[j][i] * X[j];

                    XX[i] /= m[i][i];
                }

                for (size_t i = 0; i < n; i++) {
                    if (std::fabs(XX[i] - X[i]) < eps)
                        good++;
                    X[i] = XX[i];
                }

                if (good == n) break;
            }

            for (size_t i = 0; i < n; i++)
                m[n][i] = X[i];
        }


        // ======================================================
        // Метод Зейделя (Gauss–Seidel)
        // ======================================================
        void seidel_method(T eps = 1e-6, int max_iter = 10000) {
            size_t n = sy;
            if (sx != n + 1)
                throw std::runtime_error("seidel_method: matrix must be n×(n+1)");

            std::vector<T> X(n, 0);
            std::vector<T> XX(n, 0);

            for (int step = 0; step < max_iter; step++) {
                int good = 0;

                for (size_t i = 0; i < n; i++) {
                    XX[i] = m[n][i]; // B[i]

                    for (size_t j = 0; j < n; j++) {
                        if (j < i) {
                            XX[i] -= m[j][i] * XX[j]; // новое значение
                        } else if (j > i) {
                            XX[i] -= m[j][i] * X[j]; // старое значение
                        }
                    }

                    XX[i] /= m[i][i];
                }

                for (size_t i = 0; i < n; i++) {
                    if (std::fabs(XX[i] - X[i]) < eps)
                        good++;
                    X[i] = XX[i];
                }

                if (good == n) break;
            }

            for (size_t i = 0; i < n; i++)
                m[n][i] = X[i];
        }

    };

}
