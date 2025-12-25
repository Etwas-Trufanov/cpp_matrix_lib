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


        matrix<T> operator*(const matrix<T>& other) const {
            if (sy != other.sx)
                throw std::runtime_error("matrix multiply: size mismatch");

            matrix<T> result(sx, other.sy);

            for (std::size_t i = 0; i < sx; i++) {
                for (std::size_t j = 0; j < other.sy; j++) {
                    T sum = 0;
                    for (std::size_t k = 0; k < sy; k++) {
                        sum += m[i][k] * other.m[k][j];
                    }
                    result.m[i][j] = sum;
                }
            }
            return result;
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

        // ======================================================
        //  Транспонирование матрицы
        // ======================================================
        matrix<T> transpose() const {
            matrix<T> result(sy, sx); // Меняем размеры местами

            for (std::size_t x = 0; x < sx; x++)
                for (std::size_t y = 0; y < sy; y++)
                    result.m[y][x] = m[x][y];

            return result;
        }

        // =======================================
        // Операции со строками
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
        // ======================================================
        //  Гаусс — с прямым и обратным ходом (полное решение)
        // ======================================================
        void gauss() {
            std::size_t n = sy;
            if (sx != n + 1)
                throw std::runtime_error("gauss: matrix must be n×(n+1)");

            const T eps = 1e-12;

            // Прямой ход: приведение к верхнетреугольному виду
            for (std::size_t k = 0; k < n; k++) {
                // Поиск ненулевого элемента в столбце k, начиная со строки k
                if (std::abs(m[k][k]) < eps) {
                    bool found = false;
                    for (std::size_t r = k + 1; r < n; r++) {
                        if (std::abs(m[k][r]) > eps) {
                            swap_rows(k, r);
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        throw std::runtime_error("gauss: matrix is singular");
                    }
                }

                // Нормировка строки k
                T pivot = m[k][k];
                for (std::size_t j = k; j < sx; j++) {
                    m[j][k] /= pivot;
                }

                // Обнуление элементов ниже диагонали в столбце k
                for (std::size_t i = k + 1; i < n; i++) {
                    T factor = m[k][i];
                    for (std::size_t j = k; j < sx; j++) {
                        m[j][i] -= m[j][k] * factor;
                    }
                }
            }

            // Обратный ход
            std::vector<T> X(n);
            for (int i = n - 1; i >= 0; i--) {
                X[i] = m[n][i];  // Последний столбец (свободные члены)
                for (std::size_t j = i + 1; j < n; j++) {
                    X[i] -= m[j][i] * X[j];
                }
            }

            // Запись решения в последний столбец матрицы
            for (std::size_t i = 0; i < n; i++) {
                m[n][i] = X[i];
            }

            // оставим только диагональные единицы и решение в последнем столбце
            for (std::size_t i = 0; i < n; i++) {
                for (std::size_t j = i + 1; j < n; j++) {
                    m[j][i] = 0;
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
        //  Метод прогонки для трёхдиагональной системы
        //  Матрица должна быть n x (n+1)
        // ======================================================
        void tridiagonal_solve() {
            std::size_t n = sy;
            if (sx != n + 1)
                throw std::runtime_error("tridiagonal_solve: matrix must be n×(n+1)");

            // Векторы коэффициентов
            std::vector<T> a(n);  // нижняя диагональ (l) - индексы 1..n-1
            std::vector<T> b(n);  // главная диагональ (c) - индексы 0..n-1
            std::vector<T> c(n);  // верхняя диагональ (r) - индексы 0..n-2
            std::vector<T> f(n);  // правая часть - индексы 0..n-1

            // Извлечение коэффициентов из матрицы
            for (std::size_t i = 0; i < n; i++) {
                // Главная диагональ
                b[i] = m[i][i];

                // Правая часть (последний столбец)
                f[i] = m[n][i];

                // Верхняя диагональ (r)
                if (i < n - 1) {
                    c[i] = m[i+1][i];
                }

                // Нижняя диагональ (l)
                if (i > 0) {
                    a[i] = m[i-1][i];
                }
            }

            // Прямой ход метода прогонки
            std::vector<T> alpha(n), beta(n);

            // Первый шаг
            alpha[0] = -c[0] / b[0];
            beta[0] = f[0] / b[0];

            // Промежуточные шаги
            for (std::size_t i = 1; i < n - 1; i++) {
                T denominator = b[i] + a[i] * alpha[i-1];
                alpha[i] = -c[i] / denominator;
                beta[i] = (f[i] - a[i] * beta[i-1]) / denominator;
            }

            // Последний шаг (для нахождения beta[n-1])
            T denominator = b[n-1] + a[n-1] * alpha[n-2];
            beta[n-1] = (f[n-1] - a[n-1] * beta[n-2]) / denominator;

            // Обратный ход
            std::vector<T> X(n);
            X[n-1] = beta[n-1];

            for (int i = n - 2; i >= 0; i--) {
                X[i] = alpha[i] * X[i+1] + beta[i];
            }

            // Запись решения в последний столбец матрицы
            for (std::size_t i = 0; i < n; i++) {
                m[n][i] = X[i];
            }

            // Для наглядности: обнулим ненужные элементы и оставим только решение
            for (std::size_t y = 0; y < n; y++) {
                for (std::size_t x = 0; x < n; x++) {
                    m[x][y] = 0;
                }
                m[y][y] = 1;  // Единицы на диагонали для красоты
            }
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

        // ======================================================
        //  Обращение матрицы методом присоединённой матрицы
        //  (Точная реализация по методу из документа)
        // ======================================================
        matrix<T> inverse_by_adjoint() {
            if (sx != sy)
                throw std::runtime_error("inverse_by_adjoint: matrix must be square");

            std::size_t n = sx;
            const T eps = 1e-12;

            // Создаём расширенную матрицу [A|E]
            matrix<T> aug(2 * n, n);

            // Заполняем левую часть исходной матрицей
            for (std::size_t x = 0; x < n; x++)
                for (std::size_t y = 0; y < n; y++)
                    aug.m[x][y] = m[x][y];

            // Заполняем правую часть единичной матрицей
            for (std::size_t x = n; x < 2 * n; x++)
                for (std::size_t y = 0; y < n; y++)
                    aug.m[x][y] = (x - n == y) ? static_cast<T>(1) : static_cast<T>(0);

            // Применяем метод Жордана-Гаусса с правилом прямоугольника
            for (std::size_t k = 0; k < n; k++) {
                // Если диагональный элемент близок к нулю
                if (std::abs(aug.m[k][k]) < eps) {
                    bool found = false;
                    // Ищем строку ниже с ненулевым элементом в том же столбце
                    for (std::size_t i = k + 1; i < n; i++) {
                        if (std::abs(aug.m[k][i]) > eps) {
                            // Прибавляем строку i к строке k
                            for (std::size_t j = 0; j < 2 * n; j++) {
                                aug.m[j][k] += aug.m[j][i];
                            }
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        throw std::runtime_error("inverse_by_adjoint: matrix is singular");
                    }
                }

                // 1. Пересчёт элементов вне ведущей строки и столбца по правилу прямоугольника
                for (std::size_t i = 0; i < n; i++) {
                    if (i == k) continue;  // Пропускаем ведущую строку

                    for (std::size_t j = k + 1; j < 2 * n; j++) {
                        // Правило прямоугольника: (i,j) = (i,j)*a_kk - (i,k)*a_kj / a_kk
                        aug.m[j][i] = (aug.m[j][i] * aug.m[k][k] - aug.m[j][k] * aug.m[k][i]) / aug.m[k][k];
                    }
                }

                // 2. Обнуление элементов ведущего столбца (кроме ведущего элемента)
                for (std::size_t i = 0; i < n; i++) {
                    if (i != k) {
                        aug.m[k][i] = 0;
                    }
                }

                // 3. Деление ведущей строки на ведущий элемент (по убыванию индекса)
                T pivot = aug.m[k][k];
                for (int j = 2 * n - 1; j >= static_cast<int>(k); j--) {
                    aug.m[j][k] = aug.m[j][k] / pivot;
                }
            }

            // Извлекаем обратную матрицу из правой части
            matrix<T> result(n, n);
            for (std::size_t x = 0; x < n; x++)
                for (std::size_t y = 0; y < n; y++)
                    result.m[x][y] = aug.m[n + x][y];

            return result;
        }

        // ======================================================
        //  Обращение матрицы методом LU-разложения (Холецкого)
        // ======================================================
        matrix<T> inverse_by_lu() {
            if (sx != sy)
                throw std::runtime_error("inverse_by_lu: matrix must be square");

            std::size_t n = sx;
            const T eps = 1e-12;

            // Создаём копию исходной матрицы для вычислений
            matrix<T> A_copy = *this;

            // Шаг 1: LU-разложение (A = L * U)
            // L - нижняя треугольная, U - верхняя треугольная с единицами на диагонали

            matrix<T> L(n, n);
            matrix<T> U(n, n);

            // Инициализация U: единичная диагональ
            for (std::size_t i = 0; i < n; i++) {
                U.m[i][i] = static_cast<T>(1);
            }

            // Вычисление первого столбца L
            for (std::size_t i = 0; i < n; i++) {
                L.m[i][0] = A_copy.m[i][0];
            }

            // Вычисление первой строки U
            for (std::size_t j = 1; j < n; j++) {
                if (std::abs(L.m[0][0]) < eps)
                    throw std::runtime_error("inverse_by_lu: matrix is singular");
                U.m[0][j] = A_copy.m[0][j] / L.m[0][0];
            }

            // Вычисление остальных элементов L и U
            for (std::size_t k = 1; k < n; k++) {
                // Вычисление столбца k матрицы L
                for (std::size_t i = k; i < n; i++) {
                    L.m[i][k] = A_copy.m[i][k];
                    for (std::size_t m = 0; m < k; m++) {
                        L.m[i][k] -= L.m[i][m] * U.m[m][k];
                    }
                }

                // Вычисление строки k матрицы U (только элементы справа от диагонали)
                if (std::abs(L.m[k][k]) < eps)
                    throw std::runtime_error("inverse_by_lu: matrix is singular");

                for (std::size_t j = k + 1; j < n; j++) {
                    U.m[k][j] = A_copy.m[k][j];
                    for (std::size_t m = 0; m < k; m++) {
                        U.m[k][j] -= L.m[k][m] * U.m[m][j];
                    }
                    U.m[k][j] /= L.m[k][k];
                }
            }

            // Шаг 2: Вычисление обратных матриц L⁻¹ (Y) и U⁻¹ (X)

            // L⁻¹ (обозначаем Y) - нижняя треугольная матрица
            matrix<T> Y(n, n);

            for (std::size_t i = 0; i < n; i++) {
                for (std::size_t j = 0; j < n; j++) {
                    if (j > i) {
                        Y.m[j][i] = static_cast<T>(0);
                    } else if (j == i) {
                        if (std::abs(L.m[i][i]) < eps)
                            throw std::runtime_error("inverse_by_lu: matrix is singular");
                        Y.m[j][i] = static_cast<T>(1) / L.m[i][i];
                    } else { // j < i
                        T sum = 0;
                        for (std::size_t m = j; m <= i - 1; m++) {
                            sum += L.m[i][m] * Y.m[j][m];
                        }
                        if (std::abs(L.m[i][i]) < eps)
                            throw std::runtime_error("inverse_by_lu: matrix is singular");
                        Y.m[j][i] = -sum / L.m[i][i];
                    }
                }
            }

            // U⁻¹ (обозначаем X) - верхняя треугольная матрица
            matrix<T> X(n, n);

            // Вычисляем в обратном порядке (от n-1 до 0)
            for (int i = n - 1; i >= 0; i--) {
                for (int j = n - 1; j >= 0; j--) {
                    if (j < i) {
                        X.m[j][i] = static_cast<T>(0);
                    } else if (j == i) {
                        X.m[j][i] = static_cast<T>(1);
                    } else { // j > i
                        T sum = 0;
                        for (int m = i + 1; m <= j; m++) {
                            sum += U.m[i][m] * X.m[j][m];
                        }
                        X.m[j][i] = -sum;
                    }
                }
            }

            // Шаг 3: Вычисление A⁻¹ = U⁻¹ * L⁻¹ = X * Y
            matrix<T> result = X * Y;

            return result;
        }

        // ======================================================
        //  Обращение матрицы методом квадратных корней (Холецкого)
        //  Матрица должна быть симметрической и положительно определённой
        // ======================================================
        matrix<T> inverse_by_sqrt() {
            if (sx != sy)
                throw std::runtime_error("inverse_by_sqrt: matrix must be square");

            std::size_t n = sx;
            const T eps = 1e-12;

            // Проверка симметричности матрицы (опционально)
            for (std::size_t i = 0; i < n; i++) {
                for (std::size_t j = i + 1; j < n; j++) {
                    if (std::abs(m[i][j] - m[j][i]) > eps) {
                        throw std::runtime_error("inverse_by_sqrt: matrix must be symmetric");
                    }
                }
            }

            // Шаг 1: Разложение A = L * L^T (метод Холецкого)
            matrix<T> L(n, n);

            // Вычисление первого столбца L
            if (m[0][0] <= eps)
                throw std::runtime_error("inverse_by_sqrt: matrix is not positive definite");

            L.m[0][0] = std::sqrt(m[0][0]);

            for (std::size_t i = 1; i < n; i++) {
                L.m[i][0] = m[i][0] / L.m[0][0];
            }

            // Вычисление остальных элементов L
            for (std::size_t k = 1; k < n; k++) {
                // Диагональный элемент l[k][k]
                T sum_diag = 0;
                for (std::size_t m2 = 0; m2 < k; m2++) {
                    sum_diag += L.m[k][m2] * L.m[k][m2];
                }

                T diag = m[k][k] - sum_diag;
                if (diag <= eps)
                    throw std::runtime_error("inverse_by_sqrt: matrix is not positive definite");

                L.m[k][k] = std::sqrt(diag);

                // Элементы ниже диагонали l[i][k] для i > k
                for (std::size_t i = k + 1; i < n; i++) {
                    T sum = 0;
                    for (std::size_t m2 = 0; m2 < k; m2++) {
                        sum += L.m[i][m2] * L.m[k][m2];
                    }

                    L.m[i][k] = (m[i][k] - sum) / L.m[k][k];
                }
            }

            // Шаг 2: Вычисление обратной матрицы L⁻¹ (обозначаем Y)
            matrix<T> Y(n, n);

            for (std::size_t i = 0; i < n; i++) {
                for (std::size_t j = 0; j < n; j++) {
                    if (j > i) {
                        // Выше диагонали - нули
                        Y.m[j][i] = static_cast<T>(0);
                    }
                    else if (j == i) {
                        // Диагональ: y[i][i] = 1 / l[i][i]
                        if (std::abs(L.m[i][i]) < eps)
                            throw std::runtime_error("inverse_by_sqrt: division by zero");
                        Y.m[j][i] = static_cast<T>(1) / L.m[i][i];
                    }
                    else { // j < i
                        // Ниже диагонали: y[i][j] = -1/l[i][i] * sum(l[i][m] * y[m][j])
                        T sum = 0;
                        for (std::size_t m = j; m <= i - 1; m++) {
                            sum += L.m[i][m] * Y.m[j][m];
                        }

                        if (std::abs(L.m[i][i]) < eps)
                            throw std::runtime_error("inverse_by_sqrt: division by zero");

                        Y.m[j][i] = -sum / L.m[i][i];
                    }
                }
            }

            // Шаг 3: Вычисление A⁻¹ = Y^T * Y
            matrix<T> result = Y.transpose() * Y;

            return result;
        }

        // ======================================================
        // QR-разложение (модифицированный Грама-Шмидта)
        // ======================================================
        void qr_decomposition(matrix<T>& Q, matrix<T>& R) const {
            if (sx != sy)
                throw std::runtime_error("qr_decomposition: matrix must be square");

            std::size_t n = sx;
            Q = *this;
            R = matrix<T>(n, n);

            const T eps = 1e-12;

            for (std::size_t j = 0; j < n; j++) {
                // Норма столбца j
                T norm = 0;
                for (std::size_t i = 0; i < n; i++)
                    norm += Q.m[i][j] * Q.m[i][j];

                if (norm < eps) {
                    R.m[j][j] = 0;
                    for (std::size_t i = 0; i < n; i++)
                        Q.m[i][j] = 0;
                    Q.m[j][j] = 1;
                    continue;
                }

                norm = std::sqrt(norm);
                R.m[j][j] = norm;

                for (std::size_t i = 0; i < n; i++)
                    Q.m[i][j] /= norm;


                // Ортогонализация остальных столбцов
                for (std::size_t k = j + 1; k < n; k++) {
                    T dot = 0;
                    for (std::size_t i = 0; i < n; i++) {
                        dot += Q.m[i][j] * Q.m[i][k];
                    }

                    R.m[j][k] = dot;

                    for (std::size_t i = 0; i < n; i++) {
                        Q.m[i][k] -= Q.m[i][j] * dot;
                    }
                }
            }
        }

        // ======================================================
        // QR-алгоритм для поиска собственных значений
        // ======================================================
        std::vector<T> qr_method(int iterations = 3000, T tol = 1e-10) const {
            if (sx != sy)
                throw std::runtime_error("qr_method: matrix must be square");

            std::size_t n = sx;
            matrix<T> A = *this;
            matrix<T> Q(n, n), R(n, n);

            for (int iter = 0; iter < iterations; iter++) {
                // Сдвиг — последний диагональный элемент
                T mu = A.m[n-1][n-1];

                // A - mu*I
                matrix<T> shifted = A;
                for (std::size_t i = 0; i < n; i++)
                    shifted.m[i][i] -= mu;

                // QR-разложение
                shifted.qr_decomposition(Q, R);

                // A = R*Q + mu*I
                A = R * Q;
                for (std::size_t i = 0; i < n; i++)
                    A.m[i][i] += mu;

                // Проверка сходимости (внедиагональные элементы)
                T off = 0;
                for (std::size_t i = 1; i < n; i++)
                    off += std::abs(A.m[i][i-1]);

                if (off < tol)
                    break;
            }

            std::vector<T> eigenvalues(n);
            for (std::size_t i = 0; i < n; i++)
                eigenvalues[i] = A.m[i][i];

            return eigenvalues;
        }


    };

}
