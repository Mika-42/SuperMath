#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <stdexcept>
#include <string>
#include <functional>
#include <complex>

namespace smath
{
    using namespace std::literals;
    using namespace std::complex_literals;

    typedef std::complex<double> t_cmplx;
    typedef std::vector<std::vector<t_cmplx>> t_2dVector;

    class SMathException : public std::exception
    {
    public:
        SMathException(const std::string &msg) : m_msg{msg}
        {
        }

        const char *what() const noexcept override
        {
            return m_msg.c_str();
        }

    private:
        const std::string m_msg;
    };

    class SMathDimensionsError : public SMathException
    {
    public:
        SMathDimensionsError(void)
            : SMathException("Error::Dimensions: Wrong matrix dimensions."s)
        {
        }
    };

    class SMathMultiplicityError : public SMathException
    {
    public:
        SMathMultiplicityError(void)
            : SMathException("Error::Multiplicity: columns mismatch with rows."s)
        {
        }
    };

    class SMathSquaredDimensionError : public SMathException
    {
    public:
        SMathSquaredDimensionError(void)
            : SMathException("Error::Power::Dimension: Only squared matrices allowed."s)
        {
        }
    };

    class SMathPowerExponentError : public SMathException
    {
    public:
        SMathPowerExponentError(void)
            : SMathException("Error::Power::Exponent: The exponant must be greater than 0."s)
        {
        }
    };

    class SMathSingularMatrixError : public SMathException
    {
    public:
        SMathSingularMatrixError(void)
            : SMathException("Error::Inverse: A singular matrix cannot be inverted."s)
        {
        }
    };

    class Matrix
    {
    public:
        Matrix(std::size_t rows, std::size_t cols, t_cmplx scalar = 0.0 + 0i) : m_rows{rows},
                                                                                m_cols{cols},
                                                                                m_array(m_rows, std::vector<t_cmplx>(m_cols, scalar))
        {
            if (m_rows < 1 or m_cols < 1)
            {
                throw SMathDimensionsError();
            }
        }

        constexpr t_cmplx &operator[](std::size_t row, std::size_t col)
        {
            return m_array[row][col];
        }

        Matrix &operator=(const Matrix &other) = default;
        
        void set(const t_2dVector &matrix) //? mabe use another semantic
        {
            if (not areVector2dDimensionsEqual(matrix, m_array))
            {
                throw SMathDimensionsError();
            }

            m_array = matrix;
        }

        void print(void) //? mabe use another semantic
        {
            for (auto &k : m_array)
            {
                std::cout << '[';

                for (bool isFirst{true}; auto &j : k)
                {
                    std::cout << (isFirst ? isFirst = false, "" : ", ") << j;
                }

                std::cout << ']' << std::endl;
            }
            std::cout << std::endl;
        }

        t_cmplx det(void)
        {
            return det(*this);
        }

        Matrix T(void)
        {
            Matrix ret(m_cols, m_rows, 0);

            for (std::size_t row{0}; row < m_rows; ++row)
            {
                for (std::size_t col{0}; col < m_cols; ++col)
                {
                    ret[col, row] = (*this)[row, col];
                }
            }

            return ret;
        }

        Matrix diag(void)
        {
            // todo
        }

        Matrix inv(void)
        {
            if (m_cols != m_rows)
            {
                throw SMathSquaredDimensionError();
            }

            if (det() == 0.0)
            {
                throw SMathSingularMatrixError();
            }

            return (1.0 + 0i / det()) * com().T();
        }

        Matrix com(void)
        {
            if (m_cols != m_rows)
            {
                throw SMathSquaredDimensionError();
            }

            if (m_cols == 1 and m_rows == 1)
            {
                return identity(1);
            }

            if (m_cols == 2 and m_rows == 2)
            {
                Matrix comatrix(m_rows, m_cols);
                comatrix.set(
                    {{(*this)[1, 1], -(*this)[1, 0]},
                     {-(*this)[0, 1], (*this)[0, 0]}});

                return comatrix;
            }

            Matrix comatrix(m_rows, m_cols);
            auto sign{1.0 + 0i};

            for (std::size_t row{0}; row < m_rows; ++row)
            {
                for (std::size_t col{0}; col < m_cols; ++col, sign = -sign)
                {
                    auto subMatrix{decreaseDegrees(*this, row, col)};

                    comatrix[row, col] = sign * det(subMatrix);
                }
            }
            return comatrix;
        }

        Matrix pow(unsigned int exponent)
        {
            if (m_cols != m_rows)
            {
                throw SMathSquaredDimensionError();
            }

            if (exponent <= 0)
            {
                throw SMathPowerExponentError();
            }

            Matrix ret = identity(m_rows);

            for (std::size_t k{0}; k < exponent; ++k)
            {
                ret = ret * (*this);
            }

            return ret;
        }

        static Matrix identity(std::size_t order)
        {
            if (order < 1)
            {
                throw SMathDimensionsError();
            }

            Matrix ret(order, order, 0);

            for (std::size_t k{0}; k < order; ++k)
            {
                ret[k, k] = 1.0+0i;
            }
            return ret;
        }

        friend Matrix operator*(const Matrix &lhs, const Matrix &rhs)
        {
            if (lhs.m_cols != rhs.m_rows)
            {
                throw SMathMultiplicityError();
            }

            Matrix ret(lhs.m_rows, rhs.m_cols, 0);

            for (std::size_t row{0}; row < ret.m_rows; ++row)
            {
                for (std::size_t col{0}; col < ret.m_cols; ++col)
                {
                    for (std::size_t k{0}; k < lhs.m_cols; ++k)
                    {
                        ret[row, col] += lhs.m_array[row][k] * rhs.m_array[k][col];
                    }
                }
            }

            return ret;
        }

        friend Matrix operator*(const t_cmplx &scalar, const Matrix &matrix)
        {
            return matrix * scalar;
        }

        friend Matrix operator*(const Matrix &matrix, const t_cmplx &scalar)
        {
            Matrix ret(matrix.m_rows, matrix.m_cols, scalar);

            return superOperator(matrix, ret, [](t_cmplx lhs, t_cmplx rhs)
                                 { return lhs * rhs; });
        }

        friend Matrix operator+(const Matrix &lhs, const Matrix &rhs)
        {
            return superOperator(lhs, rhs, [](t_cmplx lhs, t_cmplx rhs)
                                 { return lhs + rhs; });
        }

        friend Matrix operator-(const Matrix &lhs, const Matrix &rhs)
        {
            return superOperator(lhs, rhs, [](t_cmplx lhs, t_cmplx rhs)
                                 { return lhs - rhs; });
        }

    private:
        static Matrix superOperator(const Matrix &lhs, const Matrix &rhs, std::function<t_cmplx(t_cmplx, t_cmplx)> op)
        {
            if (not(lhs.m_cols == rhs.m_cols and lhs.m_rows == rhs.m_rows))
            {
                throw SMathDimensionsError();
            }

            Matrix ret(lhs.m_rows, lhs.m_cols);

            for (std::size_t k{0}; k < lhs.m_rows; ++k)
            {
                for (std::size_t j{0}; j < lhs.m_cols; ++j)
                {
                    ret[k, j] = op(lhs.m_array[k][j], rhs.m_array[k][j]);
                }
            }

            return ret;
        }

        static Matrix decreaseDegrees(const Matrix &matrix, std::size_t row, std::size_t col)
        {
            Matrix subMatrix = matrix;
            subMatrix.m_array.erase(subMatrix.m_array.begin() + row);
            subMatrix.m_rows = subMatrix.m_cols = matrix.m_cols - 1;
            for (auto &k : subMatrix.m_array)
            {
                k.erase(k.begin() + col);
            }
            return subMatrix;
        }

        static t_cmplx det(const Matrix &matrix)
        {
            if (matrix.m_cols != matrix.m_rows)
            {
                throw SMathSquaredDimensionError();
            }

            if (matrix.m_cols == 1 and matrix.m_rows == 1) // scalar
            {
                return matrix.m_array[0][0];
            }

            if (matrix.m_cols == 2 and matrix.m_rows == 2) // 2×2
            {
                auto &coef{matrix.m_array};
                return coef[0][0] * coef[1][1] - coef[1][0] * coef[0][1];
            }

            t_cmplx ret{0.0};

            for (std::size_t k{0}; k < matrix.m_cols; ++k)
            {
                auto subMatrix{decreaseDegrees(matrix, 0, k)};

                if (subMatrix.m_cols == 2)
                {
                    auto &coef{subMatrix.m_array};
                    t_cmplx subDet{coef[0][0] * coef[1][1] - coef[1][0] * coef[0][1]};
                    ret += ((k % 2) ? -1.0 + 0i : 1.0 + 0i) * matrix.m_array[0][k] * subDet;
                }
                else
                {
                    ret += ((k % 2) ? -1.0 + 0i : 1.0 + 0i) * matrix.m_array[0][k] * det(subMatrix);
                }
            }
            return ret;
        }

        static bool areVector2dDimensionsEqual(const t_2dVector &lhs, const t_2dVector &rhs)
        {
            if (lhs.size() != rhs.size())
            {
                return false;
            }

            for (std::size_t k{0}; k < lhs.size(); ++k)
            {
                if (lhs[k].size() != rhs[k].size())
                {
                    return false;
                }
            }

            return true;
        }

    private:
        std::size_t m_rows;
        std::size_t m_cols;
        t_2dVector m_array;
    };
} // namespace smath

#endif // MATRIX_H