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
        Matrix(std::size_t rows, std::size_t cols, t_cmplx scalar = 0.0 + 0i);

        t_cmplx &operator[](std::size_t row, std::size_t col);
        
        Matrix &operator=(const Matrix &other) = default;
        
        void set(const t_2dVector &matrix);

        void print(void);

        t_cmplx det(void);

        Matrix T(void);

        Matrix diag(void);

        Matrix inv(void);

        Matrix com(void);

        Matrix pow(unsigned int exponent);

        static Matrix identity(std::size_t order);

        friend Matrix operator*(const Matrix &lhs, const Matrix &rhs);

        friend Matrix operator*(const t_cmplx &scalar, const Matrix &matrix);

        friend Matrix operator*(const Matrix &matrix, const t_cmplx &scalar);

        friend Matrix operator+(const Matrix &lhs, const Matrix &rhs);

        friend Matrix operator-(const Matrix &lhs, const Matrix &rhs);

        bool operator==(const Matrix &matrix);

        bool operator!=(const Matrix &matrix);

    private:
        static Matrix superOperator(const Matrix &lhs, const Matrix &rhs, std::function<t_cmplx(t_cmplx, t_cmplx)> op);

        static Matrix decreaseDegrees(const Matrix &matrix, std::size_t row, std::size_t col);

        static t_cmplx det(const Matrix &matrix);

    private:
        std::size_t m_rows;
        std::size_t m_cols;
        t_2dVector m_array;
    };
} // namespace smath

#endif // MATRIX_H