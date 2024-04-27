/**
 * @file sandbox.cpp
 * @author lost_in_nowhere
 * @brief sandbox to try the code
 * @version 0.1
 * @date 2024-04-23
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "matrix.hpp"

int main()
{

    smath::Matrix A(2, 2, 0);
    smath::Matrix B(3, 3, 0);
    smath::Matrix C(4, 4, 0);

    A.set({{-2.0/3.0, 1.0/6.0}, {1.0/2.0, 0}});
    B.set({{1, 2, 3}, {4, 5, 6},{7, 8, 9}});
    C.set({{1,42,4,6}, {1,3,5,7}, {6,1,2,3}, {9,8,7,6}});

    // smath::Matrix C = 3 * A; //! do not work

    // smath::Matrix C(0, 0);
    // A.pow(2).print();

    // smath::Matrix::identity(3).print();

    //(A * B).print();

    // A.det(); // 3
    //A.print();
    smath::Matrix::identity(2).inv().print();
    A.com();
    A.det();
    //A.diag();
    A.inv();
    A.pow(2);
    A.print();
    A.T();

    //std::cout << "A[0, 0] = " << A.det() << std::endl; // 0
    
    //B.print();
    //std::cout << "|B| = " << B.det() << std::endl << std::endl; // -12

    //C.print();
    //std::cout << "|C| = " << C.det() << std::endl << std::endl; // 0
    
    // smath::Matrix::A, 1, 0).print();

    return 0;
}