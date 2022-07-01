// eigen useful commands

#include <iostream>    // for std::cout, std::endl
#include <Eigen/Dense> // for Eigen::

int main(){
    // 3*3 matrix A, Bと3*1 vector cの宣言と初期化
    Eigen::Matrix3d A;
    Eigen::Matrix3d B;
    Eigen::Vector3d c;
    Eigen::Vector3d d;
    Eigen::Matrix3d id = Eigen::Matrix3d::Identity(3, 3); // 単位行列で初期化
    A << 2,3,1, 4,-5,-1, 6,-7,1;
    B << 1,3,4, -5,-6,3, 6,-1,0;
    c << 7, 3, 5;
    d << -3, 1, 8;

    // matrix A, vecotr c, matrix idの出力
    std::cout << A << std::endl;
    std::cout << c << std::endl;
    std::cout << id << std::endl;

    // 要素へのアクセス(0行0列から始まることに注意)
    std::cout << A(2, 2) << std::endl;        // matrix Aの(2,2)成分
    std::cout << A.row(0) << std::endl;       // matrix Aの0行目
    std::cout << B.col(1) << std::endl;       // matrix Bの1列目
    std::cout << c.segment(0,2) << std::endl; //vector cの0要素目から1要素まで
    
    // 演算
    std::cout << A + B << std::endl; // matrix Aとmatrix Bの和  
    std::cout << A - B << std::endl; // matrix Aとmatrix Bの差
    std::cout << 2 * A << std::endl; // matrix Aの各要素を2倍
    std::cout << A * B << std::endl; // matrix Aとmatrix Bの積

    // 転置, 随伴, 逆行列
    std::cout << A.transpose() << std::endl; // matrix Aの転置行列
    std::cout << B.adjoint() << std::endl;   // matrix Bの随伴行列
    std::cout << A.inverse() << std::endl;   // matrix Aの逆行列

    // vectorの内積と外積
    std::cout << c.dot(d) << std::endl;   // vector cとvector dの内積
    std::cout << c.cross(d) << std::endl; //vector cとvector dの外積

    // Ax = cを解く
    Eigen::ColPivHouseholderQR <Eigen::Matrix3d> dec(A);
    Eigen::Vector3d x = dec.solve(c);
    std::cout << x << std::endl;

    // LU分解でAx = cを解く
    Eigen::FullPivLU <Eigen::Matrix3d> lu(A); //Eigen::PartialPivLUでも良い
    Eigen::Vector3d y = lu.solve(c);
    std::cout << y << std::endl;

    // 固有値分解
    Eigen::EigenSolver <Eigen::Matrix3d> eigensolver(A);
    std::cout << eigensolver.eigenvalues() << std::endl;  // 固有値
    std::cout << eigensolver.eigenvectors() << std::endl; // 固有ベクトル
}
