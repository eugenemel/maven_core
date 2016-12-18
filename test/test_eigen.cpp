#include <Eigen/Core>
//#include <Eigen/LU>
#include <Eigen/SVD>
#include <iostream>
//#include <Eigen/Array>

using namespace Eigen;
using namespace std;

int main(int, char *[])
{

int plex = 10;
MatrixXf TMT(plex,7);
TMT<<                   //p=0,j==3
0,	0,	0.0,	0, 4.69,	0.0,	0,
0,	0,	0.4,	0, 6.5,	0.0,	0,
0,	0,	0.2,	0, 4.6,	0.3,	0,
0,	0,	0.9,	0, 4.7,	0.2,	0,
0,	0.1,	0.53,	0, 2.59,	0.0,	0,
0,	0,	0.73,	0, 2.49,	0.0,	0,
0,	0,	1.3,	0, 2.5,	0.0,	0,
0,	0,	1.2,	0, 2.8,	2.7,	0,
0,	0.1,	2.9,	0, 2.9,	0.0,	0,
0,	0,	2.36,	0, 1.43,	0.0,	0;

TMT /= 100;
VectorXf diag = TMT.transpose().colwise().sum();
diag  =  VectorXf::Ones(plex)-diag;
MatrixXf IMP= diag.asDiagonal();

for(int p=0; p<plex; p++) {
    int j=0;
    for(int c=-4; c<3; c++) {
        int aC = c+p+1;
        if (aC > 0 and aC < plex) IMP(p,aC) += TMT(p,j);
        j++;
    }
}

//x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
cerr << IMP << endl;

MatrixXf A = IMP.transpose();
VectorXf b = VectorXf::Ones(plex);
VectorXf x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);


 /*
  Matrix3f m3;
  m3 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  Matrix4f m4 = Matrix4f::Identity();
  Vector4i v4(1, 2, 3, 4);

  std::cout << "m3\n" << m3 << "\nm4:\n" << m4 << "\nv4:\n" << v4 << std::endl;

  MatrixXf A(3,4); A << 1,2,0,1, 1,1,1,-1, 3,1,5,-7;
  VectorXf b(3);   b << 7,3,1;
  VectorXf x(3);

  std::cout << std::endl;
  std::cout << A << std::endl;
  std::cout << b << std::endl;
  std::cout << std::endl;

  x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
  std::cout << x;

  std::cout << std::endl;
  std::cout << A*x << endl;
  */

//  MatrixXf A = MatrixXf::Random(5,5);
//   VectorXf b(5); 
  // for(int i=0; i < b.size(); i++ ) std::cout << b(i) << std::endl;
  // std::cout << b.sum();
//   VectorXf x;
//  A.svd().solve(b, &x);
//
  // std::cout << VectorXi::Random(2) << std::endl;

  //SVD<MatrixXf>svdOfA(A);
 // svdOfA.solve(b, &x);


}

