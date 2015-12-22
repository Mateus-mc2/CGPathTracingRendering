#include <Eigen\Dense>
#include <opencv2\core\core.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2\highgui\highgui.hpp>

#include <chrono>
#include <iostream>
#include <random>

#include "quadric.h"

using Eigen::MatrixXi;
using Eigen::VectorXd;

int main() {
  MatrixXi I = MatrixXi::Identity(480, 480);
  I = 255*I;
  cv::Mat input(I.rows(), I.cols(), CV_8UC1, I.data());

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::normal_distribution<double> distribution;
  VectorXd v(10); 
  
  for (int i = 0; i < 5; ++i) {
    v(i) = distribution(generator);
  }

  try {
    util::Material mat(0.5, 0.5, 0.5, 0.5, 0.5, 1.2, 10);
    util::Quadric random_quadric(v, mat);

    std::cout << v << std::endl << "------" << std::endl;
    std::cout << random_quadric.coefficients() << std::endl;
  } catch (util::InvalidCoefficientsVectorException &e) {
    std::cout << e.what() << std::endl;
  } catch (util::InvalidMaterialCoefficientsException &e) {
    std::cout << e.what() << std::endl;
  }

  cv::imshow("Teste", input);
  cv::waitKey(0);
  return 0;
}
