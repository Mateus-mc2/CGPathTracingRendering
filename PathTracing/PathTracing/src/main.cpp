#include <Eigen\Dense>
#include <opencv2\core\core.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2\highgui\highgui.hpp>

#include <iostream>

using Eigen::MatrixXi;

int main() {
  MatrixXi I = MatrixXi::Identity(480, 480);
  I = 255*I;
  cv::Mat input(I.rows(), I.cols(), CV_8UC1, I.data());

  std::cout << "Bunda" << std::endl;
  cv::imshow("Teste", input);
  cv::waitKey(0);
  return 0;
}
