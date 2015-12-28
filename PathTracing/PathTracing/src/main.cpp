#include <Eigen\Dense>
#include <opencv2\core\core.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2\highgui\highgui.hpp>

#include <chrono>
#include <iostream>
#include <random>

#include "quadric.h"
#include "pnm_writer.h"

using Eigen::MatrixXi;
using Eigen::VectorXd;

int main() {
  MatrixXi I = MatrixXi::Identity(480, 480);
  I = 255*I;
  cv::Mat input(I.rows(), I.cols(), CV_8UC1, I.data());

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> distribution(0.0, 1000.0);
  VectorXd v(10); 
  
  for (int i = 0; i < 10; ++i) {
    v(i) = distribution(generator);
  }

  try {
    util::Material mat(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 10);
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

  double eps = 1.0e-08;

  for (int i = 0; i < 10000; ++i) {
    double x = distribution(generator);
    double y = distribution(generator);
    double z = distribution(generator);

    // Reflexive.
    assert(math::IsAlmostEqual(x, x, eps));
    assert(math::IsAlmostEqual(y, y, eps));
    assert(math::IsAlmostEqual(z, z, eps));

    // Symmetric.
    if (math::IsAlmostEqual(x, y, eps)) {
      assert(math::IsAlmostEqual(y, x, eps));
    }

    // Transitive.
    if (math::IsAlmostEqual(x, y, eps) && math::IsAlmostEqual(y, z, eps)) {
      assert(math::IsAlmostEqual(x, z, eps));
    }    
  }

  //## Testando o pnm_writer
  cv::Mat imagem = cv::imread("toreba.png");
  cv::imshow("Jubiloca", imagem);
  cv::waitKey(0);

  io::PNMWriter pnm_mgr("../../../data/output/");
  pnm_mgr.WritePNMFile(imagem);
  // Teste colocando diretorio absoluto e nome do arquivo
  //pnm_mgr.WritePNMFile(imagem, "C:/Users/rodrigo/Desktop/", "CG_do_sucesso");
  //## 
  return 0;
}
