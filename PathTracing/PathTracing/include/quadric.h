#ifndef QUADRIC_H_
#define QUADRIC_H_

#include <Eigen\Dense>

#include <string>

#include "material.h"
#include "ray.h"

namespace util {

enum Index {
  a,
  b,
  c,
  d,
  e,
  f,
  g,
  h,
  j,
  k
};

class InvalidCoefficientsVectorException : public std::exception {
  private:
    const std::string kErrorMsg;
  public:
    explicit InvalidCoefficientsVectorException(const std::string &error) : kErrorMsg(error) {}

    const char* what() const;
};

class Quadric {
  public:
    Quadric(const Eigen::VectorXd &coefficients, const Material &material);
    Quadric(const double &a, const double &b, const double &c, const double &d,
            const double &e, const double &f, const double &g, const double &h,
            const double &j, const double &k, const Material &material);
    ~Quadric() {}

    double GetIntesectionParameter(const Ray &ray);

    // Accessors.
    Eigen::VectorXd coefficients() const { return this->coefficients_; }
    Material material() const { return this->material_; }
  private:
    Eigen::VectorXd coefficients_;
    Material material_;
};

}  // namespace util

#endif  // QUADRIC_H_