#ifndef QUADRIC_H_
#define QUADRIC_H_

#include <Eigen\Dense>

#include <string>

#include "material.h"
#include "math_lib.h"
#include "ray.h"
#include "util_exception.h"

namespace util {

class InvalidCoefficientsVectorException : public UtilException {
  public:
    InvalidCoefficientsVectorException(const std::string &error_msg) : UtilException(error_msg) {}
};

class Quadric {
  public:
    // See http://mathworld.wolfram.com/QuadraticSurface.html to understand this notation.
    enum Index {kCoeffA, kCoeffB, kCoeffC,
                kCoeffF, kCoeffG, kCoeffH,
                kCoeffP, kCoeffQ, kCoeffR,
                kCoeffD};

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
    static const double kEps;

    Eigen::VectorXd coefficients_;
    Material material_;    
};

}  // namespace util

#endif  // QUADRIC_H_