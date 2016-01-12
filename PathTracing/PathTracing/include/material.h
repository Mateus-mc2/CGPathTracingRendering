#ifndef MATERIAL_H_
#define MATERIAL_H_

#include "util_exception.h"

namespace util {

class InvalidMaterialCoefficientsException : public UtilException {
  public:
    InvalidMaterialCoefficientsException(const std::string &error_msg)
        : UtilException(error_msg) {}
};

struct Material {
  Material () {}
  Material(const double &r, const double &g, const double &b, const double &refraction_coeff,
           const double &k_a, const double &k_d, const double &k_s, const double &k_t,
           const int &n, const double &lp, const double &light_sampling_step,
           const double &light_density);
  ~Material() {}
  
  double red;   // Componente vermelha.
  double green; // Componente verde.
  double blue;  // Componente azul.
  double refraction_coeff;  // Coeficiente de refração.
  
  // Relevantes aos objetos iluminados
  double k_a;  // Coeficiente de reflexão ambiente.
  double k_d;  // Coeficiente de reflexão difusa.
  double k_s;  // Coeficiente de reflexão especular.
  double k_t;  // Coeficiente de transparência.
  int n;       // Expoente de reflexão especular.

  // Relevantes aos objetos luminosos
  double lp;
  double light_sampling_step;
  double light_density;
};

}  // namespace util

#endif  // MATERIAL_H_