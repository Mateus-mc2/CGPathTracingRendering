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
  Material (const double &r, const double &g, const double &b, const double &k_a,
            const double &k_d, const double &k_s, const int &n)
      :  r(r),
         g(g),
         b(b),
         k_a(k_a),
         k_d(k_d),
         k_s(k_s),
         n(n) {}
  ~Material() {}
  
  double r;  // Componente vermelha.
  double g;  // Componente verde.
  double b;  // Componente azul.
  
  double k_a;  // Coeficiente de reflex�o ambiente.
  double k_d;  // Coeficiente de reflex�o difusa.
  double k_s;  // Coeficiente de reflex�o especular.

  int n;  // Expoente de reflex�o especular.
};

}  // namespace util

#endif  // MATERIAL_H_