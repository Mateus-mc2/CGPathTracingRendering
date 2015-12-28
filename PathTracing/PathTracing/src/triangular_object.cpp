#include "triangular_object.h"

namespace util {

TriangularObject::TriangularObject(const Material &material,
                                   const std::vector<Eigen::Vector3d> &vertices,
                                   const std::vector<Eigen::Vector3i> &faces,
                                   const bool emissive_flag)
      :  kMaterial(material),
         kVertices(vertices),
         kFaces(faces),
         kEmissive(emissive_flag) {
  
}

}  // namespace util