// Copyright (c) 2015, Jubileus
//
// Project: Sucesso do verao
// Author: Rodrigo F. Figueiredo <rodrigo.figueiredo@gprt.ufpe.br>
// Creation date: 29/12/2015 (dd/mm/yyyy)

#ifndef SDL_OBJECT_H_
#define SDL_OBJECT_H_

#include <Eigen/Dense>

#include "camera.h"
#include "light.h"
#include "quadric.h"
#include "triangular_object.h"

#include <vector>
#include <string>

namespace util {

struct SDLObject {
  SDLObject(const std::string &file_name,
            const Camera &camera,
            const Eigen::Vector3d &background_color,
            const double &ambient_light_intensity,
            const std::vector<PointLight> &point_lights,
            const std::vector<TriangularObject> extense_lights,
            const int &nmbr_paths,
            const int &max_depth,
            const double &tone_mapping,
            const int &random_seed,
            const std::vector<Quadric> &quadrics_objects,
            const std::vector<TriangularObject> &triangular_objects)
      : file_name_(file_name),
        camera_(camera),
        background_color_(background_color),
        ambient_light_intensity_(ambient_light_intensity),
        point_lights_(point_lights),
        extense_lights_(extense_lights),
        nmbr_paths_(nmbr_paths),
        max_depth_(max_depth),
        tone_mapping_(tone_mapping),
        random_seed_(random_seed),
        quadrics_objects_(quadrics_objects),
        triangular_objects_(triangular_objects) {}

  ~SDLObject();

  const std::string                   file_name_;
  const Camera                        camera_;
  const Eigen::Vector3d               background_color_;
  const double                        ambient_light_intensity_;
  const std::vector<PointLight>       point_lights_;
  const std::vector<TriangularObject> extense_lights_;
  const int                           nmbr_paths_;
  const int                           max_depth_;
  const double                        tone_mapping_;
  const int                           random_seed_;
  const std::vector<Quadric>          quadrics_objects_;
  const std::vector<TriangularObject> triangular_objects_;
};

}  // namespace util

#endif  // SDL_OBJECT_H_