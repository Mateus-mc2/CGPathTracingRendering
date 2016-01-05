// Copyright (c) 2015, Jubileus
//
// Project: Sucesso do verao
// Author: Rodrigo F. Figueiredo <rodrigo.figueiredo@gprt.ufpe.br>
// Creation date: 05/01/2016 (dd/mm/yyyy)

#ifndef PT_RENDERER_H
#define PT_RENDERER_H

#include <Eigen\Dense>
#include <opencv2\core\core.hpp>

#include "ray.h"
#include "renderable_object.h"
#include "sdl_object.h"

namespace pt {

class PTRenderer {
  public:
    explicit PTRenderer(const util::SDLObject &scene) :kScene(scene) {};
    ~PTRenderer () {};

    cv::Mat RenderScene();
  private:
    util::SDLObject kScene;

    Eigen::Vector3d TracePath(const util::Ray &ray);
    void ApplyToneMapping(cv::Mat &image);
    void GetNearestObjectAndIntersection(const util::Ray &ray,
                                         util::RenderableObject **object,
                                         double *parameter);
};

}  // namespace pt

#endif  // PT_RENDERER_H