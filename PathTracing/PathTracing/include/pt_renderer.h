// Copyright (c) 2015, Jubileus
//
// Project: Sucesso do verao
// Author: Rodrigo F. Figueiredo <rodrigo.figueiredo@gprt.ufpe.br>
// Creation date: 05/01/2016 (dd/mm/yyyy)

#ifndef PT_RENDERER_H
#define PT_RENDERER_H

#include <opencv2\core\core.hpp>

#include "sdl_object.h"

namespace pt {

class PTRenderer {
  public:
    PTRenderer() {};
    ~PTRenderer () {};

    cv::Mat RenderScene(const util::SDLObject &scene);

};

}  // namespace pt

#endif  // PT_RENDERER_H