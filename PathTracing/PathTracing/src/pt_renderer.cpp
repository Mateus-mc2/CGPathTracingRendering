#include <Eigen\Dense>

#include "pt_renderer.h"

namespace pt {

cv::Mat PTRenderer::RenderScene (const util::SDLObject &scene) {
  cv::Mat rendered_image(scene.camera_.height_, scene.camera_.height_, CV_8UC3);


  return rendered_image;
}

}  // namespace pt