#include "pt_renderer.h"

using Eigen::Vector3d;

namespace pt {

Vector3d PTRenderer::TracePath(const util::Ray &ray) {
  if (ray.depth >= this->kScene.max_depth_) {
    return Vector3d(0, 0, 0);
  } else {
    util::RenderableObject *object;
    double parameter;

    this->GetNearestObjectAndIntersection(ray, &object, &parameter);
    // return Vector3d(0, 0, 0);
  }
}

void PTRenderer::ApplyToneMapping(cv::Mat &image) {
  image = image / this->kScene.nmbr_paths_;
}

void PTRenderer::GetNearestObjectAndIntersection(const util::Ray &ray,
                                                 util::RenderableObject **object,
                                                 double *parameter) {
  for (int i = 0; i < this->kScene.quadrics_objects_.size(); ++i) {
    *object = &this->kScene.quadrics_objects_[i];
    //*object = new util::Quadric(this->kScene.quadrics_objects_[i]);
  }

  for (int i = 0; i < this->kScene.triangular_objects_.size(); ++i) {

  }
}

cv::Mat PTRenderer::RenderScene() {
  cv::Mat rendered_image = cv::Mat::zeros(this->kScene.camera_.height_, this->kScene.camera_.width_, CV_8UC3);

  double pixel_w = (this->kScene.camera_.top_(0) - this->kScene.camera_.bottom_(0)) / this->kScene.camera_.width_;
  double pixel_h = (this->kScene.camera_.top_(1) - this->kScene.camera_.bottom_(1)) / this->kScene.camera_.height_;

  for (int i = 0; i < this->kScene.nmbr_paths_; ++i) {
    // Dispara um raio n vezes em um determinado pixel.
    for (int j = 0; j < rendered_image.rows; ++j) {
      for (int k = 0; k < rendered_image.cols; ++k) {
        Vector3d looking_at((this->kScene.camera_.bottom_(0) + pixel_w / 2) + k*pixel_w,
                            (this->kScene.camera_.bottom_(1) + pixel_w / 2) + j*pixel_h,
                            0.0);
        Vector3d direction = looking_at - this->kScene.camera_.eye_;
        direction = direction / direction.norm();

        util::Ray ray(this->kScene.camera_.eye_, direction, 1);
        Vector3d additional_color = this->TracePath(ray);

        rendered_image.data[j*(rendered_image.cols*3)+k*3+0] += static_cast<uchar>(255*additional_color(0));
        rendered_image.data[j*(rendered_image.cols*3)+k*3+1] += static_cast<uchar>(255*additional_color(1));
        rendered_image.data[j*(rendered_image.cols*3)+k*3+2] += static_cast<uchar>(255*additional_color(2));
      }
    }
  }

  this->ApplyToneMapping(rendered_image);

  return rendered_image;
}

}  // namespace pt