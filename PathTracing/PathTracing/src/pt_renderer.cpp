#include "pt_renderer.h"

using Eigen::Vector3d;

namespace pt {

const double PTRenderer::kEps = 1.0e-03;

Vector3d PTRenderer::TracePath(const util::Ray &ray) {
  if (ray.depth > this->scene_.max_depth_) {
    return this->scene_.background_color_;
  } else {
    util::RenderableObject *object = nullptr;
    double t;
    Vector3d normal;

    this->GetNearestObjectAndIntersection(ray, &object, &t, &normal);
    Vector3d intersection_point = ray.origin + t*ray.direction;

    if (math::IsAlmostEqual(t, -1.0, this->kEps) || object == nullptr) {
      return this->scene_.background_color_;
    }
    
    util::Material obj_material = object->material();
    Vector3d material_color(obj_material.red, obj_material.green, obj_material.blue);
    Vector3d color = this->scene_.ambient_light_intensity_*obj_material.k_a*material_color;
    Vector3d viewer = intersection_point - this->scene_.camera_.eye_;
    viewer = viewer / viewer.norm();
    
    // Calcula-se a intensidade do objeto naquele ponto influenciada pelas fontes de luz na cena.
    for (int i = 0; i < this->scene_.point_lights_.size(); ++i) {
      Vector3d origin = this->scene_.point_lights_[i].position;
      util::Ray from_light_src(origin, intersection_point - origin, 1);  // Depth here is useless.
      
      Vector3d curr_normal;
      double nearest = object->GetIntersectionParameter(ray, curr_normal);

      if (curr_normal.dot(-from_light_src.direction) > 0.0 &&
          intersection_point.isApprox(origin + nearest*from_light_src.direction, this->kEps)) {
        Vector3d unit_direction = from_light_src.direction / from_light_src.direction.norm();
        double cos_alpha = unit_direction.dot(curr_normal);
        double sin_alpha = std::sqrt(1 - cos_alpha*cos_alpha);

        Vector3d reflected =  unit_direction*sin_alpha - unit_direction*cos_alpha;

        reflected = reflected / reflected.norm();
        double cos_theta = reflected.dot(viewer);

        color += obj_material.k_d*material_color +
          obj_material.k_s*std::pow(cos_theta, obj_material.n)*material_color;
      }
    }

    // TODO(Mateus): implementar depois!
    for (int i = 0; i < this->scene_.extense_lights_.size(); ++i) {
      
    }

    util::Ray new_ray(intersection_point, std::cos(this->distribution_(this->generator_))*normal , ray.depth + 1);
    color += this->TracePath(new_ray);  //TODO(Mateus): verificar BRDF.

    return color;
  }
}

void PTRenderer::ApplyToneMapping(cv::Mat &image) {
  for (int i = 0; i < image.rows; ++i) {
    for (int j = 0; j < image.cols; ++j) {
      image.data[i*(3*image.cols)+3*j]   /= image.data[i*(3*image.cols)+3*j]    + this->scene_.tone_mapping_;
      image.data[i*(3*image.cols)+3*j+1] /= image.data[i*(3*image.cols)+3*j+1]  + this->scene_.tone_mapping_;
      image.data[i*(3*image.cols)+3*j+2] /= image.data[i*(3*image.cols)+ 3*j+2] + this->scene_.tone_mapping_;
    }
  }

  image = 255*image;
}

void PTRenderer::GetNearestObjectAndIntersection(const util::Ray &ray,
                                                 util::RenderableObject **object,
                                                 double *parameter,
                                                 Eigen::Vector3d *normal) {
  *parameter = std::numeric_limits<double>::max();

  for (int i = 0; i < this->scene_.quadrics_objects_.size(); ++i) {
    Vector3d curr_normal;
    double curr_t = this->scene_.quadrics_objects_[i].GetIntersectionParameter(ray, curr_normal);
    
    if (*parameter > curr_t) {
      *object = &this->scene_.quadrics_objects_[i];
      *parameter = curr_t;
      *normal = curr_normal;
    }
  }

  for (int i = 0; i < this->scene_.triangular_objects_.size(); ++i) {
    Vector3d curr_normal;
    double curr_t = this->scene_.triangular_objects_[i].GetIntersectionParameter(ray, curr_normal);
    
    if (*parameter > curr_t) {
      *object = &this->scene_.triangular_objects_[i];
      *parameter = curr_t;
      *normal = curr_normal;
    }
  }
}

cv::Mat PTRenderer::RenderScene() {
  cv::Mat rendered_image = cv::Mat::zeros(this->scene_.camera_.height_, this->scene_.camera_.width_, CV_8UC3);

  double pixel_w = (this->scene_.camera_.top_(0) - this->scene_.camera_.bottom_(0)) / this->scene_.camera_.width_;
  double pixel_h = (this->scene_.camera_.top_(1) - this->scene_.camera_.bottom_(1)) / this->scene_.camera_.height_;

  for (int i = 0; i < this->scene_.nmbr_paths_; ++i) {
    // Dispara um raio n vezes em um determinado pixel.
    for (int j = 0; j < rendered_image.rows; ++j) {
      for (int k = 0; k < rendered_image.cols; ++k) {
        Vector3d looking_at((this->scene_.camera_.bottom_(0) + pixel_w / 2) + k*pixel_w,
                            (this->scene_.camera_.bottom_(1) + pixel_w / 2) + j*pixel_h,
                            0.0);
        Vector3d direction = looking_at - this->scene_.camera_.eye_;
        direction = direction / direction.norm();

        util::Ray ray(this->scene_.camera_.eye_, direction, 1);
        Vector3d additional_color = this->TracePath(ray);

        rendered_image.data[j*(rendered_image.cols*3)+k*3+0] += additional_color(2);
        rendered_image.data[j*(rendered_image.cols*3)+k*3+1] += additional_color(1);
        rendered_image.data[j*(rendered_image.cols*3)+k*3+2] += additional_color(0);
      }
    }
  }

  rendered_image = rendered_image / this->scene_.nmbr_paths_;
  this->ApplyToneMapping(rendered_image);

  return rendered_image;
}

}  // namespace pt