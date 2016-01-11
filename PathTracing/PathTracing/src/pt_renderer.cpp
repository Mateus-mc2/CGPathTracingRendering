#include "pt_renderer.h"

using Eigen::Vector3d;

namespace pt {

const double PTRenderer::kEps = 1.0e-03;

Vector3d PTRenderer::TracePath(const util::Ray &ray) {
  if (ray.depth >= this->scene_.max_depth_) {
    return Vector3d(0.0, 0.0, 0.0);
  } else {
    util::RenderableObject *object = nullptr;
    double t;
    Vector3d normal;

    this->GetNearestObjectAndIntersection(ray, &object, &t, &normal);
    Vector3d intersection_point = ray.origin + t*ray.direction;

    if (math::IsAlmostEqual(t, -1.0, this->kEps) || object == nullptr) {
      return this->scene_.background_color_;
    }

    //std::cout << "Bla" << std::endl;

    util::Material obj_material = object->material();
    Vector3d material_color(obj_material.red, obj_material.green, obj_material.blue);
    Vector3d color = this->scene_.ambient_light_intensity_*obj_material.k_a*material_color;

    Vector3d viewer = this->scene_.camera_.eye_ - intersection_point;
    viewer = viewer / viewer.norm();
    
    // Calcula-se a intensidade do objeto naquele ponto influenciada pelas fontes de luz na cena.
    for (int i = 0; i < this->scene_.point_lights_.size(); ++i) {
      Vector3d light_direction =  this->scene_.point_lights_[i].position - intersection_point;
      light_direction = light_direction / light_direction.norm();
      util::Ray shadow_ray(intersection_point, light_direction, ray.ambient_objs, 1);
      
      // Here we compute how much light is blockd by opaque and transparent surfaces, and use the
      // resulting intensity to scale diffuse and specular terms of the final color.
      double light_intensity = this->ScaleLightIntensity(this->scene_.point_lights_[i],
                                                         shadow_ray);
      double cos_theta = normal.dot(light_direction);

      if (cos_theta > 0.0) {
        Vector3d reflected = 2*normal*cos_theta - light_direction;
        reflected = reflected / reflected.norm();
        
        double cos_alpha = reflected.dot(viewer);
        Vector3d light_color(this->scene_.point_lights_[i].red, 
                             this->scene_.point_lights_[i].green,
                             this->scene_.point_lights_[i].blue);
       
        color += light_intensity*(obj_material.k_d*material_color*cos_theta +
                 obj_material.k_s*light_color*std::pow(cos_theta, obj_material.n));
      }
    }

    // TODO(Mateus): implementar depois!
    for (int i = 0; i < this->scene_.extense_lights_.size(); ++i) {
      
    }

    /*Seja ktot = kd+ks+kt  e um número aleatório R entre (0,ktot). Se R < kd então dispara raio difuso, 
    senão se R < kd+ks dispara raio especular, senão dispara raio transmitido. Para a direção do raio especular utilizar R=2N(NL) - L.
    Para estabelecermos uma direção aleatória (phi, theta) para o raio difuso, precisamos de 2 números aleatórios R1 e R2 no intervalo (0,1).
    Então teremos phi=cos-1 (sqrt(R1)) e theta = 2.pi.R2.
    */
    double k_tot = obj_material.k_d + obj_material.k_s + obj_material.k_t;
    double ray_type = this->distribution_(this->generator_)*k_tot;
    Vector3d indirect_light;

    if (ray_type < obj_material.k_d) {  // Throw a diffuse ray
      // Generate a ray with random direction with origin on intesected point (using uniform sphere distribution here).
      double r_1 = this->distribution_(this->generator_);
      double r_2 = this->distribution_(this->generator_);
      double phi = std::acos(std::sqrt(r_1));
      double theta = 2*M_PI*r_2;
      
      Vector3d rand_direction(std::sin(phi)*std::cos(theta),
                              std::sin(phi)*std::sin(theta),
                              std::cos(phi));
      util::Ray new_diffuse_ray(intersection_point, rand_direction, ray.ambient_objs, ray.depth + 1);
      indirect_light = obj_material.k_d * this->TracePath(new_diffuse_ray);
    } else if (ray_type < obj_material.k_d + obj_material.k_s) {  // Throw a specular ray
      Vector3d reflected = 2*normal*normal.dot(viewer) - viewer;
      reflected = reflected / reflected.norm();
      util::Ray new_specular_ray(intersection_point, reflected, ray.ambient_objs, ray.depth + 1);
      indirect_light = obj_material.k_s * this->TracePath(new_specular_ray);
    } else {  // Throw a refracted ray
      double n_1;
      double n_2;
      std::stack<util::RenderableObject*> objs_stack(ray.ambient_objs);

      if (objs_stack.empty()) {  // Scene's ambient refraction coefficient (we're assuming n = 1.0 here).
        n_1 = 1.0;
        n_2 = obj_material.refraction_coeff;
        objs_stack.push(object);
      } else {  // Ray is getting out of current object.
        util::RenderableObject* last_obj = objs_stack.top();
        n_1 = last_obj->material().refraction_coeff;

        if (object != last_obj) {
          n_2 = obj_material.refraction_coeff;
          objs_stack.push(object);
        } else {
          objs_stack.pop();
          n_2 = objs_stack.top()->material().refraction_coeff;
        }        
      }

      double cos_theta_incident = normal.dot(-ray.direction);
      
      if (cos_theta_incident < 0.0) {
        normal = -normal;
        cos_theta_incident = -cos_theta_incident;
      }

      double sin_theta_incident = std::sqrt(1 - cos_theta_incident*cos_theta_incident);

      // Check if it's a total internal reflection.
      if (sin_theta_incident < (n_2 / n_1)) {
        // Get new refracted ray.
        double n_r = n_1 / n_2;
        Vector3d refracted = (n_r*cos_theta_incident - std::sqrt(1 - 
                              std::pow(n_r, 2)*std::pow(sin_theta_incident, 2)))*normal
                              + n_r*ray.direction;
        // Need to update the stack of objects.
        util::Ray new_refracted_ray(intersection_point, refracted, objs_stack, ray.depth + 1);
        indirect_light = obj_material.k_t * this->TracePath(new_refracted_ray);
      } else {  // Simulate total internal reflection.
        Vector3d total_internal_reflected = 2*normal*cos_theta_incident + ray.direction;
        util::Ray total_internal_reflection_ray(intersection_point,
                                                total_internal_reflected,
                                                ray.ambient_objs,
                                                ray.depth + 1);
        indirect_light = obj_material.k_s * this->TracePath(total_internal_reflection_ray);
      }
    }

    color += indirect_light;  // TODO: A recursao nao devia para apenas depois de calcular a parte direta nao? pq como ta sempre vai ter um background somando
    
    // Clamp to guarantee that every ray returns a color r,g,b <= 1
    color[0] = std::min(color[0], 1.0);
    color[1] = std::min(color[1], 1.0);
    color[2] = std::min(color[2], 1.0);

    return color;
  }
}

void PTRenderer::ApplyToneMapping(cv::Mat &image) {
  for (int i = 0; i < image.rows; ++i) {
    for (int j = 0; j < image.cols; ++j) {
      image.at<cv::Vec3d>(i, j)[0] /= image.at<cv::Vec3d>(i, j)[0] + this->scene_.tone_mapping_;
      image.at<cv::Vec3d>(i, j)[1] /= image.at<cv::Vec3d>(i, j)[1] + this->scene_.tone_mapping_;
      image.at<cv::Vec3d>(i, j)[2] /= image.at<cv::Vec3d>(i, j)[2] + this->scene_.tone_mapping_;
    }
  }
}

void PTRenderer::GetNearestObjectAndIntersection(const util::Ray &ray,
                                                 util::RenderableObject **object,
                                                 double *parameter,
                                                 Eigen::Vector3d *normal) {
  *parameter = std::numeric_limits<double>::max();
  Vector3d curr_normal;

  for (int i = 0; i < this->scene_.quadrics_objects_.size(); ++i) {
    double curr_t = this->scene_.quadrics_objects_[i].GetIntersectionParameter(ray, curr_normal);
    
    if (*parameter > curr_t && curr_t > 0.0) {
      *object = &this->scene_.quadrics_objects_[i];
      *parameter = curr_t;
      *normal = curr_normal;
    }
  }

  for (int i = 0; i < this->scene_.triangular_objects_.size(); ++i) {
    double curr_t = this->scene_.triangular_objects_[i].GetIntersectionParameter(ray, curr_normal);
    
    if (*parameter > curr_t && curr_t > 0.0) {
      *object = &this->scene_.triangular_objects_[i];
      *parameter = curr_t;
      *normal = curr_normal;
    }
  }
}


double PTRenderer::ScaleLightIntensity(const util::PointLight &curr_light,
                                       const util::Ray &shadow_ray) {
  double final_intensity = curr_light.intensity;
  Vector3d normal;
  // Pega o índice do vetor diretor tal que v(idx) não seja zero. Ficou feioso assim, mas dane-se...
  const int idx = shadow_ray.direction(0) != 0? 0 : (shadow_ray.direction(1) != 0 ? 1 : 2);
  assert(shadow_ray.direction(idx) != 0);
  const double kMaxT = (curr_light.position(idx) - shadow_ray.origin(idx)) / shadow_ray.direction(idx);

  for (int i = 0; i < this->scene_.quadrics_objects_.size(); ++i) {
    util::Material obj_material = this->scene_.quadrics_objects_[i].material();
    double t = this->scene_.quadrics_objects_[i].GetIntersectionParameter(shadow_ray, normal);

    if (t > 0.0 && t < kMaxT) {
      final_intensity *= obj_material.k_t;
    }
  }

  for (int i = 0; i < this->scene_.triangular_objects_.size(); ++i) {
    util::Material obj_material = this->scene_.triangular_objects_[i].material();
    double t = this->scene_.triangular_objects_[i].GetIntersectionParameter(shadow_ray, normal);

    if (t > 0.0 && t < kMaxT) {
      final_intensity *= obj_material.k_t;
    }
  }

  return final_intensity;
}

cv::Mat PTRenderer::RenderScene() {
  //cv::Mat rendered_image = cv::Mat::zeros(this->scene_.camera_.height_, this->scene_.camera_.width_, CV_8UC3);
  cv::Mat rendered_image = cv::Mat::zeros(this->scene_.camera_.height_, this->scene_.camera_.width_, CV_64FC3);
  cv::Mat parcial_result = cv::Mat::zeros(this->scene_.camera_.height_, this->scene_.camera_.width_, CV_64FC3);
  double pixel_w = (this->scene_.camera_.top_(0) - this->scene_.camera_.bottom_(0)) / this->scene_.camera_.width_;
  double pixel_h = (this->scene_.camera_.top_(1) - this->scene_.camera_.bottom_(1)) / this->scene_.camera_.height_;
  
  int percent;

  for (int i = 0; i < this->scene_.nmbr_paths_; ++i) {
    // Dispara um raio n vezes em um determinado pixel.
    for (int j = 0; j < rendered_image.rows; ++j) {
      for (int k = 0; k < rendered_image.cols; ++k) {
        Vector3d looking_at((this->scene_.camera_.bottom_(0) + pixel_w / 2) + k*pixel_w,
                            (this->scene_.camera_.top_(1)    - pixel_h / 2) - j*pixel_h,
                            0.0);
        Vector3d direction = looking_at - this->scene_.camera_.eye_;
        direction = direction / direction.norm();

        util::Ray ray(this->scene_.camera_.eye_, direction, 1);
        Vector3d additional_color = this->TracePath(ray);

        rendered_image.at<cv::Vec3d>(j, k)[0] += additional_color(2);
        rendered_image.at<cv::Vec3d>(j, k)[1] += additional_color(1);
        rendered_image.at<cv::Vec3d>(j, k)[2] += additional_color(0);        
      }
    }

    percent = int(100.0*(double(i) / this->scene_.nmbr_paths_));
    std::cout << "\r" << percent << "% completed: ";
    std::cout << std::string(percent/10, '@') << std::string(10 - percent/10, '=');
    std::cout.flush();
    parcial_result = rendered_image / i;
    cv::imshow("parcial_result", parcial_result);
    cv::waitKey(1);
  }

  std::cout << std::endl;
  cv::imshow("virgin image", rendered_image);
  
  rendered_image = rendered_image / this->scene_.nmbr_paths_;
  cv::imshow("divided by N_paths image", rendered_image);

  this->ApplyToneMapping(rendered_image);
  cv::imshow("tone_mapped image", rendered_image);

  std::cout << "Aperte alguma tecla para fechar as imagens" << std::endl;
  cv::waitKey(0);
  return rendered_image*255;
}

}  // namespace pt