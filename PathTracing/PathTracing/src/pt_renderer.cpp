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
      Vector3d light_origin = this->scene_.point_lights_[i].position;
      Vector3d light_direction = intersection_point - light_origin;
      light_direction = light_direction / light_direction.norm();
      util::Ray from_light_src(light_origin, light_direction, 1);  // Depth here is useless.
      
      // Here we check if the ray 'from_light_src' is hitting the same point as the 'ray' (throwed by camera), 
      // and if the light is of same side that camera eye
      Vector3d curr_normal;
      util::RenderableObject *nearest_obj;
      double nearest;
      
      this->GetNearestObjectAndIntersection(from_light_src, &nearest_obj, &nearest, &curr_normal);

      Vector3d to_light = -light_direction;
      double cos_theta = normal.dot(to_light);

      if (intersection_point.isApprox(light_origin + nearest*light_direction, this->kEps) &&
          cos_theta > 0.0) {
        Vector3d reflected = 2*normal*cos_theta - to_light;
        reflected = reflected / reflected.norm();
        
        double cos_alpha = reflected.dot(viewer);
        double light_intensity = this->scene_.point_lights_[i].intensity;
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

    // TODO: Fazer o teste de qual raio passara: difuso, especular, ou transmitido
    /*Seja ktot = kd+ks+kt  e um número aleatório R entre (0,ktot). Se R < kd então dispara raio difuso, 
    senão se R < kd+ks dispara raio especular, senão dispara raio transmitido. Para a direção do raio especular utilizar R=2N(NL) - L.
    Para estabelecermos uma direção aleatória (phi, theta) para o raio difuso, precisamos de 2 números aleatórios R1 e R2 no intervalo (0,1).
    Então teremos phi=cos-1 (sqrt(R1)) e theta = 2.pi.R2.
    */
    double k_tot = obj_material.k_d + obj_material.k_s + obj_material.k_t;
    double ray_type = this->distribution_(this->generator_)*k_tot;
    Vector3d indirect_light;

    if (ray_type < obj_material.k_d) {        // Throw a diffuse ray
      // Generate a ray with random direction with origin on intesected point.
      // TODO(Mateus): utilizar a direção aleatória descrita na especificação.
      /*util::Ray new_diffuse_ray(intersection_point, std::cos(this->distribution_(this->generator_))*normal , ray.depth + 1);
      indirect_light = obj_material.k_d * this->TracePath(new_diffuse_ray);*/
      double r_1 = this->distribution_(this->generator_);
      double r_2 = this->distribution_(this->generator_);
      double phi = std::acos(std::sqrt(r_1));
      double theta = 2*M_PI*r_2;
    } else if (ray_type < obj_material.k_d + obj_material.k_s) { // Throw a specular ray
      Vector3d reflected = 2*normal*normal.dot(viewer) - viewer;
      reflected = reflected / reflected.norm();
      util::Ray new_specular_ray(intersection_point, reflected, ray.depth + 1);
      indirect_light = obj_material.k_s * this->TracePath(new_specular_ray);

    } else {                                  // Throw a refracted ray
      util::Ray new_refracted_ray(intersection_point, -viewer, ray.depth + 1);  // TODO: colocar indice de refracao (lei de Snell)
      indirect_light = obj_material.k_t * this->TracePath(new_refracted_ray);

    }

    color += indirect_light;  // TODO: A recursao nao devia para apenas depois de calcular a parte direta nao? pq como ta sempre vai ter um background somando

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

cv::Mat PTRenderer::RenderScene() {
  //cv::Mat rendered_image = cv::Mat::zeros(this->scene_.camera_.height_, this->scene_.camera_.width_, CV_8UC3);
  cv::Mat rendered_image = cv::Mat::zeros(this->scene_.camera_.height_, this->scene_.camera_.width_, CV_64FC3);

  double pixel_w = (this->scene_.camera_.top_(0) - this->scene_.camera_.bottom_(0)) / this->scene_.camera_.width_;
  double pixel_h = (this->scene_.camera_.top_(1) - this->scene_.camera_.bottom_(1)) / this->scene_.camera_.height_;

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
        
        //rendered_image.data[j*(rendered_image.cols*3)+k*3+0] = rendered_image.data[j*(rendered_image.cols*3)+k*3+0] + additional_color(2);
        //rendered_image.data[j*(rendered_image.cols*3)+k*3+1] = rendered_image.data[j*(rendered_image.cols*3)+k*3+1] + additional_color(1);
        //rendered_image.data[j*(rendered_image.cols*3)+k*3+2] = rendered_image.data[j*(rendered_image.cols*3)+k*3+2] + additional_color(0);
      }
    }
  }
  cv::imshow("virgin image", rendered_image);
  
  rendered_image = rendered_image / this->scene_.nmbr_paths_;
  cv::imshow("divided by N_paths image", rendered_image);

  this->ApplyToneMapping(rendered_image);
  cv::imshow("tone_mapped image", rendered_image);

  // Criando a imagem de unsigned ints
  //cv::Mat round_rendered_image = cv::Mat::zeros(this->scene_.camera_.height_, this->scene_.camera_.width_, CV_8UC3);
  //for (int j = 0; j < round_rendered_image.rows; ++j) {
  //    for (int k = 0; k < round_rendered_image.cols; ++k) {
  //      round_rendered_image.at<cv::Vec3d>(j, k)[0] = static_cast<unsigned int>( rendered_image.at<cv::Vec3d>(j, k)[0]);
  //      round_rendered_image.at<cv::Vec3d>(j, k)[1] = static_cast<unsigned int>( rendered_image.at<cv::Vec3d>(j, k)[1]);
  //      round_rendered_image.at<cv::Vec3d>(j, k)[2] = static_cast<unsigned int>( rendered_image.at<cv::Vec3d>(j, k)[2]);
  //    }
  //}
  ////rendered_image.convertTo(round_rendered_image, CV_8UC3);
  //cv::imshow("round image", round_rendered_image);


  std::cout << "Aperte alguma tecla para fechar as imagens" << std::endl;
  cv::waitKey(0);
  return rendered_image*255;
}

}  // namespace pt