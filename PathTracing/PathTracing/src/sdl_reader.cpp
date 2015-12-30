#include "sdl_reader.h"

namespace io {
  void SDLReader::ReadSDL(std::string file_directory, std::string file_name,
                          util::SDLObject &sdl_object) {
    // Preparar todas as variaveis que receberao os dados para a criacao do sdl_object
    std::string                         output_name;
    Eigen::Vector3d                     eye;
    Eigen::Vector2d                     bottom;
    Eigen::Vector2d                     top;
    double                              width;
    double                              height;
    Eigen::Vector3d                     background_color;
    double                              ambient_light_intensity;
    std::vector<util::PointLight>       point_lights;
    std::vector<util::TriangularObject> extense_lights;
    int                                 nmbr_paths;
    int                                 max_depth;
    double                              tone_mapping;
    int                                 random_seed;
    std::vector<util::Quadric>          quadrics_objects;
    std::vector<util::TriangularObject> triangular_objects;

    std::ifstream sdl_file(file_directory + (file_name + ".sdl"));
    std::string word;

    // Ler todas as linhas
    if (sdl_file.is_open()) {
      while (sdl_file >> word) {
        if (word[0] == '#') {                  // Comentario
          std::string commentary;
          std::getline(sdl_file, commentary);
          std::cout << "lido #" << commentary << std::endl;

        } else if (word == "output") {         // Nome do arquivo de saida
          sdl_file >> output_name;

        } else if (word == "eye") {            // Centro da camera
          double x, y, z;
          sdl_file >> x;
          sdl_file >> y;
          sdl_file >> z;
          eye[0] = x;
          eye[1] = y;
          eye[2] = z;
          std::cout << "Lido 'eye': " << x << " " << y << " " << z << std::endl;

        } else if (word == "ortho") {          // Viewport da camera
          double bx, by, tx, ty;
          sdl_file >> bx;
          sdl_file >> by;
          sdl_file >> tx;
          sdl_file >> ty;
          bottom[0] = bx;
          bottom[1] = by;
          top[0] = tx;
          top[1] = ty;
          std::cout << "Lido 'ortho': " << bx << " " << by << " " << tx << " " << ty << std::endl;

        } else if (word == "size") {           // Quantidade de pixels na horizontal e vertical
          sdl_file >> width;
          sdl_file >> height;
          std::cout << "Lido 'size': " << width << " " << height << std::endl;

        } else if (word == "background") {     // Cor do fundo da imagem
          double r, g, b;
          sdl_file >> r;
          sdl_file >> g;
          sdl_file >> b;
          background_color[0] = r;
          background_color[1] = g;
          background_color[2] = b;
          std::cout << "Lido 'background': " << r << " " << g << " " << b << std::endl;

        } else if (word == "ambient") {        // Intensidade da componente ambiente
          sdl_file >> ambient_light_intensity;
          std::cout << "Lido 'ambient': " << ambient_light_intensity << std::endl;

        } else if (word == "light") {          // Luzes extensas
          // TODO(figueiredo) implementar direito
          std::string commentary;
          std::getline(sdl_file, commentary);
          std::cout << "@" << commentary << std::endl;

        } else if (word == "pointlight") {     // Luzes pontuais
                                               // Nao existe no formato original da especificacao
          // TODO(figueiredo) implementar direito
          std::string commentary;
          std::getline(sdl_file, commentary);
          std::cout << "@" << commentary << std::endl;

        } else if (word == "npaths") {         // Numero de raios por pixel
          sdl_file >> nmbr_paths;
          std::cout << "Lido 'npaths': " << nmbr_paths << std::endl;

        } else if (word == "maxdepth") {       // Quantidade maxima de reflexoes
                                               // Nao existe no formato original da especificacao
          sdl_file >> max_depth;
          std::cout << "Lido 'maxdepth': " << max_depth << std::endl;

        } else if (word == "tonemapping") {    // Regulador de iluminacao -para pos processamento
          sdl_file >> tone_mapping;
          std::cout << "Lido 'tonemapping': " << tone_mapping << std::endl;

        } else if (word == "seed") {           // Semente inteira do gerador randomico
          sdl_file >> random_seed;
          std::cout << "Lido 'seed': " << random_seed << std::endl;

        } else if (word == "objectquadric") {  // Objetos baseados em surpecifies parametricas
          double a, b, c, d, e, f, g, h, j, k, red, green, blue, ka, kd, ks, kt, n;
          sdl_file >> a;
          sdl_file >> b;
          sdl_file >> c;
          sdl_file >> d;
          sdl_file >> e;
          sdl_file >> f;
          sdl_file >> g;
          sdl_file >> h;
          sdl_file >> j;
          sdl_file >> k;
          sdl_file >> red;
          sdl_file >> green;
          sdl_file >> blue;
          sdl_file >> ka;
          sdl_file >> kd;
          sdl_file >> ks;
          sdl_file >> kt;
          sdl_file >> n;

          util::Material new_material(red, green, blue, ka, kd, ks, kt, n);
          util::Quadric  new_quadric(a, b, c, d, e, f, g, h, j, k, new_material);
          std::cout << "Lido 'objectquadric': " << "te dana, mt coisa pra imprimir :P" <<std::endl;
          sdl_object.quadrics_objects_.push_back(new_quadric);
        } else if (word == "object") {         // Objetos baseados em malhas trianguladas
          // TODO(figueiredo) implementar direito
          std::string commentary;
          std::getline(sdl_file, commentary);
          std::cout << "@" << commentary << std::endl;

        } else {
          std::cout << "  BORA BOY! tolken nao suportado: " << word << std::endl;
          std::cout << "    Leitura interrompida." << std::endl;
          sdl_file.close();
          return;  // BORA BOY! tolken nao suportado
        }
      }

      std::cout << "## Arquivo SDL lido com sucesso" << std::endl;
      sdl_file.close();
    }
  }
}  // namespace io