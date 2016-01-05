// Copyright (c) 2015, Jubileus
//
// Project: Sucesso do verao
// Author: Rodrigo F. Figueiredo <rodrigo.figueiredo@gprt.ufpe.br>
// Creation date: 05/01/2016 (dd/mm/yyyy)

#include <iostream>

#include <opencv2\core\core.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2\highgui\highgui.hpp>

#include "pt_renderer.h"
#include "pnm_writer.h"
#include "sdl_reader.h"

int main(int argc, char* argv[]) {
  // Leitura da cena
  std::cout << "\n## Come�o da leitura do arquivo." << std::endl;
  io::SDLReader sdl_reader;
  util::SDLObject sdl_object = sdl_reader.ReadSDL(argv[1], argv[2]);

  // Processamento
  std::cout << "\n## Come�o da renderiza��o." << std::endl;
  pt::PTRenderer pt_renderer;
  cv::Mat rendered_img = pt_renderer.RenderScene(sdl_object);

  // Escrita da imagem renderizada
  std::cout << "\n## Come�o da exporta��o." << std::endl;
  io::PNMWriter pnm_mgr(argv[3]);
  if (argc > 4)
    pnm_mgr.WritePNMFile(rendered_img, argv[3], argv[4]);
  else
    pnm_mgr.WritePNMFile(rendered_img);

  return 0;
}
