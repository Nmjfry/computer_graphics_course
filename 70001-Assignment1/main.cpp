/*//////////////////////////////////////////////////////////////////////////
Author: Abhijeet Ghosh
Year: 2013
//////////////////////////////////////////////////////////////////////////*/
#include "loadPNM.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979323
#define uint unsigned int

#include <iostream>

using namespace std;

unsigned int width;
unsigned int height;
unsigned int numComponents;

void CreateAndSavePFM(const char *image_out) {
  width = 511; // set size of image to 511x511 pixels
  height = 511;
  numComponents = 3;

  float *img_out = new float[width * height * numComponents];

  for (uint i = 0; i < width; ++i) 
  {
    for (uint j = 0; j < height; ++j) 
    {

      for (uint k = 0; k < numComponents;
           ++k) // color channels - 3 for RGB images
      {
        uint index = i * width * numComponents + j * numComponents +
                     k; // index within the image

        // set image to black
        img_out[index] = 0.0f;

        float radius = 511.0f / 2.0f;
        float center = radius;

        float x = center - j;
        float y = center - i;
        float z = sqrt(radius * radius - x * x - y * y);
        float n_len = sqrt(x * x + y * y + z * z);

        float n[3] = {-x / n_len, y / n_len, z / n_len};

        float v[3] = {0.0f, 0.0f, 1.0f};
        float ndotv = n[0] * v[0] + n[1] * v[1] + n[2] * v[2];
        float r[3] = {2 * ndotv * n[0] - v[0], 2 * ndotv * n[1] - v[1],
                      2 * ndotv * n[2] - v[2]};
        float r_len = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

        if (sqrt((i - center) * (i - center) + (j - center) * (j - center)) <
            radius) {
          img_out[index] = (r[k] / r_len) * 0.5f + 0.5f;
        }
      }
    }
  }

  WritePFM(image_out, width, height, numComponents, img_out);
}

void gamma_correction(const char *image_out, const char *image_in,
                      float gamma) {

  float *img_in = loadPFM(image_in, width, height, numComponents);
  float *img_out = new float[width * height * numComponents];
  for (uint i = 0; i < height; ++i) // height
  {
    for (uint j = 0; j < width; ++j) // width
    {
      for (uint k = 0; k < numComponents;
           ++k) // color channels - 3 for RGB images
      {
        uint index = i * width * numComponents + j * numComponents +
                     k; // index within the image

        img_out[index] = pow(img_in[index], 1.0f / gamma);
      }
    }
  }
  WritePFM(image_out, width, height, numComponents, img_out);
}

void ReadProbeAndSavePFM(const char *probe, const char *image_out) {

  unsigned int p_width;
  unsigned int p_height;
  unsigned int p_numComponents;
  float *img_in = loadPFM(probe, p_width, p_height, p_numComponents);
  printf("width: %d, height: %d, numComponents: %d\n", p_width, p_height,
         p_numComponents);

  width = 511; // set size of image to 511x511 pixels
  height = 511;
  numComponents = 3;
  printf("width: %d, height: %d, numComponents: %d\n", width, height,
         numComponents);
  float *img_out = new float[width * height * numComponents];

  for (uint i = 0; i < width; ++i)
  {
    for (uint j = 0; j < height; ++j)
    {

      for (uint k = 0; k < numComponents;
           ++k) // color channels - 3 for RGB images
      {
        uint index = i * height * numComponents + j * numComponents +
                     k; // index within the image

        // set image to black
        img_out[index] = 0.0f;

        float radius = 511.0f / 2.0f;
        float center = radius;

        float y = center - j;
        float x = center - i;
        float z = sqrt(radius * radius - x * x - y * y);

        float n_len = sqrt(x * x + y * y + z * z);

        float n[3] = {-x / n_len, -y / n_len, z / n_len};
        float v[3] = {0.0f, 0.0f, 1.0f};

        float ndotv = n[0] * v[0] + n[1] * v[1] + n[2] * v[2];

        float r[3] = {2 * ndotv * n[0] - v[0], 2 * ndotv * n[1] - v[1],
                      2 * ndotv * n[2] - v[2]};

        float rad = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        float phi = atan2(r[1], r[0]);
        float theta = acos(r[2] / rad);

        // printf("phi: %f, theta: %f\n", phi, theta);


        uint lat = static_cast<uint>(((phi + M_PI) / (2 * M_PI)) * p_width);
        uint lon = static_cast<uint>((theta / M_PI) * p_height);
        // printf("lat: %d, lon: %d\n", lat, lon);

        uint latlong_index =
            lon * p_width * numComponents + lat * numComponents + k;

        if (sqrt((i - center) * (i - center) + (j - center) * (j - center)) <
            radius) {
          img_out[index] = img_in[latlong_index];
        }
      }
    }
  }

  WritePFM(image_out, width, height, numComponents, img_out);
}

void LoadAndSavePPM(const char *image_in, const char *image_out) {
  unsigned char *img_in = loadPNM(image_in, width, height, numComponents);
  unsigned char *img_out = new unsigned char[width * height * numComponents];

  for (uint i = 0; i < height; ++i) // height
  {
    for (uint j = 0; j < width; ++j) // width
    {
      for (uint k = 0; k < numComponents;
           ++k) // color channels - 3 for RGB images
      {
        uint index = i * width * numComponents + j * numComponents +
                     k; // index within the image

        img_out[index] = img_in[index]; // copy all color channels of each pixel
      }
    }
  }

  WritePNM(image_out, width, height, numComponents, img_out);
}

void LoadAndSavePFM(const char *image_in, const char *image_out) {
  float *img_in = loadPFM(image_in, width, height, numComponents);
  float *img_out = new float[width * height * numComponents];

  for (uint i = 0; i < height; ++i) // height
  {
    for (uint j = 0; j < width; ++j) // width
    {
      for (uint k = 0; k < numComponents;
           ++k) // color channels - 3 for RGB images
      {
        uint index = i * width * numComponents + j * numComponents +
                     k; // index within the image

        img_out[index] = img_in[index]; // copy all color channels of each pixel
      }
    }
  }

  WritePFM(image_out, width, height, numComponents, img_out);
}

void LoadPPMAndSavePFM(const char *image_in, const char *image_out) {
  unsigned char *img_in = loadPNM(image_in, width, height, numComponents);
  float *img_out = new float[width * height * numComponents];

  for (uint i = 0; i < height; ++i) // height
  {
    for (uint j = 0; j < width; ++j) // width
    {
      for (uint k = 0; k < numComponents;
           ++k) // color channels - 3 for RGB images
      {
        uint index = i * width * numComponents + j * numComponents +
                     k; // index within the image

        // typecast 0 - 255 values to the 0.0f -> 1.0f range
        img_out[index] = static_cast<float>(img_in[index]) /
                         255.0f; // typecast all color channels of each pixel
      }
    }
  }

  WritePFM(image_out, width, height, numComponents, img_out);
}

void LoadPFMAndSavePPM(const char *image_in, const char *image_out) {
  float *img_in = loadPFM(image_in, width, height, numComponents);
  unsigned char *img_out = new unsigned char[width * height * numComponents];

  for (uint i = 0; i < height; ++i) // height
  {
    for (uint j = 0; j < width; ++j) // width
    {
      for (uint k = 0; k < numComponents;
           ++k) // color channels - 3 for RGB images
      {
        uint index = i * width * numComponents + j * numComponents +
                     k; // index within the image

        // typecast 0.0f -> 1.0f values to the 0 - 255 range
        img_out[index] = static_cast<unsigned char>(
            (img_in[index] > 1.0f ? 1.0f : img_in[index]) *
            255.0f); // typecast all color channels of each pixel
      }
    }
  }

  WritePNM(image_out, width, height, numComponents, img_out);
}

int main(int argc, char **argv) {

  cerr << "main invoked: arguments - <image_out (.pfm)> " << endl;
  cerr << "main invoked: arguments - <image_in (.ppm)> <image_out (.ppm)> "
       << endl;

  int count = argc;

  if (count < 2) {
    cout << "Too few arguments! .... exiting." << endl;
    return 0;
  }

  if (count == 2) {
    CreateAndSavePFM(argv[1]); // Creates and saves a PFM
  } else if (count == 3) {
    //  LoadAndSavePPM(argv[1], argv[2]); //Loads and saves a PPM file
    //  LoadAndSavePFM(argv[1], argv[2]); //Loads and saves a PFM file
    //	LoadPPMAndSavePFM(argv[1], argv[2]); //Loads PPM and saves a PFM file
    // LoadPFMAndSavePPM(argv[1], argv[2]); // Loads PFM and saves a PPM file
    ReadProbeAndSavePFM(argv[1],
                        argv[2]); // loads light probe and saves a PFM file
    // gamma_correction(argv[2], argv[1], 1.7f);
  } else {
    cout << "Too many arguments! .... exiting." << endl;
    return 0;
  }

  return 0;
}
