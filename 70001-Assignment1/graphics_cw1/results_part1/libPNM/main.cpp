/*//////////////////////////////////////////////////////////////////////////
  Author: Abhijeet Ghosh
  Year: 2013
//////////////////////////////////////////////////////////////////////////*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "PNM.h"
#include <cmath>
#include <algorithm>

// #define PI 3.14159265358979323
#define uint unsigned int

#include <iostream>

using namespace std;

unsigned int width;
unsigned int height;
unsigned int numComponents;

void CreateAndSavePFM(const char *image_out)
{
  width = 511; // set size of image to 511x511 pixels
  height = 511;
  numComponents = 3;

  float *img_out = new float[width * height * numComponents];

  for (uint i = 0; i < height; ++i) // height
  {
    for (uint j = 0; j < width; ++j) // width
    {
      for (uint k = 0; k < numComponents; ++k) // color channels - 3 for RGB images
      {
        uint index = i * width * numComponents + j * numComponents + k; // index within the image

        // set image to white
        img_out[index] = 1.0f; // RGB all set to white
      }
    }
  }

  writePFM(image_out, width, height, numComponents, img_out);
}

void LoadAndSavePPM(const char *image_in, const char *image_out)
{
  unsigned char *img_in = loadPPM(image_in, &width, &height, &numComponents);
  unsigned char *img_out = new unsigned char[width * height * numComponents];

  for (uint i = 0; i < height; ++i) // height
  {
    for (uint j = 0; j < width; ++j) // width
    {
      for (uint k = 0; k < numComponents; ++k) // color channels - 3 for RGB images
      {
        uint index = i * width * numComponents + j * numComponents + k; // index within the image

        img_out[index] = img_in[index]; // copy all color channels of each pixel
      }
    }
  }

  writePPM(image_out, width, height, numComponents, img_out);
}

void LoadAndSavePFM(const char *image_in, const char *image_out)
{
  float *img_in = loadPFM(image_in, &width, &height, &numComponents);
  float *img_out = new float[width * height * numComponents];

  for (uint i = 0; i < height; ++i) // height
  {
    for (uint j = 0; j < width; ++j) // width
    {
      for (uint k = 0; k < numComponents; ++k) // color channels - 3 for RGB images
      {
        uint index = i * width * numComponents + j * numComponents + k; // index within the image

        img_out[index] = img_in[index]; // copy all color channels of each pixel
      }
    }
  }

  writePFM(image_out, width, height, numComponents, img_out);
}

void LoadPPMAndSavePFM(const char *image_in, const char *image_out)
{
  unsigned char *img_in = loadPPM(image_in, &width, &height, &numComponents);
  float *img_out = new float[width * height * numComponents];

  for (uint i = 0; i < height; ++i) // height
  {
    for (uint j = 0; j < width; ++j) // width
    {
      for (uint k = 0; k < numComponents; ++k) // color channels - 3 for RGB images
      {
        uint index = i * width * numComponents + j * numComponents + k; // index within the image

        // typecast 0 - 255 values to the 0.0f -> 1.0f range
        img_out[index] = static_cast<float>(img_in[index]) / 255.0f; // typecast all color channels of each pixel
      }
    }
  }

  writePFM(image_out, width, height, numComponents, img_out);
}

void LoadPFMAndSavePPM(const char *image_in, const char *image_out)
{
  float *img_in = loadPFM(image_in, &width, &height, &numComponents);
  unsigned char *img_out = new unsigned char[width * height * numComponents];

  for (uint i = 0; i < height; ++i) // height
  {
    for (uint j = 0; j < width; ++j) // width
    {
      for (uint k = 0; k < numComponents; ++k) // color channels - 3 for RGB images
      {
        uint index = i * width * numComponents + j * numComponents + k; // index within the image

        // typecast 0.0f -> 1.0f values to the 0 - 255 range
        img_out[index] = static_cast<unsigned char>(img_in[index] * 255.0f); // typecast all color channels of each pixel
      }
    }
  }

  writePPM(image_out, width, height, numComponents, img_out);
}

// // Hat weighting function
long double weight(long double x)
{
  return 2 * max(0.0l, 0.5l - abs(x - 0.5l));
}

void GenerateHDR(std::vector<float *> images, std::vector<float> exposureTime, long double *img_out, unsigned int width, unsigned int height, unsigned int numComp)
{
  // Take sums
  long double *weighted_values = new long double[width * height * numComp];
  long double *weight_sum = new long double[width * height * numComp];
  for (int img_idx = 0; img_idx < images.size(); img_idx++)
  {
    for (uint i = 0; i < height; ++i) // height
    {
      for (uint j = 0; j < width; ++j) // width
      {
        bool valid = true;
        for (uint k = 0; k < numComp; ++k) // color channels - 3 for RGB images
        {
          uint index = i * width * numComp + j * numComp + k; // index within the image
          if (images[img_idx][index] < 0.005 || images[img_idx][index] > 0.92)
            valid = false;
        }

        // Skip if invalid pixels
        if (!valid)
          continue;

        for (uint k = 0; k < numComp; ++k) // color channels - 3 for RGB images
        {
          uint index = i * width * numComp + j * numComp + k; // index within the image

          weighted_values[index] += weight(images[img_idx][index]) 
          * logl((1.0L / ((long double)exposureTime[img_idx])) * images[img_idx][index]);

          weight_sum[index] += weight(images[img_idx][index]);
        }
      }
    }
  }

  long double max_value = 0.0L;
  long double min_value = 1.0L;

  // Compute final HDR image
  for (int img_idx = 0; img_idx < images.size(); img_idx++)
  {
    for (uint i = 0; i < height; ++i) // height
    {
      for (uint j = 0; j < width; ++j) // width
      {
        for (uint k = 0; k < numComp; ++k) // color channels - 3 for RGB images
        {
          uint index = i * width * numComp + j * numComp + k; // index within the image
          img_out[index] = expl(weighted_values[index] / weight_sum[index]);


          if (img_out[index] > max_value)
          {
            max_value = img_out[index];
          }
          if (img_out[index] < min_value)
          {
            min_value = img_out[index];
          }
        }
      }
    }
  }
  // Write to outpath provided from command line
  cout << "Max: " << max_value << endl;
  cout << "Min: " << min_value << endl;
  cout << "Dynamic Range: " << max_value / min_value << endl;

  delete[] weighted_values;
  delete[] weight_sum;

  return;
}

void PFMsavePPM(long double *img_in, int width, int height, int numComponents, const char *image_out)
{
  unsigned char *img_out = new unsigned char[width * height * numComponents];

  for (uint i = 0; i < height; ++i) // height
  {
    for (uint j = 0; j < width; ++j) // width
    {
      for (uint k = 0; k < numComponents; ++k) // color channels - 3 for RGB images
      {
        uint index = i * width * numComponents + j * numComponents + k; // index within the image

        // typecast 0.0f -> 1.0f values to the 0 - 255 range
        img_out[index] = static_cast<unsigned char>(img_in[index] * 255.0f); // typecast all color channels of each pixel
      }
    }
  }

  writePPM(image_out, width, height, numComponents, img_out);
}

void simpleToneMapper(long double *img_in, unsigned int width, unsigned int height, unsigned int numComponents, long double *img_out)
{

  // Get max and min values
  float max_value = 0.0f;
  for (uint i = 0; i < height; ++i) // height
  {
    for (uint j = 0; j < width; ++j) // width
    {
      for (uint k = 0; k < numComponents; ++k) // color channels - 3 for RGB images
      {
        uint index = i * width * numComponents + j * numComponents + k; // index within the image

        if (img_in[index] > max_value)
        {
          max_value = img_in[index];
        }
      }
    }
  }

  for (uint i = 0; i < height; ++i) // height
  {
    for (uint j = 0; j < width; ++j) // width
    {
      for (uint k = 0; k < numComponents; ++k) // color channels - 3 for RGB images
      {
        uint index = i * width * numComponents + j * numComponents + k; // index within the image

        // Scales from 0 -> Max to 0 -> 1
        img_out[index] = (img_in[index]) / (max_value);
      }
    }
  }
}

void exposureMapper(long double *img_in, unsigned int width, unsigned int height, unsigned int numComponents, long double *img_out, int stops)
{
  for (uint i = 0; i < height; ++i) // height
  {
    for (uint j = 0; j < width; ++j) // width
    {
      for (uint k = 0; k < numComponents; ++k) // color channels - 3 for RGB images
      {
        uint index = i * width * numComponents + j * numComponents + k; // index within the image

        img_out[index] = img_in[index] * (long double)powl(2.0L, (long double)stops);

        // Was having issues with the c++ min function for clamping so reimplimented it myself
        if (img_out[index] > 1.0L)
        {
          img_out[index] = 1.0L;
        }
      }
    }
  }
}

void gammaMapper(long double *img_in, unsigned int width, unsigned int height, unsigned int numComponents, long double *img_out, float gamma)
{

  for (uint i = 0; i < height; ++i) // height
  {
    for (uint j = 0; j < width; ++j) // width
    {
      for (uint k = 0; k < numComponents; ++k) // color channels - 3 for RGB images
      {
        uint index = i * width * numComponents + j * numComponents + k; // index within the image

        img_out[index] = powl(img_in[index], 1.0L / gamma); // typecast all color channels of each pixel
      }
    }
  }
}

int main(int argc, char **argv)
{
  string outpath = argv[1];

  std::vector<float *> images;
  std::vector<unsigned int> widths, heights, numComponents;
  std::vector<float> exposureTime;

  for (int i = 1; i <= 7; ++i)
  {
    std::string filename = "Office/Office" + std::to_string(i) + ".pfm";

    unsigned int width, height, numComp;
    float *img = loadPFM(filename.c_str(), &width, &height, &numComp);

    // Need to check they are in the  correct range of pixel colours
    if (img != nullptr)
    {
      images.push_back(img);

      for (uint i = 0; i < height; ++i) // height
      {
        for (uint j = 0; j < width; ++j) // width
        {
          for (uint k = 0; k < numComp; ++k) // color channels - 3 for RGB images
          {
            uint index = i * width * numComp + j * numComp + k; // index within the image
            if (img[index] < 0 || img[index] > 1)
              cout << "pixel value:" << img[index] << endl;
          }
        }
      }

      // Metadata
      widths.push_back(width);
      heights.push_back(height);
      numComponents.push_back(numComp);

      // exposureTime
      exposureTime.push_back(powf(2, (2 * (i - 1))));
    }
    else
    {
      std::cerr << "Failed to load image: " << filename << std::endl;
    }
  }

  // verify that the widths are the same for all images
  for (unsigned int i = 1; i < widths.size(); ++i)
  {
    if (widths[i] != widths[0])
    {
      std::cerr << "Widths are not the same for all images" << std::endl;
      return -1;
    }
    if (heights[i] != heights[0])
    {
      std::cerr << "Heights are not the same for all images" << std::endl;
      return -1;
    }
    if (numComponents[i] != numComponents[0])
    {
      std::cerr << "Compontent count is not the same for all images" << std::endl;
      return -1;
    }
  }
  unsigned int numComp = numComponents[0];
  unsigned int width = widths[0];
  unsigned int height = heights[0];

  long double *hdr_img = new long double[width * height * numComp];
  GenerateHDR(images, exposureTime, hdr_img, width, height, numComp);

  // Apply simple tone mapping
  long double *simple_hdr_img = new long double[width * height * numComp];
  simpleToneMapper(hdr_img, width, height, numComp, simple_hdr_img);
  PFMsavePPM(simple_hdr_img, width, height, numComp, (outpath + "/simple.ppm").c_str());

  // Apply Exposure Mapping
  long double *exposure_hdr_img_1 = new long double[width * height * numComp];
  long double *exposure_hdr_img_2 = new long double[width * height * numComp];
  long double *exposure_hdr_img_3 = new long double[width * height * numComp];
  long double *exposure_hdr_img_4 = new long double[width * height * numComp];
  long double *exposure_hdr_img_5 = new long double[width * height * numComp];
  long double *exposure_hdr_img_6 = new long double[width * height * numComp];
  long double *exposure_hdr_img_7 = new long double[width * height * numComp];
  long double *exposure_hdr_img_8 = new long double[width * height * numComp];
  long double *exposure_hdr_img_9 = new long double[width * height * numComp];
  long double *exposure_hdr_img_10 = new long double[width * height * numComp];
  

  exposureMapper(simple_hdr_img, width, height, numComp, exposure_hdr_img_1, 1);
  exposureMapper(simple_hdr_img, width, height, numComp, exposure_hdr_img_2, 2);
  exposureMapper(simple_hdr_img, width, height, numComp, exposure_hdr_img_3, 3);
  exposureMapper(simple_hdr_img, width, height, numComp, exposure_hdr_img_4, 4);
  exposureMapper(simple_hdr_img, width, height, numComp, exposure_hdr_img_5, 5);
  exposureMapper(simple_hdr_img, width, height, numComp, exposure_hdr_img_6, 6);
  exposureMapper(simple_hdr_img, width, height, numComp, exposure_hdr_img_7, 7);
  exposureMapper(simple_hdr_img, width, height, numComp, exposure_hdr_img_8, 8);
  exposureMapper(simple_hdr_img, width, height, numComp, exposure_hdr_img_9, 9);
  exposureMapper(simple_hdr_img, width, height, numComp, exposure_hdr_img_10, 10);


  PFMsavePPM(exposure_hdr_img_1, width, height, numComp, (outpath + "/exposure_1.ppm").c_str());
  PFMsavePPM(exposure_hdr_img_2, width, height, numComp, (outpath + "/exposure_2.ppm").c_str());
  PFMsavePPM(exposure_hdr_img_3, width, height, numComp, (outpath + "/exposure_3.ppm").c_str());
  PFMsavePPM(exposure_hdr_img_4, width, height, numComp, (outpath + "/exposure_4.ppm").c_str());
  PFMsavePPM(exposure_hdr_img_5, width, height, numComp, (outpath + "/exposure_5.ppm").c_str());
  PFMsavePPM(exposure_hdr_img_6, width, height, numComp, (outpath + "/exposure_6.ppm").c_str());
  PFMsavePPM(exposure_hdr_img_7, width, height, numComp, (outpath + "/exposure_7.ppm").c_str());
  PFMsavePPM(exposure_hdr_img_8, width, height, numComp, (outpath + "/exposure_8.ppm").c_str());
  PFMsavePPM(exposure_hdr_img_9, width, height, numComp, (outpath + "/exposure_9.ppm").c_str());
  PFMsavePPM(exposure_hdr_img_10, width, height, numComp, (outpath + "/exposure_10.ppm").c_str());


  // Apply Gamma Mapping
  long double *gamma_hdr_img_1_5 = new long double[width * height * numComp];
  long double *gamma_hdr_img_1_8 = new long double[width * height * numComp];
  long double *gamma_hdr_img_2_2 = new long double[width * height * numComp];
  long double *gamma_hdr_img_2_5 = new long double[width * height * numComp];
  long double *gamma_hdr_img_3_0 = new long double[width * height * numComp];

  gammaMapper(exposure_hdr_img_1, width, height, numComp, gamma_hdr_img_1_5 , 1.5);
  gammaMapper(exposure_hdr_img_1, width, height, numComp, gamma_hdr_img_1_8 , 1.8);
  gammaMapper(exposure_hdr_img_1, width, height, numComp, gamma_hdr_img_2_2 , 2.2);
  gammaMapper(exposure_hdr_img_1, width, height, numComp, gamma_hdr_img_2_5 , 2.5);
  gammaMapper(exposure_hdr_img_1, width, height, numComp, gamma_hdr_img_3_0 , 3.0);

  PFMsavePPM(gamma_hdr_img_1_5, width, height, numComp, (outpath + "/gamma_1_5.ppm").c_str());
  PFMsavePPM(gamma_hdr_img_1_8, width, height, numComp, (outpath + "/gamma_1_8.ppm").c_str());
  PFMsavePPM(gamma_hdr_img_2_2, width, height, numComp, (outpath + "/gamma_2_2.ppm").c_str());
  PFMsavePPM(gamma_hdr_img_2_5, width, height, numComp, (outpath + "/gamma_2_5.ppm").c_str());
  PFMsavePPM(gamma_hdr_img_3_0, width, height, numComp, (outpath + "/gamma_3_0.ppm").c_str());

  // Remember to free the memory for each loaded image when done
  for (float *img : images)
  {
    delete[] img;
  }

  return 0;
}
