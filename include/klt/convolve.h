/*
 * convolve.h
 *
 *  Created on: Feb 12, 2013
 *      Author: luzdora
 */

#ifndef CONVOLVE_H_
#define CONVOLVE_H_

#include "klt.h"
//#include "klt_util.h"
//#include "base.h"
#include "error.h"

#ifdef __cplusplus
extern "C" {
#endif

void _KLTToFloatImage(
  KLT_PixelType *img,
  int ncols, int nrows,
  _KLT_FloatImage floatimg);

void _KLTComputeGradients(
  _KLT_FloatImage img,
  float sigma,
  _KLT_FloatImage gradx,
  _KLT_FloatImage grady);

void _KLTGetKernelWidths(
  float sigma,
  int *gauss_width,
  int *gaussderiv_width);

void _KLTComputeSmoothedImage(
  _KLT_FloatImage img,
  float sigma,
  _KLT_FloatImage smooth);

#ifdef __cplusplus
}
#endif


#endif /* CONVOLVE_H_ */
