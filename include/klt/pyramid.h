/*
 * pyramid.h
 *
 *  Created on: Feb 9, 2013
 *      Author: luzdora
 */

#ifndef PYRAMID_H_
#define PYRAMID_H_

#include "klt_util.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct  {
  int subsampling;
  int nLevels;
  _KLT_FloatImage *img;
  int *ncols, *nrows;
}  _KLT_PyramidRec, *_KLT_Pyramid;


_KLT_Pyramid _KLTCreatePyramid(
  int ncols,
  int nrows,
  int subsampling,
  int nlevels);

void _KLTComputePyramid(
  _KLT_FloatImage floatimg,
  _KLT_Pyramid pyramid,
  float sigma_fact);

void _KLTFreePyramid(
  _KLT_Pyramid pyramid);

#ifdef __cplusplus
}
#endif


#endif /* PYRAMID_H_ */
