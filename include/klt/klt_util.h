/*
 * klt_util.h
 *
 *  Created on: Feb 9, 2013
 *      Author: luzdora
 */

#ifndef KLT_UTIL_H_
#define KLT_UTIL_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct  {
  int ncols;
  int nrows;
  float *data;
}  _KLT_FloatImageRec, *_KLT_FloatImage;

_KLT_FloatImage _KLTCreateFloatImage(
  int ncols,
  int nrows);

void _KLTFreeFloatImage(
  _KLT_FloatImage);

void _KLTPrintSubFloatImage(
  _KLT_FloatImage floatimg,
  int x0, int y0,
  int width, int height);

void _KLTWriteFloatImageToPGM(
  _KLT_FloatImage img,
  char *filename);

/* for affine mapping */
void _KLTWriteAbsFloatImageToPGM(
  _KLT_FloatImage img,
  char *filename,float scale);

#ifdef __cplusplus
}
#endif

#endif /* KLT_UTIL_H_ */
