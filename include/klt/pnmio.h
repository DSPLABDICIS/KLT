/*
 * pnmio.h
 *
 *  Created on: Feb 9, 2013
 *      Author: luzdora
 */

#ifndef PNMIO_H_
#define PNMIO_H_

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif


/**********
 * With pgmReadFile and pgmRead, setting img to NULL causes memory
 * to be allocated
 */

/**********
 * used for reading from/writing to files
 */
unsigned char* pgmReadFile(
  char *fname,
  unsigned char *img,
  int *ncols,
  int *nrows);
void pgmWriteFile(
  char *fname,
  unsigned char *img,
  int ncols,
  int nrows);
void ppmWriteFileRGB(
  char *fname,
  unsigned char *redimg,
  unsigned char *greenimg,
  unsigned char *blueimg,
  int ncols,
  int nrows);

/**********
 * used for communicating with stdin and stdout
 */
unsigned char* pgmRead(
  FILE *fp,
  unsigned char *img,
  int *ncols, int *nrows);
void pgmWrite(
  FILE *fp,
  unsigned char *img,
  int ncols,
  int nrows);
void ppmWrite(
  FILE *fp,
  unsigned char *redimg,
  unsigned char *greenimg,
  unsigned char *blueimg,
  int ncols,
  int nrows);
#ifdef __cplusplus
}
#endif

#endif /* PNMIO_H_ */
