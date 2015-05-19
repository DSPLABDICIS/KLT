/*
 * kltTest.cpp
 *
 *  Created on: Feb 9, 2013
 *      Author: luzdora
 */

/** This test function is based on the example 3 of the KLT lib, but adds a lot of dynamic parameters, and enable
fine tunings such as in the example 5 of the KLT lib. It should enable to do a lot of tests in a tcl shell
without having to recompile, or to make a long tcl macro ...
 */

#include "kltTest.hpp"
//#include "kernel/timingTools.hpp"
//namespace klt {
using namespace klt;
void initTrackContext(KLT_TrackingContext tc,
		int affineConsistencyCheck,
		int window_size,
		int nPyramidLevels,
		int subsampling,
		int mindist,
		double max_residue,
		int max_iterations,
		double min_determinant)
{

	///**** Tuning

	if (affineConsistencyCheck < 0) affineConsistencyCheck = -1;
	if (window_size < 0) window_size = 11;
	if (nPyramidLevels < 0) nPyramidLevels = 3;
	if (subsampling < 0) subsampling = 2;
	if (mindist < 0) mindist = 10;
	if (max_residue < 0) max_residue = 24;
	if (max_iterations < 0) max_iterations = 10;
	if (min_determinant < 0) min_determinant = 0.01;

	tc->affineConsistencyCheck = affineConsistencyCheck; // -1
	tc->window_width = window_size; //en realidad 11  tiene 7 desde la libreria
	tc->window_height = window_size; //en realidad 11  tiene 7 desde la libreria
	tc->subsampling = subsampling; //2
	//    KLTUpdateTCBorder(tc);
	tc->nPyramidLevels = nPyramidLevels; //5
	tc->mindist = mindist; // 10
	tc->max_residue = max_residue; //24
	tc->max_iterations = max_iterations;//10
	tc->min_determinant = min_determinant; //0.01
	tc->sequentialMode = FALSE; // to improve speed (do not recompute gradient of i in (i-1->i),(i,i+1)
	tc->writeInternalImages = FALSE;
	tc->smoothBeforeSelecting = TRUE;   // TRUE
	//KLTStopSequentialMode(tc); // don't do that if no tracking has been done because it will try to free non allocated mem
	KLTUpdateTCBorder(tc);
	//   std::cout << "tc->w_width " <<  tc->window_width << " tc->wheight "<<  tc->window_height <<::std::endl;
}

int countMovingTrackingFeatures( KLT_FeatureList fl)
{
	return KLTCountMovingFeatures(fl);

}

int countMovingFeaturesTotal(KLT_FeatureList fl, int modu)
{
  int count = 0;
  int i;

  for (i = 0 ; i < fl->nFeatures ; i++)
    if (fl->feature[i]->trail > modu-1)
      count++;
  return count;
}

//In this function sMprob is the size of the Probability Mask that has
// the same value during all the process
uchar* reserveOccupationMask (int ncols, int nrows)
{
	uchar *featuremap;

	featuremap = (uchar *) malloc(ncols * nrows * sizeof(uchar));
	//  memset(featuremap, 1, ncols*nrows);
	return featuremap;
}


// Pyramidal mask
/*
	  void initialiseProbaMask (uchar * maskProb, int size, int mdist)
	  {
	    // initialise probability mask
	    int i, j, vi, vj, incp;
	    incp = (255-127)/mdist;

	    for (i = 0 ; i < size; i++){
	      vi = (i<mdist)? incp*i:(2*mdist-i)*incp;
	      for (j=0; j < size; j++) {
		vj = (j<mdist)? incp*j:(2*mdist-j)*incp;
		maskProb[j*size+i]= (vi<vj)? vi :vj;
		//	printf (" M(%d,%d)=%d ", i, j, maskProb[j*size+i]);
	      }
	      //   printf ("\n ");
	    }

	  }
 */

// Gaussian  mask

void initialiseProbaMask (unsigned char * maskProb, int sizeM, int mdist)
{
	// initialise probability mask
	int i, j, i1, j1;
	// Parameters of the bivariate normal distribution
	float mux = 0, muy = 0, sigx = 6, sigy = 6;
	int fair =  floor(255/2);
	double aux;

	for ( i=0, i1 =-1*mdist ; i < sizeM; i++, i1++){
		for (j=0, j1=-1*mdist; j < sizeM; j++, j1++) {
			aux = (-0.5)*( ((i1-mux)/sigx)*((i1-mux)/sigx) + ((j1-muy)/sigy)*((j1-muy)/sigy));
			maskProb[j*sizeM+i] = fair*exp(aux);
	//   printf (" M(%d,%d)=%d ", i, j, maskProb[j*size+i]);
		}
	//  printf ("\n ");
	}

}

void initialiseOccupationMask (uchar * featuremap, uchar *maskP, int cols, int rows)
{
//	int i, j, mindist = 9;
//	int x1=60, y1=10;
	memset(featuremap, 127, cols*rows);

	/*   for(j=0; j <= (int) (310/mindist); j++)
	      for (i=0; i<= (int) (190/mindist) ; i++)
		_fillFeaturemap(x1 + (mindist*j), y1 + (mindist*i), featuremap, maskP, mindist, cols, rows);

	 */
}

void setOccupationZoneInMask (uchar * mask, uchar *maskP,  KLT_FeatureList fl, int ncols, int nrows, int mindist, int flag)
{
	int indx, x, y;
	//  mindist--;

	if (flag)
	{ // Solo los puntos que no se movieron seran colocados en la mascara
		for (indx = 0 ; indx < fl->nFeatures ; indx++){
			//	  if (fl->feature[indx]->val == 1)  {
			if ( (fl->feature[indx]->val == 0)  &&  ((fabs(fl->feature[indx]->vx) < LIM ) || ( fabs(fl->feature[indx]->vy) < LIM))) {
				x   = (int) fl->feature[indx]->x;
				y   = (int) fl->feature[indx]->y;
				_fillFeaturemap(x, y, mask, maskP, mindist, ncols, nrows);
				fl->feature[indx]->val = -1; //para que sea eliminado de la lista
			}
		}
	}
	else
	{ // todos son colocados en la mascara y eliminados de la lista
		for (indx = 0 ; indx < fl->nFeatures ; indx++){
			x   = (int) fl->feature[indx]->x;
			y   = (int) fl->feature[indx]->y;
			_fillFeaturemap(x, y, mask, maskP, mindist, ncols, nrows);
			fl->feature[indx]->val = -1; //para que sea eliminado de la lista
		}
	}
}

void freeOccupationMask(uchar *featuremap)
{
	free(featuremap);
}

void KLTwriteImageToImg (cv::Mat &imGrey, cv::Mat & imColor)
{
	WriteImageToImg(imGrey,imColor);
}

//void featuresToWorldCoordinate (KLT_FeatureList fl, jafar::camera::CameraPinhole  sensor, jafar::geom::T3DEuler cameraToOrigin, int imaH)
//{
//	boost_incl::vec point(2);
//	boost_incl::vec pointCF(3);
//	boost_incl::vec pointW(3);
//	boost_incl::vec pointCam(3);
//	boost_incl::vec pWCam(3);
//
//	pointCam(0) = 0;
//	pointCam(1) = 0;
//	pointCam(2) = 0;
//	int i;
//	// Camera World Position
//	geom::t3d::pointFromFrame(cameraToOrigin, pointCam, pWCam);
//
//	for (i = 0; i <fl->nFeatures ; i++)
//	{
//		if (fl->feature[i]->val >= 0)  {
//			point(0)  =  fl->feature[i]->x;
//			point(1)  =  fl->feature[i]->y;
//
//			sensor.imageToCameraFrame(point, pointCF);
//			// Infinito projection : p0 (point)
//
//			/*	  if ( fl->feature[i]->y < imaH/3 ){
//		    pointCF(0)*=20;
//		    pointCF(1)*=20;
//		    pointCF(2) =20;
//		    } */
//
//			geom::t3d::pointFromFrame(cameraToOrigin, pointCF, pointW);
//
//			fl->feature[i]->xem = pointW(0);
//			fl->feature[i]->yem = pointW(1);
//			//	  fl->feature[i]->yem = 100;
//			//	  fl->feature[i]->zem = 10000;
//			fl->feature[i]->zem = pointW(2);
//
//			// Floor projection : p0 (point) , p1 (camera)
//			if ( fl->feature[i]->y > imaH/2 ){ // al suelo
//				fl->feature[i]->xpro = ( (-1*pointW(2))*(pointW(0)-pWCam(0))/(pointW(2)-pWCam(2)) )+pointW(0);
//				fl->feature[i]->ypro = ( (-1*pointW(2))*(pointW(1)-pWCam(1))/(pointW(2)-pWCam(2)) )+pointW(1);
//				fl->feature[i]->zpro = 0.0;
//			} //end projection
//			else
//			{
//				fl->feature[i]->xpro = -1.0;
//				fl->feature[i]->ypro = -1.0;
//				fl->feature[i]->zpro = -1.0;
//			}
//		} // en if
//
//	} // end for
//
//} // end function


//void featuresToNextImage (KLT_FeatureList fl, jafar::camera::CameraPinhole & sensor, jafar::geom::T3DEuler cameraToOrigin)
//{
//	boost_incl::vec point(2);
//	boost_incl::vec pointWF(3);
//	boost_incl::vec pointW(3);
//	int i;
//
//	for (i = 0; i <fl->nFeatures ; i++)
//	{
//		if (fl->feature[i]->val >= 0)  {
//			if ( fl->feature[i]->xpro*fl->feature[i]->ypro*fl->feature[i]->zpro == -1 ){
//				pointW(0)   =  fl->feature[i]->xem;
//				pointW(1)   =  fl->feature[i]->yem;
//				pointW(2)   =  fl->feature[i]->zem;
//			}
//			else
//			{
//				pointW(0)   =  fl->feature[i]->xpro;
//				pointW(1)   =  fl->feature[i]->ypro;
//				pointW(2)   =  fl->feature[i]->zpro;
//			}
//
//			geom::t3d::pointToFrame(cameraToOrigin, pointW, pointWF);
//			sensor.project(pointWF, point);
//
//			fl->feature[i]->xem = point(0);
//			fl->feature[i]->yem = point(1);
//			fl->feature[i]->zem = pointW(2);
//
//		} // en if
//		else
//		{
//			fl->feature[i]->xem = -1.0;
//			fl->feature[i]->yem = -1.0;
//			fl->feature[i]->zem = -1.0;
//		}
//
//	} // end for
//
//} // end function


void projectAlSuelo( boost_incl::vec3 & pointW,  boost_incl::vec3 & pWCam,  boost_incl::vec3 & pS )
{
	pS(0) = ( (-1*pointW(2))*(pointW(0)-pWCam(0))/(pointW(2)-pWCam(2)) )+pointW(0);
	pS(1) = ( (-1*pointW(2))*(pointW(1)-pWCam(1))/(pointW(2)-pWCam(2)) )+pointW(1);
	pS(2) = 0.0;
}
//} // end namespace klt

