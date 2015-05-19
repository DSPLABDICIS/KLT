/*
 * klt_.cpp
 *
 *  Created on: Feb 9, 2013
 *      Author: luzdora
 */

/*********************************************************************
 * klt.cpp
 * interface jafar pour l'impl�mentation de Birchfield du
 * Kanade-Lucas-Tomasi tracker
 * http://www.ces.clemson.edu/~stb/klt/
 * attention la fonctionnalit� de ne choisir les features que dans une fen�tre a �t� ajout�e,
 * un patch devrait �tre dispo pour permettre de mettre � jour la lib
 * contact : cyril.roussillon@laas.fr
 *********************************************************************/

#include <iostream>
#include "klt.hpp"
//#include "kernel/jafarException.hpp"
#include "error.h"

namespace klt {
using namespace klt;

enum rect {x1p=0,y1p,x2p,y2p};

inline double sqr(double x) { return x*x; }
// klt::CompatTreatment treatment;
//  extern int KLT_verbose;

/***************************************************************************
          ################################################################################
           CalcVelocity components X and Y (pix/sec)

          ################################################################################
 *****************************************************************************/

void calcVelocity(KLT_FeatureList featurelist,
		KLT_FeatureTable ft,
		int currIma, double t, int flag, bool _3d)
{
	float dx, dy, dz;
	int nfeat;

	if (flag){
		for ( nfeat=0; nfeat < featurelist->nFeatures; nfeat++)
		{
			if ((featurelist->feature[nfeat]->val == KLT_TRACKED) && (ft->feature[nfeat][currIma-1]->val >= KLT_TRACKED))  {
				//   if ( featurelist->feature[nfeat]->y > 480/2 ){
				/*
		dx = featurelist->feature[nfeat]->x - featurelist->feature[nfeat]->xem;
		dy = featurelist->feature[nfeat]->y - featurelist->feature[nfeat]->yem;
		}
		  else
		  { */
				if(_3d){
					if(featurelist->feature[nfeat]->_3Dlost == 1){
						if(ft->feature[nfeat][currIma-1]->_3Dlost == 0){        //3D image Derp, give 1 chance and set a constant speed
							dx = ft->feature[nfeat][currIma-1]->vx * t;
							dy = ft->feature[nfeat][currIma-1]->vy * t;
							dz = ft->feature[nfeat][currIma-1]->vz * t;
							featurelist->feature[nfeat]->_3Dx = ft->feature[nfeat][currIma-1]->_3Dx + dx;   // 3D Position estimation in with constant speed
							featurelist->feature[nfeat]->_3Dy = ft->feature[nfeat][currIma-1]->_3Dy + dy;
							featurelist->feature[nfeat]->_3Dz = ft->feature[nfeat][currIma-1]->_3Dz + dz;
						}
						else{													//3D If derp for a second time, NO GO...
							dx = 0;
							dy = 0;
							dz = 0;
							featurelist->feature[nfeat]->val = KLT_NOT_FOUND;
						}
					}
					else{
						dx=featurelist->feature[nfeat]->_3Dx - ft->feature[nfeat][currIma-1]->_3Dx;
						dy=featurelist->feature[nfeat]->_3Dy - ft->feature[nfeat][currIma-1]->_3Dy;
						dz=featurelist->feature[nfeat]->_3Dz - ft->feature[nfeat][currIma-1]->_3Dz;
					}
				}
				else{
					dx=featurelist->feature[nfeat]->x - ft->feature[nfeat][currIma-1]->x;
					dy=featurelist->feature[nfeat]->y - ft->feature[nfeat][currIma-1]->y;
					dz = 0;
				}
				//	  } */
				featurelist->feature[nfeat]->vx = dx/t;
				featurelist->feature[nfeat]->vy = dy/t;
				featurelist->feature[nfeat]->vy = dz/t;

				if ((fabs(featurelist->feature[nfeat]->vx) > LIM ) || ( fabs(featurelist->feature[nfeat]->vy) > LIM) /*|| (fabs(featurelist->feature[nfeat]->vz) > LIM )*/)
					//  if (fabs(featurelist->feature[nfeat]->vx) > LIM )
					featurelist->feature[nfeat]->trail++;
				//   std::cout << "feat " << nfeat  << " x " << featurelist->feature[nfeat]->x << " y " <<featurelist->feature[nfeat]->y << " currIma: "<< currIma << " trail " << featurelist->feature[nfeat]->trail << std::endl;
				else
					featurelist->feature[nfeat]->trail = 0;  // just it is not a moving point
			}  // In the case of it was a moving point before

			else
				featurelist->feature[nfeat]->trail = 0;  // just it is not a moving point
			//	  std::cout << "NO MOV  feat " << nfeat  << " x " << featurelist->feature[nfeat]->x << " y " <<featurelist->feature[nfeat]->y << " currIma: "<< currIma << " trail " << featurelist->feature[nfeat]->trail << std::endl;

		}
	}
	else
	{
		for ( nfeat=0; nfeat < featurelist->nFeatures; nfeat++)
		{
			if ((featurelist->feature[nfeat]->val == KLT_TRACKED) && (ft->feature[nfeat][currIma-1]->val >= KLT_TRACKED)){
				if(_3d){
					if(featurelist->feature[nfeat]->_3Dlost == 1){
						if(ft->feature[nfeat][currIma-1]->_3Dlost == 0){        //3D image Derp, give 1 chance and set a constant speed
							dx = ft->feature[nfeat][currIma-1]->vx * t;
							dy = ft->feature[nfeat][currIma-1]->vy * t;
							dz = ft->feature[nfeat][currIma-1]->vz * t;
							featurelist->feature[nfeat]->_3Dx = ft->feature[nfeat][currIma-1]->_3Dx + dx;   // 3D Position estimation in with constant speed
							featurelist->feature[nfeat]->_3Dy = ft->feature[nfeat][currIma-1]->_3Dy + dy;
							featurelist->feature[nfeat]->_3Dz = ft->feature[nfeat][currIma-1]->_3Dz + dz;
						}
						else{													//3D If derp for a second time, NO GO...
							dx = 0;
							dy = 0;
							dz = 0;
							featurelist->feature[nfeat]->val = KLT_NOT_FOUND;
						}
					}
					else{
						dx=featurelist->feature[nfeat]->_3Dx - ft->feature[nfeat][currIma-1]->_3Dx;
						dy=featurelist->feature[nfeat]->_3Dy - ft->feature[nfeat][currIma-1]->_3Dy;
						dz=featurelist->feature[nfeat]->_3Dz - ft->feature[nfeat][currIma-1]->_3Dz;
					}
				}
				else{
					dx=featurelist->feature[nfeat]->x - ft->feature[nfeat][currIma-1]->x;
					dy=featurelist->feature[nfeat]->y - ft->feature[nfeat][currIma-1]->y;
					dz = 0;
				}
				featurelist->feature[nfeat]->vx   = dx/t;
				featurelist->feature[nfeat]->vy   = dy/t;
				featurelist->feature[nfeat]->vz   = dz/t;

				if ((fabs(featurelist->feature[nfeat]->vx) > LIM ) || ( fabs(featurelist->feature[nfeat]->vy) > LIM) || (fabs(featurelist->feature[nfeat]->vz) > LIM ))
					//if (fabs(featurelist->feature[nfeat]->vx) > LIM )
					featurelist->feature[nfeat]->trail++;
				//   std::cout << "feat " << nfeat  << " x " << featurelist->feature[nfeat]->x << " y " <<featurelist->feature[nfeat]->y << " currIma: "<< currIma << " trail " << featurelist->feature[nfeat]->trail << std::endl;
				else
					featurelist->feature[nfeat]->trail = 0;  // just it is not a moving point

			}  // In the case of it was a moving point before

			else
				featurelist->feature[nfeat]->trail = 0;  // just it is not a moving point
			//	  std::cout << "NO MOV  feat " << nfeat  << " x " << featurelist->feature[nfeat]->x << " y " <<featurelist->feature[nfeat]->y << " currIma: "<< currIma << " trail " << featurelist->feature[nfeat]->trail << std::endl;
		}
	} // end else
} // end of function


void Klt::SelectGoodFeatures(KLT_TrackingContext tc,
		cv::Mat &jImg,
		KLT_FeatureList fl,
		double timeStamp, int r, uchar *featmap, uchar *maskProb)
{
	if (timeStamp < 0) timeStamp = 0;
	KLTSelectGoodFeatures(tc, jImg.data, jImg.cols, jImg.rows, fl, r, featmap, maskProb);
	for (int i=0; i < fl->nFeatures; i++){
		fl->feature[i]->vx=0.0;
		fl->feature[i]->vy=0.0;
		fl->feature[i]->vz=0.0;
	}
	frame5 = 0;

	//   StoreFeatureList(fl, ft, curFrame);
	//   StoreFeatureList(fl, ft5, frame5);
	// Found features were already added to the probability map
}

// Extrae puntos moviles y los estaticos los activa en la mascara comoimage::Image
// zona no interesante para busqueda de nuevos puntos

int Klt::ExtractMovingPoints(boost_incl::mat &Xp, int fima, int mod, KLT_FeatureList fl, bool _3D)
{
	int i, j=0, l, aux, nima, dim, extradim=0, count=0;
	dim = Xp.size2();

	// Count the number of total points
	for ( i=0; i<fl->nFeatures; i++)
		if (fl->feature[i]->trail>=mod-1)
			count+=fl->feature[i]->trail;
	//	 std::cout << " feat num " << i << " trail " << fl->feature[i]->trail << std::endl;

	// numFeat = count;
	//std::cout << "count " << count << std::endl;

	Xp.resize(count, dim, 0);

	// Extract the feature points in a jblas matrix
	nima = fima+mod;
	for (l=fima; l< nima; l++){
		if (l==ft->nFrames)
		{
			aux = (fima+mod)-l;
			nima = aux;
			l = 0;
		}

		if(_3D)
			for ( i=0; i<fl->nFeatures; i++)
			{
				if (fl->feature[i]->trail>=mod-1)
					if ( ft->feature[i][l]->trail > 0 )
					{
						Xp(j,0)  = ft->feature[i][l]->x;
						Xp(j,1)  = ft->feature[i][l]->y;
						Xp(j,2)  = ft->feature[i][l]->vx;
						Xp(j,3)  = ft->feature[i][l]->vy;
						Xp(j,4)  = 0;
						Xp(j,5)  = ft->feature[i][l]->auxval;
						Xp(j,6)  = ft->feature[i][l]->_3Dx;
						Xp(j,7)  = ft->feature[i][l]->_3Dy;
						Xp(j,8)  = ft->feature[i][l]->_3Dz;
						Xp(j,9) = ft->feature[i][l]->vz;
						//	std::cout << " feat num " << i << "En imagen " << l ;
						//	std::cout << " En cartesianas " << 	Xp(j,0) << " " << Xp(j,1) << " " << Xp(j,2) << " " << Xp(j,3) << "\n";
						//std::cout<<Xp(j,0)<<" "<<Xp(j,1)<<" "<<Xp(j,2)<<" "<<Xp(j,3)<<" "<<Xp(j,6)<<" "<<Xp(j,7)<<" "<<Xp(j,8)<<" "<<Xp(j,9)<<std::endl;
						j++;

					}
				//	  else
				//   ft->feature[i][l]->val = 1; //Bandera para que utilise el lugar
			}
		else
			for ( i=0; i<fl->nFeatures; i++)
			{
				if (fl->feature[i]->trail>=mod-1)
					if ( ft->feature[i][l]->trail > 0 )
					{
						Xp(j,0) = ft->feature[i][l]->x;
						Xp(j,1) = ft->feature[i][l]->y;
						Xp(j,2) = ft->feature[i][l]->vx;
						Xp(j,3) = ft->feature[i][l]->vy;
						Xp(j,4) = 0;
						Xp(j,5) = ft->feature[i][l]->auxval;
						//	std::cout << " feat num " << i << "En imagen " << l ;
						//	std::cout << " En cartesianas " << 	Xp(j,0) << " " << Xp(j,1) << " " << Xp(j,2) << " " << Xp(j,3) << "\n";
						j++;

					}
				//	  else
				//   ft->feature[i][l]->val = 1; //Bandera para que utilise el lugar
			}
		//    std::cout << std::endl;

	}
	//std::cout << "Size of data to cluster " << j << std::endl;
	//Before to leave, reset value in feature list for new time of trail
	for ( i=0; i<fl->nFeatures; i++)
		fl->feature[i]->trail = 0;

	return count;
}

/*
      void Klt::ExtractMovingPoints( jblas::mat &Xp, int fima, int mod, int numFeat)
      {
	    int i, j=0, l, aux, nima, dim;
	    dim = Xp.size2();
	    Xp.resize(numFeat, dim, 0);

	    // Extract the feature points in a jblas matrix
	    nima = fima+mod;

	    for (l=fima; l< nima; l++){
	      if (l==ft->nFrames)
		{
		  aux = (fima+mod)-l;
		  nima = aux;
		  l = 0;
		}
	      for ( i=0; i<ft->nFeatures; i++)
		{
		  if ( (ft->feature[i][l]->val == 0)  &&  ((fabs(ft->feature[i][l]->vx) > LIM ) || ( fabs(ft->feature[i][l]->vy) > LIM)))
		    {

		      Xp(j,0) = ft->feature[i][l]->x;
		      Xp(j,1) = ft->feature[i][l]->y;
		      Xp(j,2) = ft->feature[i][l]->vx;
		      Xp(j,3) = ft->feature[i][l]->vy;
		      Xp(j,4) = 0;
		      Xp(j,5) = ft->feature[i][l]->auxval;
		      //   std::cout << " j " << j ;
		      j++;
		      //	      ft->feature[i][l]->val = 1; //Bandera para que utilise el lugar
		    }
		  //	  else
		  //   ft->feature[i][l]->val = 1; //Bandera para que utilise el lugar
		}
	      //      std::cout << "numero de Imagen = " << l << "f= "<< fima << "j = "<< j << std::endl;
	    }

	  }
 */

void Klt::TrackFeatures(KLT_TrackingContext tc,
		cv::Mat &jImg1, cv::Mat &jImg2,
		KLT_FeatureList fl, bool replace,uchar *featmap, uchar *maskProb, double timeStamp )
{
	//  tDate dt = 0.1; //1.0
	//	double dt = timeStamp;
	//	    std::cout << " dt " << dt << "timeStamp " << timeStamp  << std::endl;
	//	klt::CompatTreatment treatment1, treatment2;
	KLT_PixelType *img1 = NULL, *img2 = NULL;
	int i, nLostFeatures=0, x,y;

	//	if (!tc->sequentialMode)
	img1 = jImg1.data;
	img2 = jImg2.data;
	KLTTrackFeatures(tc, img1, img2, jImg2.cols, jImg2.rows, fl);

	if (replace)
	{
		//  nLostFeatures = fl->nFeatures - KLTCountRemainingFeatures(fl);
		for (i = 0 ; i < fl->nFeatures ; i++){
			x   = (int) fl->feature[i]->errorx;
			y   = (int) fl->feature[i]->errory;
			// FIRST: Invert probability for all the points in its t-1 position (tracked or not)
			_invertFeaturemap(x, y, featmap, tc->mindist, jImg2.cols, jImg2.rows);
			if (fl->feature[i]->val >= 0)  {
				// SECOND: Update probability for all the tracked points in its current position
				x   = (int) fl->feature[i]->x;
				y   = (int) fl->feature[i]->y;
				_fillFeaturemap(x, y, featmap, maskProb, tc->mindist, jImg2.cols, jImg2.rows );
			}
			else nLostFeatures++;
		}

		if ( 1)  {
			fprintf(stderr,  "(KLT) Attempting to replace %d features "
					"in a %d by %d image...  \n", nLostFeatures, jImg2.cols, jImg2.rows);
			fflush(stderr);
		}

		/* If there are any lost features, replace them */
		if (nLostFeatures > 0) {

			KLTReplaceLostFeaturesWindow(tc, img2,jImg2.cols, jImg2.rows, fl,
					0, 0,
					jImg2.cols, jImg2.rows, featmap, maskProb);

		}
	}

	curFrame++;
	frame5++;

	//  calcVelocity(fl, ft, curFrame, dt);
	//  if (curFrame >= ft->nFrames) curFrame = 0;
	//	    std::cout << " ESTOY TRACKEANDO\n " << std::endl;
	//  StoreFeatureList(fl, ft, curFrame);
	//   StoreFeatureList(fl, ft5, frame5);
	//	if (!tc->sequentialMode)
	//	destroyCompatImage(img1, treatment1);
	//    free(index);
}


void Klt::featuresVelocity (KLT_FeatureList fl, double dt, int flag, bool _3d)
{
	calcVelocity(fl, ft, curFrame, dt, flag, _3d);
	if (curFrame >= ft->nFrames) curFrame = 0;
}


/****************************************************************************************
		#########################################################################################
		write shapes to img
		#########################################################################################
 ****************************************************************************************/

// draw a rect frame
void WriteRectToImg(cv::Mat &jRGBimg, int x1, int y1, int x2, int y2, int R, int G, int B)
{
	int ncols = jRGBimg.cols;
	int nrows = jRGBimg.rows;
	int width = jRGBimg.step/3;
	int offset;
	unsigned char* RGBimg = (unsigned char *)jRGBimg.data;
	int xx, yy;

	xx = x1;
	if (xx >= 0 && xx < ncols)
		for (yy = y1; yy <= y2; yy++)
			if (yy >= 0 && yy < nrows)
			{
				offset = (yy * width + xx)*3;
				RGBimg[offset+0] = B;
				RGBimg[offset+1] = G;
				RGBimg[offset+2] = R;
			}

	xx = x2;
	if (xx >= 0 && xx < ncols)
		for (yy = y1; yy <= y2; yy++)
			if (yy >= 0 && yy < nrows)
			{
				offset = (yy * width + xx)*3;
				RGBimg[offset+0] = B;
				RGBimg[offset+1] = G;
				RGBimg[offset+2] = R;
			}

	yy = y1;
	if (yy >= 0 && yy < nrows)
		for (xx = x1+1; xx < x2; xx++)
			if (xx >= 0 && xx < ncols)
			{
				offset = (yy * width + xx)*3;
				RGBimg[offset+0] = B;
				RGBimg[offset+1] = G;
				RGBimg[offset+2] = R;
			}

	yy = y2;
	if (yy >= 0 && yy < nrows)
		for (xx = x1+1; xx < x2; xx++)
			if (xx >= 0 && xx < ncols)
			{
				offset = (yy * width + xx)*3;
				RGBimg[offset+0] = B;
				RGBimg[offset+1] = G;
				RGBimg[offset+2] = R;
			}
}

// draw a point
void WritePointToImg(cv::Mat &jRGBimg, int x, int y, int R, int G, int B)
{
	int ncols = jRGBimg.cols;
	int nrows = jRGBimg.rows;
	int width = jRGBimg.step/3;
	int offset;
	unsigned char* RGBimg = (unsigned char *)jRGBimg.data;

	for (int yy = y - 1 ; yy <= y + 1 ; yy++)
		for (int xx = x - 1 ; xx <= x + 1 ; xx++)
			if (xx >= 0 && yy >= 0 && xx < ncols && yy < nrows)
			{
				offset = (yy * width + xx)*3;
				RGBimg[offset+0] = B;
				RGBimg[offset+1] = G;
				RGBimg[offset+2] = R;
			}
}

// draw a dashed point
void WriteHalfPointToImg(cv::Mat &jRGBimg, int x, int y, int R, int G, int B)
{
	int ncols = jRGBimg.cols;
	int nrows = jRGBimg.rows;
	int width = jRGBimg.step/3;
	int offset;
	unsigned char* RGBimg = (unsigned char *)jRGBimg.data;

	for (int yy = y - 1 ; yy <= y + 1 ; yy++)
		for (int xx = x - 1 ; xx <= x + 1 ; xx++)
			if (xx >= 0 && yy >= 0 && xx < ncols && yy < nrows)
			{
				if (yy==y || xx==x) continue;
				offset = (yy * width + xx)*3;
				RGBimg[offset+0] = B;
				RGBimg[offset+1] = G;
				RGBimg[offset+2] = R;
			}
}


// Bresenham's algorithm for drawing lines, adapted from
// http://www.cs.unc.edu/~mcmillan/comp136/Lecture6/Lines.html
void WriteLineToImg(cv::Mat &jRGBimg, int x0, int y0, int x1, int y1, int R, int G, int B)
{
	//	int ncols = jRGBimg.cols;
	int nrows = jRGBimg.rows;
	int width = jRGBimg.step/3;
	int size = 3*nrows*width;
	int offset;
	unsigned char* RGBimg = (unsigned char *)jRGBimg.data;

	int dy = y1 - y0;
	int dx = x1 - x0;
	int stepx, stepy;

	if (dy < 0) { dy = -dy;	stepy = -width; } else { stepy = width; }
	if (dx < 0) { dx = -dx;	stepx = -1; } else { stepx = 1; }
	dy <<= 1;
	dx <<= 1;

	y0 *= width;
	y1 *= width;
	offset = 3*(x0+y0);
	if (offset>=0 && offset<size)
	{
		RGBimg[offset+0] = B;
		RGBimg[offset+1] = G;
		RGBimg[offset+2] = R;
	}
	if (dx > dy)
	{
		int fraction = dy - (dx >> 1);
		while (x0 != x1)
		{
			if (fraction >= 0)
			{
				y0 += stepy;
				fraction -= dx;
			}
			x0 += stepx;
			fraction += dy;
			offset = 3*(x0+y0); if (offset<0 || offset>size) continue;
			if (offset>=0 && offset<size)
			{
				RGBimg[offset+0] = B;
				RGBimg[offset+1] = G;
				RGBimg[offset+2] = R;
			}
		}
	}
	else
	{
		int fraction = dx - (dy >> 1);
		while (y0 != y1)
		{
			if (fraction >= 0)
			{
				x0 += stepx;
				fraction -= dy;
			}
			y0 += stepy;
			fraction += dx;
			offset = 3*(x0+y0); if (offset<0 || offset>size) continue;
			if (offset>=0 && offset<size)
			{
				RGBimg[offset+0] = B;
				RGBimg[offset+1] = G;
				RGBimg[offset+2] = R;
			}
		}
	}
}

/****************************************************************************************
		#########################################################################################
		WriteImageToImg
		#########################################################################################
 ****************************************************************************************/
void WriteImageToImg(cv::Mat &jGreyimg, cv::Mat &jRGBimg)
{

	//	PRECOND(jRGBimg.colorSpace()==JfrImage_CS_BGR || jRGBimg.colorSpace()==JfrImage_CS_RGB, "jRGBimg is not in the RGB color space !");
	//	PRECOND(jRGBimg.depth()==IPL_DEPTH_8U, "jRGBimg has not u8 depth !");
	PRECOND(jRGBimg.channels()==3, "jRGBimg has not 3 channels !");
	PRECOND(jRGBimg.rows==jGreyimg.rows && jRGBimg.cols==jGreyimg.cols, "jRGBimg and jGreyimage must have same dimensions !");


	int ncols = jGreyimg.cols;
	int nrows = jGreyimg.rows;
	unsigned char* RGBimg = (unsigned char *)jRGBimg.data;
	unsigned char* greyimg = (unsigned char *)jGreyimg.data;

	int step = jGreyimg.step - jGreyimg.cols;
	//  std::cout<< "Step:  " << step << std::cout;

	//Verify if jGreyimg is really a grey image
	//	if (jGreyimg.colorSpace()==JfrImage_CS_BGR || jGreyimg.colorSpace()==JfrImage_CS_RGB)
	if (jGreyimg.channels()==3)
	{
		step=0;
		if (sizeof(KLT_PixelType) != 1)
			KLTWarning("(KLTWriteFeaturesToPPM)	KLT_PixelType is not uchar");
		//  RGBimg= greyimg;
		for (int x=0, offset=0 ; x<ncols; x++, offset+=3*step)
			for (int y=0; y<nrows; y++, offset+=3)
			{
				RGBimg[offset+0] = greyimg[offset+0];
				RGBimg[offset+1] = greyimg[offset+1];
				RGBimg[offset+2] = greyimg[offset+2];
			}
	}
	else
	{
		//  int step = jGreyimg.step() - jGreyimg.width();
		// Copy grey image to color image
		if (sizeof(KLT_PixelType) != 1)
			KLTWarning("(KLTWriteFeaturesToPPM)	KLT_PixelType is not uchar");
		for (int x=0, offset=0, i=0; x<ncols; x++, offset+=3*step, i+=step)
			for (int y=0; y<nrows; y++, offset+=3, i++)
			{
				RGBimg[offset+0] = greyimg[i];
				RGBimg[offset+1] = greyimg[i];
				RGBimg[offset+2] = greyimg[i];
			}
	}
}

#if 0
void Klt::WriteFeatureListToImg(KLT_FeatureTable ft, KLT_FeatureList fl, image::Image &jRGBimg, int im)
{
	JFR_PRECOND(jRGBimg.colorSpace()==JfrImage_CS_BGR || jRGBimg.colorSpace()==JfrImage_CS_RGB, "jRGBimg is not in the RGB color space !");
	JFR_PRECOND(jRGBimg.depth()==IPL_DEPTH_8U, "jRGBimg has not u8 depth !");
	JFR_PRECOND(jRGBimg.channels()==3, "jRGBimg has not 3 channels !");

	int initIm = 100;

	if (im==initIm)
	{
		// Overlay features in red
		for (int i = 0 ; i < fl->nFeatures ; i++) {//for per feature

			//Save 5 points in differents colors
			switch (i)
			{
			case 0: //RED COLOR
				WritePointToImg(jRGBimg, (int)(fl->feature[i]->x + 0.5), (int)(fl->feature[i]->y + 0.5), 255,0,0);
				break;
			case 1: // YELLOW COLOR
				WritePointToImg(jRGBimg, (int)(fl->feature[i]->x + 0.5), (int)(fl->feature[i]->y + 0.5), 255,255,0);
				break;
			case 2: // PURPLE COLOR
				WritePointToImg(jRGBimg, (int)(fl->feature[i]->x + 0.5), (int)(fl->feature[i]->y + 0.5), 150,0,255);
				break;
			case 3: //ORANGE COLOR
				WritePointToImg(jRGBimg, (int)(fl->feature[i]->x + 0.5), (int)(fl->feature[i]->y + 0.5), 255,150,0);
				break;
			case 4: //CYAN COLOR
				WritePointToImg(jRGBimg, (int)(fl->feature[i]->x + 0.5), (int)(fl->feature[i]->y + 0.5), 0,255,255);
				std::cout << "Imagen inicial " << im << std::endl;
				break;
			default:
				//points in white color
				WritePointToImg(jRGBimg, (int)(fl->feature[i]->x + 0.5), (int)(fl->feature[i]->y + 0.5), 255,255,255);

			}//end switch
		}//end for per features

	}//end if first frame
	else {
		for (int ima=initIm; ima<=im; ima++) {// for per frame
			for (int i = 0 ; i < fl->nFeatures ; i++){ //for per feature
				//Save 5 points in differents colors
				switch (i)
				{
				case 0: //RED COLOR
					if ((ft->feature[0][ima-initIm]->val == 0)||(ima==initIm))
						WritePointToImg(jRGBimg, (int)(ft->feature[0][ima-initIm]->x + 0.5), (int)(ft->feature[0][ima-initIm]->y + 0.5), 255,0,0);
					else //replace or out of bounds
						WriteHalfPointToImg(jRGBimg, (int)(ft->feature[0][ima-initIm]->x + 0.5), (int)(ft->feature[0][ima-initIm]->y + 0.5), 255,0,0);
					break;
				case 1: // YELLOW COLOR
					if ((ft->feature[1][ima-initIm]->val == 0)||(ima==initIm))
						WritePointToImg(jRGBimg, (int)(ft->feature[1][ima-initIm]->x + 0.5), (int)(ft->feature[1][ima-initIm]->y + 0.5), 255,255,0);
					else //replace or out of bounds
						WriteHalfPointToImg(jRGBimg, (int)(ft->feature[1][ima-initIm]->x + 0.5), (int)(ft->feature[1][ima-initIm]->y + 0.5), 255,255,0);
					break;
				case 2: // PURPLE COLOR
					if ((ft->feature[2][ima-initIm]->val == 0)||(ima==initIm))
						WritePointToImg(jRGBimg, (int)(ft->feature[2][ima-initIm]->x + 0.5), (int)(ft->feature[2][ima-initIm]->y + 0.5), 150,0,255);
					else //replace or out of bounds
						WriteHalfPointToImg(jRGBimg, (int)(ft->feature[2][ima-initIm]->x + 0.5), (int)(ft->feature[2][ima-initIm]->y + 0.5), 150,0,255);
					break;
				case 3: //ORANGE COLOR
					if ((ft->feature[3][ima-initIm]->val == 0)||(ima==initIm))
						WritePointToImg(jRGBimg, (int)(ft->feature[3][ima-initIm]->x + 0.5), (int)(ft->feature[3][ima-initIm]->y + 0.5), 255,150,0);
					else //replace or out of bounds
						WriteHalfPointToImg(jRGBimg, (int)(ft->feature[3][ima-initIm]->x + 0.5), (int)(ft->feature[3][ima-initIm]->y + 0.5), 255,150,0);
					break;
				case 4: //CYAN COLOR
					if ((ft->feature[4][ima-initIm]->val == 0)||(ima==initIm))
						WritePointToImg(jRGBimg, (int)(ft->feature[4][ima-initIm]->x + 0.5), (int)(ft->feature[4][ima-initIm]->y + 0.5), 0,255,255);
					else //replace or out of bounds
						WriteHalfPointToImg(jRGBimg, (int)(ft->feature[4][ima-initIm]->x + 0.5), (int)(ft->feature[4][ima-initIm]->y + 0.5), 0,255,255);
					break;
				default:
					//points in white color
					WritePointToImg(jRGBimg, (int)(ft->feature[i][ima-initIm]->x + 0.5), (int)(ft->feature[i][ima-initIm]->y + 0.5), 255,255,255);
				}//end switch
			}//end for per features

		} //end for per frame
	}   // end else


}// end function

#endif


void Klt::writeFeatureTableToImg(KLT_FeatureTable ft, int size, cv::Mat &jRGBimg)
{
	int i,j;
	int maxFrames = min(size,ft->nFrames);
	for(j = 0; j<ft->nFeatures; j++)
		for (int i = 0; i < maxFrames; i++)
		{
			if (ft->feature[j][i]->_3Dlost == 1)
				WritePointToImg(jRGBimg, (int)(ft->feature[j][i]->x + 0.5), (int)(ft->feature[j][i]->y + 0.5), 0,0,255); //Blue
			else if (ft->feature[j][i]->val == 0)
				WritePointToImg(jRGBimg, (int)(ft->feature[j][i]->x + 0.5), (int)(ft->feature[j][i]->y + 0.5), 255,0,0); //Tracked Red
			else if (ft->feature[j][i]->val > 0)
			{
				if (ft->feature[j][i]->errorval == KLT_LANDSCAPE)
				{
					WritePointToImg(jRGBimg, (int)(ft->feature[j][i]->x + 0.5), (int)(ft->feature[j][i]->y + 0.5), 250,250,0); //
					WriteHalfPointToImg(jRGBimg, (int)(ft->feature[j][i]->errorx + 0.5), (int)(ft->feature[j][i]->errory + 0.5), 255,0,255);
				}
				else
					WritePointToImg(jRGBimg, (int)(ft->feature[j][i]->x + 0.5), (int)(ft->feature[j][i]->y + 0.5), 0,200,0); // New Green
			}

		}
}


/****************************************************************************************
		#########################################################################################
		WriteFeatureListToImg
		#########################################################################################
 ****************************************************************************************/
void Klt::WriteFeatureListToImg(KLT_FeatureList fl, cv::Mat &jRGBimg)
{
	//	PRECOND(jRGBimg.colorSpace()==JfrImage_CS_BGR || jRGBimg.colorSpace()==JfrImage_CS_RGB, "jRGBimg is not in the RGB color space !");
	//	PRECOND(jRGBimg.depth()==IPL_DEPTH_8U, "jRGBimg has not u8 depth !");
	//	PRECOND(jRGBimg.channels()==3, "jRGBimg has not 3 channels !");

	// Overlay features in red
	for (int i = 0 ; i < fl->nFeatures ; i++)
	{
		if (fl->feature[i]->_3Dlost == 1)
			WritePointToImg(jRGBimg, (int)(fl->feature[i]->x + 0.5), (int)(fl->feature[i]->y + 0.5), 0,0,255); //Blue
		else if (fl->feature[i]->val == 0)
			WritePointToImg(jRGBimg, (int)(fl->feature[i]->x + 0.5), (int)(fl->feature[i]->y + 0.5), 255,0,0); //Tracked Red
		else if (fl->feature[i]->val > 0)
		{
			if (fl->feature[i]->errorval == KLT_LANDSCAPE)
			{
				WritePointToImg(jRGBimg, (int)(fl->feature[i]->x + 0.5), (int)(fl->feature[i]->y + 0.5), 250,250,0); //
				WriteHalfPointToImg(jRGBimg, (int)(fl->feature[i]->errorx + 0.5), (int)(fl->feature[i]->errory + 0.5), 255,0,255);
			}
			else
				WritePointToImg(jRGBimg, (int)(fl->feature[i]->x + 0.5), (int)(fl->feature[i]->y + 0.5), 0,200,0); // New Green
		}

	}
}

/****************************************************************************************
		#########################################################################################
		WriteMascaraToImg Mascara de Probabilidad
		#########################################################################################
 ****************************************************************************************/
void Klt::WriteMascaraToImg(uchar *mask, cv::Mat &jRGBimg)
{
	int i, j, val;

	//	PRECOND(jRGBimg.colorSpace()==JfrImage_CS_BGR || jRGBimg.colorSpace()==JfrImage_CS_RGB, "jRGBimg is not in the RGB color space !");
	PRECOND(jRGBimg.depth()==IPL_DEPTH_8U, "jRGBimg has not u8 depth !");
	PRECOND(jRGBimg.channels()==3, "jRGBimg has not 3 channels !");

	for (i=0; i <jRGBimg.cols; i++ )
		for (j=0; j < jRGBimg.rows; j++){
			val = mask[j*jRGBimg.cols+i];
			WritePointToImg(jRGBimg, i, j, val,val,val);
		}

}
}// end namespaces
