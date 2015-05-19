/*
 * klp.hpp
 *
 *  Created on: Feb 9, 2013
 *      Author: luzdora
 */

#ifndef KLT_HPP_
#define KLT_HPP_

#include "klt.h"
#include "error.h"
#include <cv.h>

//#include "image/Image.hpp"
#include "imageConversion.hpp"
#include "../boost_incl.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

namespace klt {

using namespace boost_incl;
typedef double tDate;

void WriteImageToImg(cv::Mat &jGreyimg, cv::Mat &jRGBimg);
void WriteRectToImg(cv::Mat &jRGBimg, int x1, int y1, int x2, int y2, int R, int G, int B);
void WritePointToImg(cv::Mat &jRGBimg, int x, int y, int R, int G, int B);
void WriteHalfPointToImg(cv::Mat &jRGBimg, int x, int y, int R, int G, int B);
void WriteLineToImg(cv::Mat &jRGBimg, int x0, int y0, int x1, int y1, int R, int G, int B);


class Klt
{
protected :
	int workingPos[2];
	int prevWorkingPos[2];
	int workingDims[2];
	int prevWorkingDims[2];

	//int trackingZone[4]; // = workingZone - borders
	int searchingZone[4];

	int securityMargin;
	double lastTimeStamp;

public:

	KLT_FeatureTable ft; //! table with an historical of all features, limited to 96 frames
	//	(circular buffer so that you can let the algorithm running)

	KLT_FeatureTable ft5;
	int curFrame; //! the current frame number in the feature table (frame modulo 96)
	int frame5;
	Klt(int securityMargin=10):securityMargin(securityMargin),lastTimeStamp(0),ft(NULL),ft5(NULL)
	{
		curFrame = 0;
		frame5 = 0;

		// std::cout << "Crea clase Klt " << std::endl;
	}
	~Klt()
	{
		if (ft != NULL) KLTFreeFeatureTable(ft);
	}
	KLT_TrackingContext CreateTrackingContext(void)
	{ return KLTCreateTrackingContext(); }

	KLT_FeatureList CreateFeatureList(int nFeatures)
	{
		ft = CreateFeatureTable(96, nFeatures);
		ft5 = CreateFeatureTable(5, nFeatures);
		//   curFrame = 0;
		//  std::cout << "Antes de salir Create Feature List" << std::endl;
		return KLTCreateFeatureList(nFeatures);
	}

	KLT_FeatureList ResizeFeatureList(KLT_FeatureList fl, int nP)
	{
		int i;
		FreeFeatureList(fl);
		fl =  KLTCreateFeatureList(nP);
		for (i = 0; i<nP; i++){
			fl->feature[i]->val = KLT_NOT_FOUND;
			fl->feature[i]->_3Dlost = 0;
		}
		return fl;
	}

	KLT_FeatureHistory CreateFeatureHistory(int nFrames)
	{ return KLTCreateFeatureHistory(nFrames); }

	KLT_FeatureTable CreateFeatureTable(int nFrames, int nFeatures)
	{ return KLTCreateFeatureTable(nFrames, nFeatures); }

	// Free
	void FreeTrackingContext(KLT_TrackingContext tc)
	{ return KLTFreeTrackingContext(tc); }

	void FreeFeatureList(KLT_FeatureList fl)
	{ return KLTFreeFeatureList(fl); }

	void FreeFeatureHistory(KLT_FeatureHistory fh)
	{ return KLTFreeFeatureHistory(fh); }

	void FreeFeatureTable(KLT_FeatureTable ft)
	{ return KLTFreeFeatureTable(ft); }

	// Processing
	void SelectGoodFeatures(KLT_TrackingContext tc,
			cv::Mat &jImg, KLT_FeatureList fl,
			double timeStamp, int r, uchar *featmap, uchar *maskProb);


	void TrackFeatures(KLT_TrackingContext tc,
			cv::Mat &jImg1, cv::Mat &jImg2,
			KLT_FeatureList fl, bool replace,uchar *featmap, uchar *maskProb,
			double timeStamp);

	void featuresVelocity (KLT_FeatureList fl, double dt, int flag , bool _3d);


	int ExtractMovingPoints(boost_incl::mat &Xp, int fima, int mod, KLT_FeatureList fl, bool _3D);

	// Utilities
	int CountRemainingFeatures(KLT_FeatureList fl)
	{ return KLTCountRemainingFeatures(fl); }

	void PrintTrackingContext(KLT_TrackingContext tc)
	{ return KLTPrintTrackingContext(tc); }

	void ChangeTCPyramid(KLT_TrackingContext tc, int search_range)
	{ return KLTChangeTCPyramid(tc, search_range); }

	void UpdateTCBorder(KLT_TrackingContext tc)
	{ return KLTUpdateTCBorder(tc); }

	void StopSequentialMode(KLT_TrackingContext tc)
	{ return KLTStopSequentialMode(tc); }

	void SetVerbosity(int verbosity)
	{ return KLTSetVerbosity(verbosity); }

	float ComputeSmoothSigma(KLT_TrackingContext tc)
	{ return _KLTComputeSmoothSigma(tc); }

	// Storing/Extracting Features
	//	  void StoreFeatureList(KLT_FeatureList fl, KLT_FeatureTable ft, int frame)
	void StoreFeatureList(KLT_FeatureList fl)
	{ return KLTStoreFeatureList(fl, ft, curFrame); }

	void StoreFeatureList_Trail(KLT_FeatureList fl)
	{ return KLTStoreFeatureList(fl, ft5, frame5); }


	void ExtractFeatureList(KLT_FeatureList fl, KLT_FeatureTable ft, int frame)
	{ return KLTExtractFeatureList(fl, ft, frame); }

	void StoreFeatureHistory(KLT_FeatureHistory fh, KLT_FeatureTable ft, int feat)
	{ return KLTStoreFeatureHistory(fh, ft, feat); }

	void ExtractFeatureHistory(KLT_FeatureHistory fh, KLT_FeatureTable ft, int feat)
	{ return KLTExtractFeatureHistory(fh, ft, feat); }

	// Writing/Reading

	void WriteFeatureListToPPM(KLT_FeatureList fl, KLT_PixelType *greyimg,
			int ncols, int nrows, char *filename)
	{ return KLTWriteFeatureListToPPM(fl, greyimg, ncols, nrows, filename); }


	void writeFeatureTableToImg(KLT_FeatureTable ft, int size, cv::Mat &jrgbimg);
	void WriteFeatureListToImg(KLT_FeatureList fl, cv::Mat &jRGBimg);
	void WriteMascaraToImg(uchar *mask, cv::Mat &jRGBimg);

	void WriteFeatureList(KLT_FeatureList fl, char *filename, char *fmt, int _3D)
	{ return KLTWriteFeatureList(fl, filename, fmt, _3D); }

	void WriteFeatureHistory(KLT_FeatureHistory fh, char *filename, char *fmt, int _3D)
	{ return KLTWriteFeatureHistory(fh, filename, fmt,_3D); }

	void WriteFeatureTable(char *filename, const char *fmt, bool _3D)
	{ return KLTWriteFeatureTable(ft, filename, fmt, _3D); }

	void WriteFeatureTableReduced(char *filename, const char *fmt, bool _3D)
	{ return KLTWriteFeatureTableReduced(ft, filename, fmt, _3D); }


	void WriteFeatureTableForClusters(char *filename, char *fmt)
	{ return writeFeatureTableForClusters(ft, filename,fmt); }

	void WriteFeatureTableForSLAM(char *filename, char *fmt)
	{ return writeFeatureTableForSLAM(ft5, filename,fmt); }


	KLT_FeatureList ReadFeatureList(KLT_FeatureList fl, char *filename)
	{ return KLTReadFeatureList(fl, filename); }

	KLT_FeatureHistory ReadFeatureHistory(KLT_FeatureHistory fh, char *filename)
	{ return KLTReadFeatureHistory(fh, filename); }

	KLT_FeatureTable ReadFeatureTable(KLT_FeatureTable ft, char *filename)
	{ return KLTReadFeatureTable(ft, filename); }
}; //end Klt



}//end klt

#endif /* KLT_HPP_ */
