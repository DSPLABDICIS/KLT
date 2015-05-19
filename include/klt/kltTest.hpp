/*
 * kltTest.hpp
 *
 *  Created on: Feb 9, 2013
 *      Author: luzdora
 */

#ifndef KLTTEST_HPP_
#define KLTTEST_HPP_

#include "klt.hpp"

#include <math.h>
#include <cstring>
//#include "image/Image.hpp"

//#include "camera/cameraPinhole.hpp"
//#include "geom/t3d.hpp"
//#include "geom/t3dPointTools.hpp"

//#include "geom/t3dEuler.hpp"
//namespace klt {
     void initTrackContext( KLT_TrackingContext tc,
			    int affineConsistencyCheck=-1,
			    int window_size=11,
			    int nPyramidLevels=3,
			    int subsampling=2,
			    int mindist=10,
			    double max_residue=24,
			    int max_iterations=10,
			    double min_determinant=0.01);


    int countMovingTrackingFeatures( KLT_FeatureList fl);
    int countMovingFeaturesTotal(KLT_FeatureList fl, int modu);
    void KLTwriteImageToImg (cv::Mat &imGrey, cv::Mat & imColor);
    uchar* reserveOccupationMask (int ncols, int nrows);
    //   void initialiseOccupationMask (uchar & featuremap, int cols, int rows, int flag);
    void initialiseOccupationMask (uchar * featuremap, uchar *maskP, int cols, int rows);
    void initialiseProbaMask (uchar * maskProb, int sizeM, int mdist);
    void freeOccupationMask(uchar *featuremap);
    void setOccupationZoneInMask (uchar * mask, uchar *maskP, KLT_FeatureList fl, int ncols, int nrows, int mindist, int flag);
//    void featuresToWorldCoordinate (KLT_FeatureList fl, jafar::camera::CameraPinhole  sensor, jafar::geom::T3DEuler  cameraToOrigin, int imaH);
//    void featuresToNextImage (KLT_FeatureList fl, jafar::camera::CameraPinhole & sensor, jafar::geom::T3DEuler cameraToOrigin);
    void projectAlSuelo( boost_incl::vec3 & pointW,  boost_incl::vec3 & pWCam, boost_incl::vec3 & pS );

//  }


#endif /* KLTTEST_HPP_ */
