/*
 * imageConversion.cpp
 *
 *  Created on: Feb 19, 2013
 *      Author: luzdora
 */

#include "imageConversion.hpp"

namespace klt{
cv::Mat * convertImage(cv::Mat & image)
{
	int imsize[] = {image.rows,image.cols};
	cv::Mat *gray_image, *tmp;

	gray_image = new cv::Mat(2,imsize, CV_8UC1);
	tmp = new cv::Mat(2,imsize, CV_8UC1);
	double vmax, vmin;
	double alpha = 1, beta = 0;
	//cv::cvtColor(image,*tmp,cv::COLOR_RGB2GRAY);

	if((int)image.step1(1)!=1)
		cv::cvtColor(image,*tmp,CV_RGB2GRAY);
	else
		*tmp = image;

	cv::minMaxLoc(*tmp,&vmin,&vmax,NULL,NULL);
	if ((vmax-vmin) != 0)
	{
		alpha = 255.0/(vmax-vmin);
		beta = -vmin*alpha;
	}
	cv::convertScaleAbs(*tmp,*gray_image,alpha,beta);
	delete tmp;

	return(gray_image);
}
}//end namespace klt
