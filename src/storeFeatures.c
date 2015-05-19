/*
 * storeFeatures.c
 *
 *  Created on: Feb 26, 2013
 *      Author: luzdora
 */

/*********************************************************************
 * storeFeatures.c
 *
 *********************************************************************/

/* Our includes */
#include "error.h"
#include "klt.h"


/*********************************************************************
 *
 */

void KLTStoreFeatureList(
		KLT_FeatureList fl,
		KLT_FeatureTable ft,
		int frame)
{
	int feat, aux;

	if (frame < 0 || frame >= ft->nFrames)
		KLTError("(KLTStoreFeatures) Frame number %d is not between 0 and %d",
				frame, ft->nFrames - 1);

	if (fl->nFeatures != ft->nFeatures)
	{
		//   KLTError("(KLTStoreFeatures) FeatureList and FeatureTable must "
		//       "have the same number of features");
		for (feat = 0 ; feat < fl->nFeatures ; feat++)  {
			ft->feature[feat][frame]->x   = fl->feature[feat]->x;
			ft->feature[feat][frame]->y   = fl->feature[feat]->y;
			ft->feature[feat][frame]->xem   = fl->feature[feat]->xem;
			ft->feature[feat][frame]->yem   = fl->feature[feat]->yem;
			ft->feature[feat][frame]->zem   = fl->feature[feat]->zem;
			ft->feature[feat][frame]->vx  = fl->feature[feat]->vx;
			ft->feature[feat][frame]->vy  = fl->feature[feat]->vy;
			ft->feature[feat][frame]->vz = fl->feature[feat]->vz;
			ft->feature[feat][frame]->val = fl->feature[feat]->val;
			ft->feature[feat][frame]->_3Dlost = fl->feature[feat]->_3Dlost;
			ft->feature[feat][frame]->auxval = fl->feature[feat]->auxval;
			ft->feature[feat][frame]->trail = fl->feature[feat]->trail;
			ft->feature[feat][frame]->_3Dx = fl->feature[feat]->_3Dx;
			ft->feature[feat][frame]->_3Dy = fl->feature[feat]->_3Dy;
			ft->feature[feat][frame]->_3Dz = fl->feature[feat]->_3Dz;

		}
		for (aux = feat ; aux < ft->nFeatures; aux++)
		{
			ft->feature[aux][frame]->x   = 0.0;
			ft->feature[aux][frame]->y   = 0.0;
			ft->feature[feat][frame]->xem= 0.0;
			ft->feature[feat][frame]->yem= 0.0;
			ft->feature[feat][frame]->zem= 0.0;
			ft->feature[aux][frame]->vx  = 0.0;
			ft->feature[aux][frame]->vy  = 0.0;
			ft->feature[feat][frame]->vz = 0.0;
			ft->feature[aux][frame]->val = -9;
			ft->feature[feat][frame]->_3Dlost = 0;
			ft->feature[feat][frame]->_3Dx = 0.0;
			ft->feature[feat][frame]->_3Dy = 0.0;
			ft->feature[feat][frame]->_3Dz = 0.0;
		}
	}
	else
		for (feat = 0 ; feat < fl->nFeatures ; feat++)  {
			ft->feature[feat][frame]->x   = fl->feature[feat]->x;
			ft->feature[feat][frame]->y   = fl->feature[feat]->y;
			ft->feature[feat][frame]->xem   = fl->feature[feat]->xem;
			ft->feature[feat][frame]->yem   = fl->feature[feat]->yem;
			ft->feature[feat][frame]->zem   = fl->feature[feat]->zem;
			ft->feature[feat][frame]->vx  = fl->feature[feat]->vx;
			ft->feature[feat][frame]->vy  = fl->feature[feat]->vy;
			ft->feature[feat][frame]->vz = fl->feature[feat]->vz;
			ft->feature[feat][frame]->val = fl->feature[feat]->val;
			ft->feature[feat][frame]->_3Dlost = fl->feature[feat]->_3Dlost;
			ft->feature[feat][frame]->auxval = fl->feature[feat]->auxval;
			ft->feature[feat][frame]->trail = fl->feature[feat]->trail;
			ft->feature[feat][frame]->_3Dx = fl->feature[feat]->_3Dx;
			ft->feature[feat][frame]->_3Dy = fl->feature[feat]->_3Dy;
			ft->feature[feat][frame]->_3Dz = fl->feature[feat]->_3Dz;
		}
}


/*********************************************************************
 *
 */

void KLTExtractFeatureList(
		KLT_FeatureList fl,
		KLT_FeatureTable ft,
		int frame)
{
	int feat;

	if (frame < 0 || frame >= ft->nFrames)
		KLTError("(KLTExtractFeatures) Frame number %d is not between 0 and %d",
				frame, ft->nFrames - 1);

	/*  if (fl->nFeatures != ft->nFeatures)
    KLTError("(KLTExtractFeatures) FeatureList and FeatureTable must "
		"have the same number of features");

	 */


	for (feat = 0 ; feat < fl->nFeatures ; feat++)  {
		fl->feature[feat]->x   = ft->feature[feat][frame]->x;
		fl->feature[feat]->y   = ft->feature[feat][frame]->y;
		fl->feature[feat]->val = ft->feature[feat][frame]->val;
		fl->feature[feat]->_3Dlost = ft->feature[feat][frame]->_3Dlost;
		fl->feature[feat]->_3Dx    = ft->feature[feat][frame]->_3Dx;
		fl->feature[feat]->_3Dy    = ft->feature[feat][frame]->_3Dy;
		fl->feature[feat]->_3Dz    = ft->feature[feat][frame]->_3Dz;
	}
}


/*********************************************************************
 *
 */

void KLTStoreFeatureHistory(
		KLT_FeatureHistory fh,
		KLT_FeatureTable ft,
		int feat)
{
	int frame;

	if (feat < 0 || feat >= ft->nFeatures)
		KLTError("(KLTStoreFeatureHistory) Feature number %d is not between 0 and %d",
				feat, ft->nFeatures - 1);

	if (fh->nFrames != ft->nFrames)
		KLTError("(KLTStoreFeatureHistory) FeatureHistory and FeatureTable must "
				"have the same number of frames");

	for (frame = 0 ; frame < fh->nFrames ; frame++)  {
		ft->feature[feat][frame]->x   = fh->feature[frame]->x;
		ft->feature[feat][frame]->y   = fh->feature[frame]->y;
		ft->feature[feat][frame]->val = fh->feature[frame]->val;
		ft->feature[feat][frame]->_3Dlost 	= fh->feature[frame]->_3Dlost;
		ft->feature[feat][frame]->_3Dx 		= fh->feature[frame]->_3Dx;
		ft->feature[feat][frame]->_3Dy 		= fh->feature[frame]->_3Dy;
		ft->feature[feat][frame]->_3Dz 		= fh->feature[frame]->_3Dz;
	}
}


/*********************************************************************
 *
 */

void KLTExtractFeatureHistory(
		KLT_FeatureHistory fh,
		KLT_FeatureTable ft,
		int feat)
{
	int frame;

	if (feat < 0 || feat >= ft->nFeatures)
		KLTError("(KLTExtractFeatureHistory) Feature number %d is not between 0 and %d",
				feat, ft->nFeatures - 1);

	if (fh->nFrames != ft->nFrames)
		KLTError("(KLTExtractFeatureHistory) FeatureHistory and FeatureTable must "
				"have the same number of frames");

	for (frame = 0 ; frame < fh->nFrames ; frame++)  {
		fh->feature[frame]->x   = ft->feature[feat][frame]->x;
		fh->feature[frame]->y   = ft->feature[feat][frame]->y;
		fh->feature[frame]->val = ft->feature[feat][frame]->val;
		fh->feature[frame]->_3Dlost = ft->feature[feat][frame]->_3Dlost;
		fh->feature[frame]->_3Dx   	= ft->feature[feat][frame]->_3Dx;
		fh->feature[frame]->_3Dy 	= ft->feature[feat][frame]->_3Dy;
		fh->feature[frame]->_3Dz   	= ft->feature[feat][frame]->_3Dz;

	}
}



