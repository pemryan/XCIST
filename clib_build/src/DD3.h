// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

#ifndef __DD3__
#define __DD3__

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))


// GE Proprietary

void DD3FlatWBackRow(float imgX,
		float imgXstep,
		int nrcols,
		float imgZstart,
		float imgZstep,
		int nrplanes,
		float *pImg,
		float *detX,
		int increment,
		float *detZ,
		float *scales,
		float z0,
		float *view,
		 int nrdetrows);
void DD3FlatWBackView(float x0,
		 float y0,
		 float z0,
		 int nrdetcols,
		 int nrdetrows,
		 int vertical,
		 float *xdi,
		 float *ydi,
		 float *detX, // empty distance array (nrxdist + 2)
		 float *detZ,
		 float *scales, // empty array (nrdetcols + 2)
		 float dzdx,
		 float *sinogram,
		 float *view, // empty view with 1 pixel margin all around
		 int nrcols,
		 int nrrows,       // images
		 int nrplanes,     //     do NOT
		 float *pOrig,     //        contain
		 float *pTrans);    //             a dummy 1 pixel frame
void DD3FlatWBack(float x0,
	     float y0,
	     float z0,
	     int nrdetcols,
	     int nrdetrows,
	     float *xds,
	     float *yds,
	     float *zds,
	     float dzdx,
	     float imgXoffset,
	     float imgYoffset,
	     float imgZoffset,
	     float *viewangles,
	     float *zshifts,
	     int nrviews,
	     float *sinogram,
	     int nrcols,         // image
	     int nrrows,         //    does NOT
	     int nrplanes,       //        contain a dummy 1 pixel frame
	      float *pOrig);
void DD3CurvedWBackRow(float imgX,
		float imgXstep,
		int nrcols,
		float imgZstart,
		float imgZstep,
		int nrplanes,
		float *pImg,
		float *detX,
		int increment,
		float *detZ,
		float *scales,
		float z0,
		float *view,
		 int nrdetrows);
void DD3CurvedWBackView(float x0,
		 float y0,
		 float z0,
		 int nrdetcols,
		 int nrdetrows,
		 int vertical,
		 float *xdi,
		 float *ydi,
		 float *detX, // empty distance array (nrxdist + 2)
		 float *detZ,
		 float *scales, // empty array (nrdetcols + 2)
		 float dzdx,
		 float *sinogram,
		 float *view, // empty view with 1 pixel margin all around
		 int nrcols,
		 int nrrows,       // images
		 int nrplanes,     //     do NOT
		 float *pOrig,     //        contain
		 float *pTrans);    //             a dummy 1 pixel frame
void DD3CurvedWBack(float x0,
	     float y0,
	     float z0,
	     int nrdetcols,
	     int nrdetrows,
	     float *xds,
	     float *yds,
	     float *zds,
	     float dzdx,
	     float imgXoffset,
	     float imgYoffset,
	     float imgZoffset,
	     float *viewangles,
	     float *zshifts,
	     int nrviews,
	     float *sinogram,
	     int nrcols,         // image
	     int nrrows,         //    does NOT
	     int nrplanes,       //        contain a dummy 1 pixel frame
	      float *pOrig);


/*
 * DD3 transpose
 */
void DD3Transpose(int nrcols,
		  int nrrows,
		  int nrplanes,
		  float *pOrig,
		  float *pTrans);




/*
 * DD3 boundaries
 */
void DD3Boundaries(int nrBoundaries,
		   float *pCenters,
		   float *pBoundaries);

/*
 * DD3 add transpose
 */
void DD3AddTranspose(int nrcols,
		     int nrrows,
		     int nrplanes,
		     float *pOrig,
		     float *pTrans);

/*
 * DD3 row projector
 */
void DD3ProjRow(float imgX,
		float imgXstep,
		int nrcols,
		float imgZstart,
		float imgZstep,
		int nrplanes,
		float *pImg,
		float *detX,
		int increment,
		float *detZ,
		float *scales,
		float z0,
		float *view,
		int nrdetrows);

void DD3ProjRow_mm(float imgX,
				   float imgXstep,
				   int nrcols,
				   float imgZstart,
				   float imgZstep,
				   int nrplanes,
				   float *pImg,
				   float *detX,
				   int increment,
				   float *detZ,
				   float *scales,
				   float z0,
				   float *view,
				   int nrdetrows,
				   int nrdetcols);



/*
 * DD3 row backprojector
 */
void DD3BackRow(float imgX,
		float imgXstep,
		int nrcols,
		float imgZstart,
		float imgZstep,
		int nrplanes,
		float *pImg,
		float *detX,
		int increment,
		float *detZ,
		float *scales,
		float z0,
		float *view,
		int nrdetrows);


void DD3BackRow_mm(float imgX,
				float imgXstep,
				int nrcols,
				float imgZstart,
				float imgZstep,
				int nrplanes,
				float *pImg,
				float *detX,
				int increment,
				float *detZ,
				float *scales,
				float z0,
				float *view,
				int nrdetrows,
				int nrdetcols);


/*
 * DD3 view projector
 */
void DD3ProjView(float x0,
		 float y0,
		 float z0,
		 int nrdetcols,
		 int nrdetrows,
		 int vertical,
		 float *xdi,
		 float *ydi,
		 float *detX,
		 float *detZ,
		 float *scales,
		 float dzdx,
		 float *view,
		 float *newproj,
		 int nrcols,
		 int nrrows,       // images
		 int nrplanes,     //     do NOT
		 float *pOrig,     //        contain
		 float *pTrans);   //             a dummy 1 pixel frame



/*
* DD3 view projector Kai's version with mm metric
*/
void DD3ProjView_mm(float x0,
				 float y0,
				 float z0,
				 int nrdetcols,
				 int nrdetrows,
				 int vertical,
				 float *xdi,
				 float *ydi,
				 float *detX,
				 float *detZ,
				 float *scales,
				 //float dzdx,
				 float *view,
				 float *newproj,
				 int nrcols,
				 int nrrows,       // images
				 int nrplanes,     //     do NOT
				 float *pOrig,     //        contain
				 float *pTrans,		//             a dummy 1 pixel frame
				 float vox_xy_size,
				 float vox_z_size);   //added fields





/*
 * DD3 view backprojector
 */
void DD3BackView(float x0,
		 float y0,
		 float z0,
		 int nrdetcols,
		 int nrdetrows,
		 int vertical,
		 float *xdi,
		 float *ydi,
		 float *detX, // empty distance array (nrxdist + 2)
		 float *detZ,
		 float *scales,		
		 float dzdx,
		 float *sinogram,
		 float *view, // empty view with 1 pixel margin all around
		 int nrcols,
		 int nrrows,       // images
		 int nrplanes,     //     do NOT
		 float *pOrig,     //        contain
		 float *pTrans);   //             a dummy 1 pixel frame

void DD3BackView_mm(float x0,
				 float y0,
				 float z0,
				 int nrdetcols,
				 int nrdetrows,
				 int vertical,
				 float *xdi,
				 float *ydi,
				 float *detX, // empty distance array (nrxdist + 2)
				 float *detZ,
				 float *scales,		
				 //float dzdx,
				 float *sinogram,
				 float *view, // empty view with 1 pixel margin all around
				 int nrcols,
				 int nrrows,       // images
				 int nrplanes,     //     do NOT
				 float *pOrig,     //        contain
				 float *pTrans,
				 float vox_xy_size,
				 float vox_z_size);   //             a dummy 1 pixel frame

/*
 * DD3 projector
 */
#ifdef __cplusplus
extern "C"{
#endif
	void DD3Proj(float x0,
	     float y0,
	     float z0,
	     int nrdetcols,
	     int nrdetrows,
	     float *xds,
	     float *yds,
	     float *zds,
	     float dzdx,
	     float imgXoffset,
	     float imgYoffset,
	     float imgZoffset,
	     float *viewangles,
	     float *zshifts,
	     int nrviews,
	     float *sinogram,
	     int nrcols,           // image
	     int nrrows,           //    does NOT
	     int nrplanes,         //        contain a dummy 1 pixel frame
	     float *pOrig);


//Kai's version, which all the measurement are ALL based on mm, optimized for helical scan
void DD3Proj_mm(float x0,
			 float y0,
			 float z0,
			 int nrdetcols,
			 int nrdetrows,
			 float *xds,
			 float *yds,
			 float *zds,
			 //float dzdx,
			 float imgXoffset,
			 float imgYoffset,
			 float imgZoffset,
			 float *viewangles,
			 float *zshifts,
			 int nrviews,
			 float *sinogram,
			 int nrcols,           // image
			 int nrrows,           //    does NOT
			 int nrplanes,         //        contain a dummy 1 pixel frame
			 float *pOrig,
			 float vox_xy_size, //added fields, mm
			 float vox_z_size); //added fields  mm


void DD3Back_mm(float x0,
			 float y0,
			 float z0,
			 int nrdetcols,
			 int nrdetrows,
			 float *xds,
			 float *yds,
			 float *zds,
			 //float dzdx,
			 float imgXoffset,
			 float imgYoffset,
			 float imgZoffset,
			 float *viewangles,
			 float *zshifts,
			 int nrviews,
			 float *sinogram,
			 int nrcols,           // image
			 int nrrows,           //    does NOT
			 int nrplanes,         //        contain a dummy 1 pixel frame
			 float *pOrig,
			 float vox_xy_size,
			 float vox_z_size);


/*
 * DD3 backprojector
 */
void DD3Back(float x0,
	     float y0,
	     float z0,
	     int nrdetcols,
	     int nrdetrows,
	     float *xds,
	     float *yds,
	     float *zds,
	     float dzdx,
	     float imgXoffset,
	     float imgYoffset,
	     float imgZoffset,
	     float *viewangles,
	     float *zshifts,
	     int nrviews,
	     float *sinogram,
	     int nrcols,           // image
	     int nrrows,           //    does NOT
	     int nrplanes,         //        contain a dummy 1 pixel frame
	     float *pOrig);


#ifdef __cplusplus
}
#endif

#endif

