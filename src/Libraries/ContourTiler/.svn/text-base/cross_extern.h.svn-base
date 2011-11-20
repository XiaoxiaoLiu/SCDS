/* cross_extern.h - prototypes for cross */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

extern int cross_exact(float x1, float y1,
		       float x2, float y2, /* endpts of first line */
		       float x3, float y3,
		       float x4, float y4, /* endpts of second line */
		       float *xi_p, float *yi_p /* RET: intercept */
		       );

extern int cross_exclusive(float x1, float y1,
			   float x2, float y2, /* endpts of first line */
			   float x3, float y3,
			   float x4, float y4, /* endpts of second line */
			   float *xi_p, float *yi_p /* RET: intercept */
			   );

extern int cross(float x1, float y1,
		 float x2, float y2, /* endpts of first line */
		 float x3, float y3,
		 float x4, float y4, /* endpts of second line */
		 float *xi_p, float *yi_p, /* RET: intercept */
		 int exact_flag	/* 1: seg endpts included, else excluded */
		 );

extern int between_exact(double a,
			 double b,
			 double c
			 );
extern int between_exclusive(double a,
			     double b,
			     double c
			     );


int
    Classify_Intersection(float from1x, float from1y, /* endpts of the 2 lines */
			  float to1x,   float to1y,
			  float from2x, float from2y,
			  float to2x,   float to2y
			  );


extern int In_Between(float iCoord, /* intersect coordinate */
		      float aCoord, /* coord from 1st node  */
		      float bCoord  /* coord from 2cd node  */
		      );
