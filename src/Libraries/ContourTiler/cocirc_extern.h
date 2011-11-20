/* cocirc_extern.h - prototypes for cocirc.c */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

extern int cocirc(float x1,float y1,
		  float x2,float y2,
		  float x3,float y3,
		  float x4,float y4);

extern int circenter_infinity(float x1,float y1,
			      float x2,float y2,
			      float x3,float y3,
			      float *xc_p,float *yc_p /* Ret: circle's center */
			      );

extern int circenter(float x1,float y1,
		     float x2,float y2,
		     float x3,float y3,
		     double *xc_p,double *yc_p /* Ret: circle's center */
		     );

extern void bisect(float x1,float y1,
		   float x2,float y2,
		   double *a_p,double *b_p,double *c_p	/* Ret: line Ax + By + C = 0 */
		   );

extern void linetween(float x1,float y1,
		      float x2,float y2,
		      float *a_p,float *b_p,float *c_p /* Ret: line Ax + By + C = 0 */
		      );

extern void crosspt(double a1,double b1,double c1, /* line 1 */
		    double a2,double b2,double c2, /* line 2 */
		    double *xi_p,double *yi_p	/* Ret: intersection pt */
		    );

