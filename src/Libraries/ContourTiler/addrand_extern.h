/* addrand_extern.h - prototypes for addrand.c */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

extern int addrand_pt(int pt,	/* pt to be shifted */
		      float maxper /* max amount that x or y may change */
		      );

extern void addrand_all (float maxper); /* max amount that x or y may change */


extern int addrand_slice(int slice, /* slice whose points are to be shifted */
			 float maxper /* max amount that x or y may change */
			 );


extern void addval_slice(int slice, float xval, float yval);

