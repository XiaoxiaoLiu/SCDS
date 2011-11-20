/* ctiler.h - prototypes for ctiler.c */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

extern int  ctiler( cti_flags_t flags,	 /* user specified parameters */
		            float perturb_dist	 /* (cm) max perturbation distance, if */
				  );	                 /* flags & CTI_PERT_DIST */

extern int setContours(int numContours, int* pointCounts,float*** contourPoints);
extern std::vector<std::vector<float> > getTiles();
extern void clearTiles();
