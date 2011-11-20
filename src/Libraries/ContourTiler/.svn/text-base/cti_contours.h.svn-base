/* cti_contours.h - types for CTI (Contour TIler) contours */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

typedef struct {
    float x,y,z;	/* perturbed coordinates of a point */
    float xo,yo;	/* original copies */

    int con;		/* to tell if 2 pts are on the same contour */

} CTI_POINT;


/* these are the defines for Phase 1 contours (massager) */

/* static array bounds */
#ifndef MAX_NNODES
#define MAS_NNODES    1500	/* max # of pts/con */
#endif
#ifndef MAS_NCONS
#define MAS_NCONS     500	/* max # of cons/volume */
#endif
/* NOTE: MAS_NCONS & MAS_NNODES should be dynamic */

#define MAS_NDIMS       3	/* data nodes (x,y,z) are 3 */
				/* dimensional; this is not adjustable */
				/* -- it's just here to give a */
				/* symbolic name */



/*
 * The defines for Phase 2 contours; since the info kept for each
 * point takes more memory in Phase 2 than in Phase 1, these numbers
 * are smaller than their Phase counterparts.
 */

#ifndef CTI_MAX_ENDS
#define	CTI_MAX_ENDS          500 /* max # of pts per contour */
#endif

#ifndef CTI_MAX_CONTOURS
#define CTI_MAX_CONTOURS       40 /* max # of contours per slice */
#endif

#ifndef CTI_MAX_SLICES	
#define CTI_MAX_SLICES	      200 /* max # of slices */
#endif

#ifndef CTI_MAX_PTS_PER_SLICE
#define CTI_MAX_PTS_PER_SLICE 1500 /* max # of pts on all contours on a slice */
#endif

#ifndef CTI_MAX_PTS
#define CTI_MAX_PTS         150000 /* realistic count of all pts */
				  /* (not really MAX_SLICES * MAX_CONTOURS * MAX_ENDS) */
#endif





