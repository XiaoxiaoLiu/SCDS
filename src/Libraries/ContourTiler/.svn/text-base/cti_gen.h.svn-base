#ifndef CTI_GEN_H
#define CTI_GEN_H

/* cti_gen.h - Contour TIler GENeral usage:
 *   standard I/O, NULL
 *   memory operations
 *   string operations
 *   conveniences: NITEMS, MIN,MAX, FUZZY_ZERO, ABS
 *   debuggery,
 *   TRUE,FALSE
 *   DOT_PRODUCT
 *   OUT_THERE
 */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#include "stdio.h"
#ifndef NULL
#define NULL 0
#endif

#ifndef NITEMS			/* tells # of entries in an array */
#define NITEMS(array) (sizeof(array) / sizeof(array[0]))
#endif

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif
#ifndef ABS
#define ABS(a) ((a)<0 ? -(a) : (a))
#endif

#define  PRINTF if (cti_debug)  printf
#define FPRINTF if (cti_debug) fprintf
//static int cti_debug = 0;
extern int cti_debug;

/*
 * FUZZY_ZERO tests for a float being near enough to zero to count as being
 * equal to zero.
 */
#define FUZZY_ZERO(a,epsilon) ( (-epsilon < (a)) && ((a) < epsilon) )

#ifndef TRUE
#define TRUE	1
#endif
#ifndef FALSE
#define	FALSE	0
#endif

/* The dot product is the projection of one 2d vector onto another
 * vector; it results in a scalar which is:
 *   > 0   angle between vectors is acute
 *   = 0  vectors are perpendicular
 *   < 0   angle between vectors is obtuse
 */
#define DOT_PRODUCT(x1,x2,y1,y2) ((x1)*(x2) + (y1)*(y2))


#define OUT_THERE 5000		/* my concept of infinity */

/* OUT_THERE is a threshold which determines how far out (and thus how
 * long) the sides of the tet's will be allowed to be. It also
 * determines where infinity is, when the circum-center of a set of 3
 * co-linear points is to be found: the center of the circle is placed
 * at infinity, with a radius of the distance to the points.
 */


#define DEFAULT_PERTURB_DIST (.01) /* (cm) distance CTI can move a */
				   /* point and still consider it the */
				   /* same point */

#endif
