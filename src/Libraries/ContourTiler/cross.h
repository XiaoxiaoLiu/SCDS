/* cross.h - general purpose line crossing detection defines */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

enum CROSS_E			/* these are the ways in which 2 line segments
				 * might intersect each other */
{
    NO_CROSS,			/* don't cross & are skew */
    PARALLEL,			/* don't cross & slopes are equal */
    COLINEAR,			/* don't cross & slopes are equal & colinear */
    COINCIDE,			/* cross in many points & are colinear */
    CROSS,			/* cross at a single point */
    ENDPTCROSS 			/* one of the endpoints is shared */
};


/*
 * Much of the rest has been copied from cti_gen.h, because you might not have
 * that file...
 */

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/*
 * FUZZY_ZERO tests for a float being near enough to zero to count as being
 * equal to zero.
 */
#ifndef FUZZY_ZERO
#define FUZZY_ZERO(a,epsilon) ( (-epsilon < (a)) && ((a) < epsilon) )
#endif

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
