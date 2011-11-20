/* cross.c - determines if and where 2 lines segments intercept. There
 * are 2 flavors, one which considers the endpoints of the segments as
 * part of the segments (cross_exact) and one which excludes the
 * endpts (cross_exclusive).  (gst 5-Aug-91)
 */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#include "cti_gen.h"		/* for PRINTF */
#include "cross.h"
#include <stdexcept>

/* * * * * * helper routines * * * * * */

/* check that c is between a & b (or equal) */
int between_exact(double a,
		  double b,
		  double c
		  )
{
#define EPSILON (.0001)
  if ( (a-EPSILON) <= c && c <= (b+EPSILON)) return TRUE;
  if ( (b-EPSILON) <= c && c <= (a+EPSILON)) return TRUE;
  return FALSE;
}

/* check that c is between a & b (but not equal) */
int between_exclusive(double a,
		      double b,
		      double c
		      )
{
  if ( a < c && c < b) return TRUE;
  if ( b < c && c < a) return TRUE;
  return FALSE;
}


int In_Between(float iCoord, /* intersect coordinate */
	       float aCoord, /* coord from 1st node  */
	       float bCoord  /* coord from 2cd node  */
	       )
{
    if ((aCoord-bCoord) < 0)
	return((aCoord <= iCoord) && (iCoord <= bCoord));
    else
	return((bCoord <= iCoord) && (iCoord <= aCoord));
}


	    /* * * * * * * * CLASSIFY BEFORE CALC * * * * * * */

/* This is based on identifying the differences in the slopes of the
 * line segs before calculating their intersection, to detect special
 * cases BEFORE the intersection calc is done, which can simplify the
 * formulae a bit. */

/* With absolutely NO assumptions made about the connectivity of the 2
 * line segments, the manner in which the 2 line segments touch is
 * classified.  Vector-based math, from Graphics Gems: code from which
 * this is derived starts on p 700, theory on p 127 (Appendix
 * concerning intersection calc) */


int Classify_Intersection(float from1x, float from1y, /* endpts of the 2 lines */
				   float to1x,   float to1y,
				   float from2x, float from2y,
				   float to2x,   float to2y
				   )
{
    float a,b,c;
    /* These 3 dot products can be explained in many ways; here's a few:
     * (1) move point from1 onto point from2; then the projection of
     *     to1 onto the perpendicular of to2 will tell which vector points
     *     more clockwise:
     *       positive: 1
     *       negative: 2
     *       0: 1 & 2 have the same slope (neither is more clockwise)
     * (2) which side of vector 1 that point from2 is on
     * (3) which side of vector 2 that point from1 is on
     */

    a = (to1y - from1y)*(to2x - from2x) - (to1x - from1x)*(to2y - from2y);
    b = (from1x - from2x)*(to2y - from2y) - (from1y - from2y)*(to2x - from2x);
    c = (from1x - from2x)*(to1y - from1y) - (from1y - from2y)*(to1x - from1x);

    if (a == 0)			/* parallel? */
    {
	if (b == 0)		/* colinear? */
	{
	    if (In_Between(from1x, from2x,to2x) && In_Between(from1y, from2y,to2y) ||
		In_Between(  to1x, from2x,to2x) && In_Between(  to1y, from2y,to2y) ||
		In_Between(from2x, from1x,to1x) && In_Between(from2y, from1y,to1y) ||
		In_Between(  to2x, from1x,to1x) && In_Between(  to2y, from1y,to1y))
				/* overlapping? */
	    {
		if (from1x == from2x && from1y == from2y
		    && to1x-from1x>0 != to2x-from2x>0 ||
		    from1x ==   to2x && from1y ==   to2y
		    && to1x-from1x>0 == to2x-from2x>0 ||
		      to1x == from2x &&   to1y == from2y
		    && to1x-from1x>0 == to2x-from2x>0 ||
		      to1x ==   to2x &&   to1y ==   to2y
		    && to1x-from1x>0 != to2x-from2x>0)
		    return CROSS; /* one point cross */
		else
		    return COINCIDE; /* multi-point cross */
	    }
	    else
		return COLINEAR; /* don't cross */
	}
	else
	    return PARALLEL;	/* don't cross & are parallel */
    }

    else if (a > 0)		/* vector 1's slope is more clockwise? */
    {
	if (b >= 0 && b <= a &&
	    c >= 0 && c <= a)
	    return CROSS;
	else
	    return NO_CROSS;
    }

    else
    {
	if (b <= 0 && b >= a &&
	    c <= 0 && c >= a)
	    return CROSS;
	else
	    return NO_CROSS;
    }
}	


/* calc's where 2 line segments intercept, if they cross.  The endpts
 * of the lines are NOT valid touching pts, that is, the line segments are
 * "closed."
 *
 * Ret: 1=cross was found.
 *      0=cross was not found.
 *      crossing pt is returned in (xi_p,yi_p)
 */

int cross(float x1, float y1, float x2, float y2, /* endpts of first line */
	  float x3, float y3, float x4, float y4, /* endpts of second line */
	  float *xi_p, float *yi_p, /* RET: intercept */
	  int exact_flag	/* 1: seg endpts included, else excluded */
	  )
{
    double a1,b1,c1,		/* coefficients of line seg 1 */
           a2,b2,c2,		/* coefficients of line seg 2 */
           xi,yi;		/* intersection */
    

    /* inefficient -- only calc these if they are going to be used */
    a1 = y1-y2;
    a2 = y3-y4;
    b1 = x2-x1;
    b2 = x4-x3;
    c1 = (double)x1*(double)y2-(double)x2*(double)y1;
    c2 = (double)x3*(double)y4-(double)x4*(double)y3;

    /* Calc the intersection pt of the 2 LINES, then test if the
     * intersection pt is on both LINE SEGMENTS. */

    if (a1 == 0)
    {
	/* line 1 is horiz */
	*yi_p = yi = - c1/b1;
	if (a2 == 0)
	{
	    /* line 2 is horiz, too! test for line overlap */
	    *xi_p = xi = x1;		/* any x on the line seg is valid */
	    if (yi == (-c2 / b2))
	    {
		if (exact_flag)
		{
		    if (   between_exact((double)y1,(double)y2,yi)
			&& between_exact((double)y3,(double)y4,yi))
			return TRUE;
		    else
			return FALSE;
		}
		else
		{
		    if (   between_exclusive((double)y1,(double)y2,yi)
			&& between_exclusive((double)y3,(double)y4,yi))
			return TRUE;
		    else
			return FALSE;
		}
	    }
	}
	else
	    *xi_p = xi = - (b2*yi + c2)/a2; /* intersection of line 2 */
	                                    /* with a horiz line 1 */
    }
    
    if (b2*a1 == a2*b1)		/* does this need an EPSILON? */
    {
	/* lines are parallel! */

	if (c1 == c2)
	{
	    /* parallel lines are at same height */
	    if (exact_flag)
	    {
		if (   between_exact((double)x1,(double)x2,(double)x3)
		    || between_exact((double)x1,(double)x2,(double)x4)
		    || between_exact((double)x3,(double)x4,(double)x1)
		    || between_exact((double)x3,(double)x4,(double)x2) )
		    return TRUE;
		else
		    return FALSE;
	    }
	    else
	    {
		if (   between_exclusive((double)x1,(double)x2,(double)x3)
		    || between_exclusive((double)x1,(double)x2,(double)x4)
		    || between_exclusive((double)x3,(double)x4,(double)x1)
		    || between_exclusive((double)x3,(double)x4,(double)x2) )
		    return TRUE;
		else
		    return FALSE;
	    }
	}
	return FALSE;
    }
    
    if (a1 != 0)
    {
	/* line seg 1 is not horiz -- use standard intercept formula */
	*yi_p = yi = (a2*c1 - c2*a1) / (b2*a1 - a2*b1);
	*xi_p = xi = - (b1*yi + c1) / a1;
    }
    
    /* test if intersection pt is on both line segs */
    if (exact_flag)
    {
	if (   between_exact((double)x1,(double)x2,xi)
	    && between_exact((double)x3,(double)x4,xi)
	    && between_exact((double)y1,(double)y2,yi)
	    && between_exact((double)y3,(double)y4,yi) )
	    return TRUE;
	else
	    return FALSE;
    }
    else
    {
	if (   between_exclusive((double)x1,(double)x2,xi)
	    && between_exclusive((double)x3,(double)x4,xi)
	    && between_exclusive((double)y1,(double)y2,yi)
	    && between_exclusive((double)y3,(double)y4,yi) )
	    return TRUE;
	else
	    return FALSE;
    }
}

int cross_exact(float x1, float y1, float x2, float y2, /* endpts of first line */
		float x3, float y3, float x4, float y4, /* endpts of second line */
		float *xi_p, float *yi_p /* RET: intercept */
		)
{
  switch (Classify_Intersection(x1,y1, x2,y2, x3,y3, x4,y4))
    {
    case NO_CROSS:
    case PARALLEL:
    case COLINEAR:
      return 0;
    case ENDPTCROSS:
    case CROSS:
    case COINCIDE:
      return cross(x1,y1, x2,y2, x3,y3, x4,y4, xi_p,yi_p, 1);
    default:
      throw std::runtime_error("unknown case");
    }
}

int cross_exclusive(float x1, float y1, float x2, float y2, /* endpts of first line */
		    float x3, float y3, float x4, float y4, /* endpts of second line */
		    float *xi_p, float *yi_p /* RET: intercept */
		    )
{
    float xi,yi;
    if (cross(x1,y1, x2,y2, x3,y3, x4,y4, xi_p,yi_p, 1))
    {

	/* even if cross finds an intersection, return fail if the
	 * intersection pt is too close (identical) to either endpt of
	 * either line seg */

	xi = *xi_p;
	yi = *yi_p;
	if (   (FUZZY_ZERO(x1-xi, .001) && FUZZY_ZERO(y1-yi, .001))
	    || (FUZZY_ZERO(x2-xi, .001) && FUZZY_ZERO(y2-yi, .001))
	    || (FUZZY_ZERO(x3-xi, .001) && FUZZY_ZERO(y3-yi, .001))
	    || (FUZZY_ZERO(x4-xi, .001) && FUZZY_ZERO(y4-yi, .001)) )
	    return FALSE;
	else
	    return TRUE;
    }
    return FALSE;
}
