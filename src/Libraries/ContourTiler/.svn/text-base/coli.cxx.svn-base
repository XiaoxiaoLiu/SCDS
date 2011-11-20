/* coli - tests if 3 points are (nearly) colinear */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#include "cti_gen.h"		/* for FUZZY_ZERO */
#define EPSILON (.0001)		/* how close is "close" ? */

int coli(float x1,float y1, float x2,float y2, float x3,float y3)
{
    double a;

    /* check if the x & y deltas from pt 1 to pt 2 are similiar to */
    /* those from pt 1 to pt 3 */
    if ( FUZZY_ZERO(x3-x1, EPSILON) ) {
	if ( FUZZY_ZERO(x3-x2, EPSILON) )
	     return 1;
        else return 0;
    }

    if ( FUZZY_ZERO(y3-y1, EPSILON) ) {
	if ( FUZZY_ZERO(y3-y2, EPSILON) )
	     return 1;
	else return 0;
    }

    a = (double)(x2-x1)/(double)(x3-x1) - (double)(y2-y1)/(double)(y3-y1);
    if ( FUZZY_ZERO(a, EPSILON) )
	 return 1;
    else return 0;

}

