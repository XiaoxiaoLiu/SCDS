/* addrand.c - add a small random amount to each point in an effort to
 * make no 4 co-circular. Results are NOT guarenteed!
 */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#include <stdlib.h>
#include "math.h"
#include "cti_globals.h"
#include "unsafe_perturb_extern.h"


/* add a small random amount to the position of THIS point
 * Ret: 1 ok
 *      0 the pt could not be safely shifted
 */
int addrand_pt(int pt,		/* pt to be shifted */
	       float maxper	/* max amount that x or y may change */
	       )
{
    float tx,ty;
    int tries_left = 1000;	/* # of attempts to shift pt before */
				/* bailing out */

    /* check each point before it's position is changed */
    do
    {
	tx = savept[pt].x + maxper * ((float)(rand() % 32675) / 32676.);
	ty = savept[pt].y + maxper * ((float)(rand() % 32675) / 32676.);
	tries_left--;
    }
    while (unsafe_perturb(pt, tx, ty) && tries_left);
    if (! tries_left)
	return 0;		/* err - don't overwrite pt */
    savept[pt].x = tx;
    savept[pt].y = ty;
    return 1;
}


/* add a small random amount to the position of EVERY point */
void addrand_all (float maxper)	/* max amount that x or y may change */
{
    int pt;

    for (pt=0; pt < numpts; pt++)
    {
	addrand_pt(pt, maxper);
    }
}

/* add a small random amount to the position of each point IN THIS SLICE
 * Ret: 1 ok
 *      0 a pt in the slice could not be safely shifted
 */
int addrand_slice(int slice,	/* slice whose points are to be shifted */
		  float maxper	/* max amount that x or y may change */
		  )
{
    int pt, firstpt, con;
    
    for (con=0; con < num_contours[slice]; con++)
    {
	pt = firstpt = begin_contour[slice][con];
	do
	{
	    if (! addrand_pt(pt, maxper))
		return 0;	/* err */
	    pt = nextpt[pt];
	} while (pt != firstpt); 
    }
    return 1;
}

void addval_slice(int slice, float xval, float yval)
{
    int pt, firstpt, con;
    
    for (con=0; con < num_contours[slice]; con++)
    {
	pt = firstpt = begin_contour[slice][con];
	do
	{
	    savept[pt].x = savept[pt].xo + xval;
	    savept[pt].y = savept[pt].yo + yval;
	    pt = nextpt[pt];
	} while (pt != firstpt); 
    }
}
