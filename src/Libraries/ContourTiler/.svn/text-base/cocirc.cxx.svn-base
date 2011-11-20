/* cocirc.c - determine if 4 points are cocircular, that is, a circle
 * can be drawn which (nearly) intersects all of them */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#include "cti_gen.h"		/* for FUZZY_ZERO */
#include "coli_extern.h"

#define EPSILON (.001)		/* how close is "cocircular close" ? */


/* find the line that the 2 pts intersect */
void linetween(float x1,float y1,
	       float x2,float y2,
	       float *a_p,float *b_p,float *c_p /* Ret: line Ax + By + C = 0 */
	       )
{
    *a_p = y1-y2;
    *b_p = x2-x1;
    *c_p = x1*y2 - x2*y1;
}


/* find the perpendicular bisector of the 2 pts */
void bisect(float x1,float y1,
	    float x2,float y2,
	    double *a_p,double *b_p,double *c_p	/* Ret: line Ax + By + C = 0 */
	    )
{
    double xx1=(double)x1;
    double yy1=(double)y1;
    double xx2=(double)x2;
    double yy2=(double)y2;

    *a_p = 2.0*(xx2-xx1);
    *b_p = 2.0*(yy2-yy1);
    *c_p = (xx1*xx1 - xx2*xx2 + yy1*yy1 - yy2*yy2);
}


/* find the pt at which 2 lines interesect;
 *  the lines are known to cross and are NOT parallel.
 */
void crosspt(double a1,double b1,double c1, /* line 1 */
	     double a2,double b2,double c2, /* line 2 */
	     double *xi_p,double *yi_p	/* Ret: intersection pt */
	     )
{
  double yi;

  PRINTF("a1=%f b1=%f c1=%f  a2=%f b2=%f c2=%f\n",a1,b1,c1,a2,b2,c2);/**/

  if (! FUZZY_ZERO(a1, EPSILON))
  {
      PRINTF("crosspt: a2==0\n");
      *yi_p = yi = (a2*c1 - c2*a1) / (b2*a1 - a2*b1);
      *xi_p = - (b1* yi + c1) / a1;
      return;
  }

  if (! FUZZY_ZERO(a2, EPSILON))
  {
      PRINTF("crosspt: a1==0\n");
      *yi_p = yi = (a1*c2 - c1*a2) / (b1*a2 - a1*b2);
      *xi_p = - (b2* yi + c2) / a2;
      return;
  }

  PRINTF("crosspt: a1==a2==0 SHOULDN'T HAPPEN!\n");
}



/* find the circumcenter of the 3 points, that is, return the center
 * of the circle in (*xc_p,*yc_p) that all 3 points intersect. This
 * version of the routine tests if the circle exists. Ret: TRUE if
 * circumcenter exists, else FALSE (which means that (*xc_p,*yc_p) is
 * invalid).
 */
int circenter(float x1,float y1,
	      float x2,float y2,
	      float x3,float y3,
	      double *xc_p,double *yc_p /* Ret: circle's center */
	      )
{
    double a1,b1,c1, a2,b2,c2;

    /*
     * calc the lines formed by:
     *   - the pts equidistance to pts 1 & 2
     *   - the pts equidistance to pts 2 & 3
     * calc the intersection of these 2 lines; this is
     *   the circle's center
     */

    if ( coli(x1,y1, x2,y2, x3,y3) )
    {
	return FALSE;  /* co-linear pts CANT be co-circular */
    }

    bisect(x1,y1, x2,y2, &a1,&b1,&c1);
    bisect(x3,y3, x2,y2, &a2,&b2,&c2);

    crosspt(a1,b1,c1, a2,b2,c2, xc_p, yc_p);

    return TRUE;
}


int cocirc(float x1,float y1,
	   float x2,float y2,
	   float x3,float y3,
	   float x4,float y4)
{
    double xc,yc;
    double rad, dist;

    /*
     * if the dist from pt 4 to the center of the circumcenter of
     * the other 3 is the same as the radius of the circumcenter,
     * then pt 4 is cocircular with the other 3.
     */

    if (! circenter(x1,y1, x2,y2, x3,y3, &xc,&yc))
    {
	PRINTF("cocirc:  circenter ret FALSE\n");
	return FALSE;  /* can't find circumcenter -- pts are co-linear */
    }

    PRINTF("cocirc: xc(%f) yc(%f)\n",xc,yc);

    /*rad =  sqrt((x1-xc)*(x1-xc) + (y1-yc)*(y1-yc));*/
    rad =  ((x1-xc)*(x1-xc) + (y1-yc)*(y1-yc));
    /*dist = sqrt((x4-xc)*(x4-xc) + (y4-yc)*(y4-yc));*/
    dist = ((x4-xc)*(x4-xc) + (y4-yc)*(y4-yc));
    if (FUZZY_ZERO(dist-rad, EPSILON))
	return TRUE;
    return FALSE;
}




