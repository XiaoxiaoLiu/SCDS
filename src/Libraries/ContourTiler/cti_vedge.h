/* cti_vedge.h - CTIler Voronoi edges */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

/* History:
    1-OCT-90 gst  added CTI_ to all capitalized names, to avoid
                  name conflicts.
*/

typedef struct vedge_t
{
    float xs,ys,xe,ye;
    int npts; /* 0,1, or 2 */
    int s,e;     /* pts dual to this edge */
    int ts,te;   /* pts whose bisectors intersect this bisector */
} CTI_VEDGE;
