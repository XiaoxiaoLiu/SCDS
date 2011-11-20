/* tri.h - CTIler Delaunay Triangles */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

/* History:
    1-OCT-90 gst  added CTI_ to all capitalized names, to avoid
                  name conflicts.
*/
 
/* a single Del Tri */
typedef struct tri_t
{
    int p[3];			/* the 3 pts which form the triangle;
				 * the endpoints of edge i are the two
				 * pts which are NOT p[i],
				 * specifically, p[(i+1)%3] and
				 * p[(i+2)%3]. If you don't get that,
				 * try this:
				 *   edge 0 is formed by pts 1 & 2
				 *   edge 1 is formed by pts 2 & 0
				 *   edge 2 is formed by pts 0 & 1
				 */

    float xc,yc;		/* tri's circumcenter (the center of the
				 * circle which goes through all 3 pts)
				 */

    /* for edge i: npts[i] & t[i] apply */

    char npts[3];		/* for edge i: npts[i] determines if (xs,ys) &
				 * (xe,ye) are valid pts:
				 *   0: ray (start pt valid)
				 *   1: endpoint of some other lineseg (none known)
				 *   2: ray (start & end pts valid)
				 *   3: lineseg (start & end pts valid)
				 *   9: the v-edge is outside of contour
				 */

    int t[3];			/* index of tri which shares this Del
				 * edge, or -1 if there is no tri */
    
    char outside;		/* 1: the triangle is outside, 0: otherwise */

} CTI_TRI;



