/* surfacetotile.c - write a list of tiles found to lie on the surface
 * of this set of tet's. [There is some debuggery going on in this
 * file...]
 *
 * 20-sep-93 gst: define DRAW_NORMALS: additionally generates pointy
 *                tiles representing corresponding tile normal
 *                directions; use this to visualize each tiles' normal
 *                direction. 
*/

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#include "cti_globals.h"
#include "math.h"

extern void write_tile(FILE *fptile,float,float,float,float,float,float,float,float,float);

/* Draw the face of tetra_p which does NOT include the point
 * not_included (3 faces of the tet will include this point and 1 face
 * will not include this point):
 * - allow or disallow horizontal faces
 * - orient the face so that it points towards the outside of the tet,
 *   using the convention that the normal of the triangle is the
 *   front-facing side
 * 
 * Ret: TRUE: tile passed all tests and was drawn.
 */
int draw_tile(FILE *fptile,
	      CTI_TETRA *tetra_p,
	      int not_included, /* [0..3] the vertex who is NOT in this tile */
	      int allow_horiz,
	      int flip_normal,
	      int draw_normals
	      )
{

    double x1,y1,z1,		/* tile's vertices & 4th pt */
           x2,y2,z2,
           x3,y3,z3,
           x4,y4,z4;
    int p0,p1,p2,p3, vert;	/* tile's pts & 4th pt */
    double A,B,C,D;		/* coeff's of Ax + By + Cz + D = 0 */
    double planeside;		/* which side of the ABCD plane is the */
				/* 4th pt on? >0: one side <0: other */
				/* side =0: on plane */
    double xa,ya,za,len;
    // unused var // double xd,yd,zd;
    int flip;

    /* collect this tile's vertices, in order */
    for (vert = 0; vert < 4; vert++)
    {
	if (vert != not_included)
	{
	    /* ...and the little one said, "move over! move over!" */
	    p2=p1;
	    p1=p0;
	    p0=tetra_p->p[vert];  /* tricky! */
	}
    }
    p3 = tetra_p->p[not_included];

    /* Don't draw faces which lie in the contour planes (This z test
     * needs to be replaced by a test of the slicenum of each pt, so
     * that tilted slices will work.  -gst) */

    FPRINTF(stderr,"Considering face {%d %d %d}\n",p0,p1,p2);
    z1 = savept[p0].z;
    z2 = savept[p1].z;
    z3 = savept[p2].z;

    /* skip horizontal tiles if they are prohibited */
    if ( !allow_horiz && (z1 == z2) && (z2 == z3))
    {
	FPRINTF(stderr, "INFO: tet face [%d %d %d] is horizontal. Ignoring.\n",
		p0,p1,p2);
	return 0; /* return as "NOT drawn" */
    }

    /*seek_max_tile(tetra_p,not_included);  */

    x1 = savept[p0].x;		/* 0 */
    y1 = savept[p0].y;
    x2 = savept[p1].x;		/* 1 */
    y2 = savept[p1].y;
    x3 = savept[p2].x;		/* 2 */
    y3 = savept[p2].y;
    x4 = savept[p3].x;		/* 3 */
    y4 = savept[p3].y;

    z4 = savept[p3].z;		/* and the last z */

    /* enforce the tile's orientation, so that it always faces towards
     * it's front. To do this, make sure that the value of the 4th
     * point in the tet makes the equation of the plane negative,
     * because the 4th pt is on the backside of the tile. Note that I
     * could have found out sided-ness with just a dot-product but the
     * quadratic coefficients are used for drawing unit normals, too.
     */

    A = y1*(z2-z3) + y2*(z3-z1) + y3*(z1-z2);
    B = z1*(x2-x3) + z2*(x3-x1) + z3*(x1-x2);
    C = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2);
    D = - x1*(y2*z3 - y3*z2)
	- x2*(y3*z1 - y1*z3)
	- x3*(y1*z2 - y2*z1);

    planeside = A*x4 + B*y4 + C*z4 + D;	/* compute from perturbed positions */

    /* output original positions, flipping normals if desired */

    /* tricky! If not flipping normals, look for positive planeside, */
    /* or if flipping normals, look for negative */
    if ((planeside > 0) != flip_normal)
    {
	    /* order: 2 1 0 */
        write_tile(fptile,
		    savept[p2].xo, savept[p2].yo, z3,
		    savept[p1].xo, savept[p1].yo, z2,
		    savept[p0].xo, savept[p0].yo, z1);
    } else {
	    /* order: 1 2 0 */
        write_tile(fptile,
		    savept[p1].xo, savept[p1].yo, z2,
		    savept[p2].xo, savept[p2].yo, z3,
		    savept[p0].xo, savept[p0].yo, z1);
    }


    if (draw_normals)
    {
	/* draw the tile's unit normal */
	xa = (x1 + x2 + x3) / 3;
	ya = (y1 + y2 + y3) / 3;
	za = (z1 + z2 + z3) / 3;
	len  = sqrt(A*A + B*B + C*C); /* normalize length */
	flip = ((planeside > 0) != flip_normal) ? -1 : 1;
    write_tile(fptile,
		xa,ya,za,	/* geometric average */
		(xa + x1)/2,	/* in-plane, near centroid */
		(ya + y1)/2,
		(za + z1)/2,
		xa + flip*A/len, /* tip of unit normal */
		ya + flip*B/len,
		za + flip*C/len);
    }

    return 1;  /* return as "drawn" */
}



void surfacetotile(FILE *fptile,
		   CTI_TETRA tetra[],
		   int ntetra,
		   int allow_horiz, /* 1: allow horizontal tiles on surf */
		   int collect_caps, /* 1: update cap status */
		   int draw_normals
		   )
{
    int surfacesize = 0, got_one;
    int tet, neigh_tet, neigh, p, n;

    for (tet = 0; tet < ntetra; tet++)
    {
	/* don't even consider sculpted tet's */
	if (tetra[tet].type & CTI_ELIM)
	{
	    FPRINTF(stderr,"ignoring tet #%d: OUT\n", tet);
	    continue;
	}

	got_one = 0;
	for (neigh = 0; neigh < 4; neigh++)
	{
	    /* test for neighbors: if a tet has no neighbor on a side,
	     * then that tet has (at least) 1 face on the surface &
	     * the cons involved are ineligable for capping.
	     */

	    neigh_tet = tetra[tet].t[neigh];
	    if ( (neigh_tet == -1) || (tetra[neigh_tet].type & CTI_ELIM) )
	    {
		if (draw_tile(fptile, &tetra[tet], neigh, allow_horiz,
			      0, draw_normals))
		{
		    surfacesize++;
		    got_one = 1;
		}
	    }

	}

	if (collect_caps && got_one)
	{
	    /* no cap for cons touching this tile */
	    for (n=0; n < 4; n++)
	    {
		p = tetra[tet].p[n];
		/*if (caps [savept[p].con])
		    printf("NO CAP: slice=%d con=%d\n", slicenum[p], savept[p].con);*/
		caps [savept[p].con] = 0;
	    }
	}
    }
    /*seek_max_tile(NULL,0);   */
    //fprintf(stderr,"number of tri's on surface is %d\n", surfacesize);
}



