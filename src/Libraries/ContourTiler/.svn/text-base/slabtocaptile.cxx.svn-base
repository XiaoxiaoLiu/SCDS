/* slabtocaptile.c - write a list of tiles in the contours found to
 * have no surface tiles eminating from them
 */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#include "cti_globals.h"
#include "surfacetotile_extern.h"
#include "cti_cap.h"

static void draw_con_tiles(FILE *fptile,
			   CTI_TETRA tetra[],
			   int ntetra,
			   int tettype,
			   int con,
			   int slice,
			   int flip_normals,
			   int draw_normals)
{
    int surfacesize = 0;
    int tet;
    int c1,c2,c3, p0,p1,p2, not_included;
    CTI_TETRA *tetra_p;

    for (tet=0; tet < ntetra; tet++)
    {
	/* if the tet is the right type (T1 or T2) and 3 of the tets'
	 * pts are in the spec'd con, draw the 3 pts
	 */
	if ((tetra[tet].type & tettype) == tettype)
	{	
	    tetra_p = &tetra[tet];
	    for (not_included = 0; not_included < 4; not_included++)
	    {
		switch (not_included)
		{
		case 0:
		    p0 = tetra_p->p[3];
		    p1 = tetra_p->p[2];
		    p2 = tetra_p->p[1];
		    break;
		case 1:
		    p0 = tetra_p->p[3];
		    p1 = tetra_p->p[2];
		    p2 = tetra_p->p[0];
		    break;
		case 2:
		    p0 = tetra_p->p[3];
		    p1 = tetra_p->p[1];
		    p2 = tetra_p->p[0];
		    break;
		case 3:
		    p0 = tetra_p->p[2];
		    p1 = tetra_p->p[1];
		    p2 = tetra_p->p[0];
		    break;
		}
		
		// int p3 = tetra_p->p[not_included];

		c1 = savept[p0].con;
		c2 = savept[p1].con;
		c3 = savept[p2].con;

		if ((con == c1) && (c1 == c2) && (c2 == c3))
		{
		    if (draw_tile(fptile, &tetra[tet], not_included,
				  1, flip_normals, draw_normals))
			surfacesize++;
		}
	    }
	}
    }
    //fprintf(stderr,"number of tri's on cap is %d\n", surfacesize);
}


void draw_caps(FILE *fptile,
	       int slice,
	       int all_cons,	/* 1: cap all cons in slice */
				/* 0: cap cons having (caps[con]==1) */
	       int tettype,	/* get cap from a tet that points in this dir */
	       CTI_TETRA tetra[],
	       int ntetra,
	       int draw_normals)
{
    int con_index, con;

    for (con_index=0; con_index < num_contours[slice]; con_index++)
    {
	con = savept[begin_contour[slice][con_index]].con;
	if (all_cons || caps[con])
	    draw_con_tiles(fptile, tetra,ntetra, tettype, con, slice,
			   !all_cons, draw_normals);
    }
}
