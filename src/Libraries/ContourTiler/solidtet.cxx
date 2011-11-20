/* solidtet.c - for debugging the tiles, this is the file you hack to
 * select special views of the volumes.  Runtime control via TETFLAG
 * tells which tets should have their tiles visible. (jc/gst
 * 21-oct-92)
 */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#include "cti_globals.h"
#include "surfacetotile_extern.h"
#include <stdlib.h>
#include <math.h>


/* instead of outputing the surface as the tiles, output the VOLUME
 * occupied by the tets. We cut the volume into NCUTS very
 * thin slices and output the tiles for each slice. This creates a
 * very large number of tiles which completely fills the space when
 * viewed from a direction oblique to the cutting planes. Set NCUTS to
 * 1 to output a single cut plane at the base of the volumes.
 */
#define NCUTS 1

int count_t1t2;
int count_t12;

void solidtet(FILE *fptile, CTI_TETRA *tetra_p, int tetflag);
extern void write_tile(FILE *fptile,float x1,float y1,float z1,float x2,float y2,float z2,float x3,float y3,float z3);
void seek_max_tile(CTI_TETRA* tetra,int not_included);

void make_tiles(FILE *fptile,
	   CTI_TETRA tetra[],
	   int ntetra,
	   int allow_horiz,
	   int collect_caps,
	   int draw_normals
	   )
{
    int tet,neigh;
    int tetflag = atoi(getenv("TETFLAG")); /* see code for values */

    if (tetflag == 999)
	surfacetotile(fptile,tetra,ntetra,allow_horiz, collect_caps,draw_normals);
    else if (tetflag == 99)
    {
	for (tet = 0; tet < ntetra; tet++)
	{
	    for (neigh = 0; neigh < 4; neigh++)
	    {
		if (tetra[tet].t[neigh] == -1) {
		    draw_tile(fptile, &tetra[tet], neigh, allow_horiz,
			      0, 0);
		}
	    }
	}
    }
    else
    {
	printf("\n");
	count_t1t2 = count_t12 = 0;
	for (tet=0; tet < ntetra; tet++)
	{
	    solidtet(fptile, &tetra[tet], tetflag);
	}
	printf("\n");
	printf("***Number of T1 and T2: %d\n",count_t1t2);
	printf("***Number of T12      : %d\n",count_t12);
    }
}		

void solidtet(FILE *fptile, CTI_TETRA *tetra_p, int tetflag)
{
    int p0,p1,p2,p3;
    int i;
    float x0,x1,x2,x3,y0,y1,y2,y3,z0,z1,z2,z3;
    float dx0,dx1,dx2,dy0,dy1,dy2,dzz;
    float x02,x03,x12,x13,y02,y03,y12,y13,zz;
    float dx02,dx03,dx12,dx13,dy02,dy03,dy12,dy13;
 
    /* all trimed */
    if (tetflag == 0 && !(tetra_p->type & CTI_ELIM)) return;

    /* solid T1 */
    if (tetflag == 1 && (!(tetra_p->type & CTI_T1) ||
			 (tetra_p->type & CTI_ELIM))) return;

    /* solid T2 */
    if (tetflag == 2 && (!(tetra_p->type & CTI_T2) ||
			 (tetra_p->type & CTI_ELIM))) return;

    /* solid T12*/
    if (tetflag == 12 && (!(tetra_p->type & CTI_T12) ||
			   (tetra_p->type & CTI_ELIM))) return;

    /* solid T1 and T2,  1+2 */
    if (tetflag == 3 && (!(tetra_p->type & (CTI_T1 | CTI_T2)) ||
			  (tetra_p->type & CTI_ELIM))) return;

    /* all solid, 1+2+12 */
    if (tetflag == 15 && (tetra_p->type & CTI_ELIM)) return;

    /* trimed T1, append 0 */			
    if (tetflag == 10 && !(tetra_p->type & CTI_T1)) return;

    /* trimed T2 */
    if (tetflag == 20 && !(tetra_p->type & CTI_T2)) return;

    /* trimed T12 */
    if (tetflag == 120 && !(tetra_p->type & CTI_T12)) return;

    /* trimed T1 or T2 */
    if (tetflag == 30 && !(tetra_p->type & (CTI_T1 | CTI_T2))) return;
							   
    /* all T1, 1+10 */
    if (tetflag == 11 && !(tetra_p->type & CTI_T1)) return;

    /* all T2, 2+20 */
    if (tetflag == 22 && !(tetra_p->type & CTI_T2)) return;

    /* all T12, 12+120 */
    if (tetflag == 132 && !(tetra_p->type & CTI_T12)) return;

    /* all T1 or T2, 3+30 */
    if (tetflag == 33 && !(tetra_p->type & (CTI_T1 | CTI_T2))) return;

    /* specify 100 or else other than 999 and above to create all tiles */

    p0 = tetra_p->p[0];
    p1 = tetra_p->p[1];
    p2 = tetra_p->p[2];
    p3 = tetra_p->p[3];

    x0 = savept[p0].xo;
    x1 = savept[p1].xo;
    x2 = savept[p2].xo;
    x3 = savept[p3].xo;

    y0 = savept[p0].yo;
    y1 = savept[p1].yo;
    y2 = savept[p2].yo;
    y3 = savept[p3].yo;

    z0 = savept[p0].z;
    z1 = savept[p1].z;
    z2 = savept[p2].z;
    z3 = savept[p3].z;

    /* write_tile(fptile,x0,y0,z0, x1,y1,z1, x2,y2,z2); */
    write_tile(fptile,x1,y1,z1, x2,y2,z2, x3,y3,z3);
    write_tile(fptile,x2,y2,z2, x3,y3,z3, x0,y0,z0);
    write_tile(fptile,x3,y3,z3, x0,y0,z0, x1,y1,z1);

    if (tetra_p->type & (CTI_T2 | CTI_T1))
    {
	count_t1t2++;

	dx0 = (x3-x0)/NCUTS;
	dx1 = (x3-x1)/NCUTS;
	dx2 = (x3-x2)/NCUTS;

	dy0 = (y3-y0)/NCUTS;
	dy1 = (y3-y1)/NCUTS;
	dy2 = (y3-y2)/NCUTS;

	dzz = (z3-z0)/NCUTS;
	zz  = z0;

	for (i = 0; i<NCUTS; i++)
	{
	    write_tile(fptile,x0,y0,zz, x1,y1,zz, x2,y2,zz);

	    x0 += dx0;
	    y0 += dy0;
	    x1 += dx1;
	    y1 += dy1;
	    x2 += dx2;
	    y2 += dy2;
	    zz += dzz;
	}
    }

    else if (tetra_p->type & CTI_T12)
    {
	count_t12++;

	printf("***T12(%2d): vertex %d,%d; %d,%d\n", count_t12, p0,p1,p2,p3);
	/**/
	write_tile(fptile,x0,y0,z0, x1,y1,z1, x2,y2,z2); /**/

	dx02 = (x2-x0)/NCUTS;
	dx03 = (x3-x0)/NCUTS;
	dx12 = (x2-x1)/NCUTS;
	dx13 = (x3-x1)/NCUTS;
	
	dy02 = (y2-y0)/NCUTS;
	dy03 = (y3-y0)/NCUTS;
	dy12 = (y2-y1)/NCUTS;
	dy13 = (y3-y1)/NCUTS;

	dzz = (z3-z0)/NCUTS;
	zz  = z0;

	x02 = x03 = x0;
	x12 = x13 = x1;
	y02 = y03 = y0;
	y12 = y13 = y1;

	for (i = 0; i<NCUTS-1; i++)
	{
	    x02 += dx02;
	    x03 += dx03;
	    x12 += dx12;
	    x13 += dx13;

	    y02 += dy02;
	    y03 += dy03;
	    y12 += dy12;
	    y13 += dy13;
	    
	    zz  += dzz;

	    write_tile(fptile,x02,y02,zz, x03,y03,zz, x13,y13,zz);
	    write_tile(fptile,x02,y02,zz, x13,y13,zz, x12,y12,zz);
	}
    }
}


#define fmax(a,b) ((a)>(b))?(a):(b)
void seek_max_tile(CTI_TETRA* tetra,int not_included) {
    int vert;
    int pt1,pt2,pt3;
    float x1,y1,x2,y2,x3,y3,dx1,dx2,dx3,dy1,dy2,dy3,v;
    static float maxv = 0.0;
    static float maxpt1 = 0;
    static float maxpt2 = 0;
    static float maxpt3 = 0;

    if (tetra == NULL)
    {
	printf("The max tile consists of pts %f %f %f\n",
	       maxpt1,maxpt2,maxpt3);
	return;
    }

    for (vert = 0; vert < 4; vert++)
    {
	if (vert != not_included)
	{
	    pt3 = pt2;
	    pt2 = pt1;
	    pt1 = tetra->p[vert]; 
	}
    }

    x1 = savept[pt1].x;
    y1 = savept[pt1].y;

    x2 = savept[pt2].x;
    y2 = savept[pt2].y;

    x3 = savept[pt3].x;
    y3 = savept[pt3].y;

    dx1 = fabs(x1-x2);
    dx2 = fabs(x2-x3);
    dx3 = fabs(x3-x1);

    dy1 = fabs(y1-y2);
    dy2 = fabs(y2-y3);
    dy3 = fabs(y3-y1);

    v = fmax(dx1, fmax(dx2, dx3)) + fmax(dy1, fmax(dy2, dy3));
    if (v>maxv)
    {
	maxv = v;
	maxpt1 = pt1;
	maxpt2 = pt2;
	maxpt3 = pt3;
    }
}
