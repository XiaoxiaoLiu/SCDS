/* del.c - compute a Delaunay Triangulation */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

//#define DRAW/* DEBUG -- define to draw iterations*/
//#define DRAWFINAL/* DEBUG -- define to draw final del tri */

#include "math.h"
#include "cti_globals.h"
#include "cross_extern.h"
#include "insertpt_extern.h"
#include "coli_extern.h"
#include "cocirc_extern.h"
#include "unsafe_perturb_extern.h"

#include <string.h>
#include <stdlib.h>		/* for memcpy */

#define COCIR_EPS (.00001)	/* threshhold for co-circular */
#define ADD_EPS (.05)		/* minimum length of contour seg should be maintained */
//#define BREAK_LENGTH (2.0)	/* minimum length of contour seg should be broken out*/
#define MIN_OBTUSE_ANGLE -0.156 /* approximately equals to cos(0.55*M_PI) */
#define MIN_OBTUSE_TRI_RAD 0.0  /* minimum radius of obtuse triangle we deal with */


typedef struct cirtri		/* local definition of triangle used */
				/* for del() only */
{
    int p[3];			/* savept relative pt index */
    double xc,yc;
    double rad;
} CIRTRI;

typedef struct edge
{
    int s, e;			/* edges start & end indices */
    char deleted;		/* 1: duplicate edge  */
} EDGE;

typedef struct detour		/* when triaglation is not constrained, this is used */
{				/* to record the points we should go back to after   */
    int pt;			/* breaking (middle) points are processed. */
    struct detour *next;
} DETOUR;

extern void sortpts(int p[3]);
int find_obtuse(int slice,CIRTRI tri[],int ntri);
int obtuse(float x0,float y0,float x1,float y1,float x2,float y2,float &xp,float &yp);

/*

  Build up the Del Tri, by iteratively adding each pt:

  For each Del triangle (trn) which contains this pt INSIDE its
  circumcircle, remove trn from the tri list and add trns edges to the
  edge list. Duplicate edges are thrown away. The remaining edges form
  a star-shaped polygon, and are combined (thus consuming them) with
  the current pt one edge at a time to form triangles, and then merged
  into the current Del Tri. An initial Del Tri is formed by a
  rectangular bounding box of all of the points, called the "frame,"
  which is turned 45 degrees on alternate slices, or turned as a
  possible conflict resolution when something snags.

  Note that this method will NOT always find a Del Tri, because
  thresholds are used to detect when 2 points are colinear or
  incident, and there is no such mathematical restriction on the
  contour edges.  There is a limited precision available on each
  floating point processor, and it is feasible that this precision may
  be exceeded by the particular math that we have chosen.

  There are 2 special cases, handled quite elegantly by suggestion of
  Bernhard Geiger:

  (1) If the pt happens to be ON (or almost on) the circumcircle, this
  indicates that we must choose one of the many equally good Del Tris.
  We chose to perturb the point, revert both the tri list & the edge
  list to their states before any tris involving this pt were merged,
  and process the pt again. Perturbing the point may or may not move
  it off of the circumcircle, so repetition of this step may be
  necessary.  I check (in unsafe_perturb) that the perturbance does
  not either (a) cause any contour edges to cross or (b) change the
  geometry of the contour. LATER, perturbance in the normal direction
  of the tangent of the circumcircle may avoid this repetition, at
  least with this particular circumcircle.

  (2) The points are processed with respect to the contour edges, to
  assure that the triangulation is constrained. This is done by
  verifying that each edge that is deleted by addedges (an edge is
  deleted if it is added twice -- a misleading name, eh?) is not on
  the contour.  Note that we may have split the contour segments, so
  testing deleted edges against the original contour segments is not
  sufficient; we must test against the split contour segements.  If a
  *contour* edge (as opposed to a non-contour edge) has been marked
  for deletion, the triangulation could become unconstrained, so a
  point is added at the midpoint of that contour edge, spliting the
  edge into 2 colinear edges. The current points contribution to the
  triangulation is ignored, and the new mid-point is processed
  immediately (an advantage of an iterative implementation is that the
  processing order may be changed "on the fly"), after which the
  normal processing of the contour points is resumed with the point
  whose contribution was ignored. If the edge is still not sufficient
  to constrain the triangulation, the splitting process is repeated
  until either the triangulation is constrained or no more points may
  be added to the contour (points are only added if they are not
  incidental with other points, within a threshold).

  A slower alternative is checking that the edge (starting with the
  current point) is present in the edge list before merging the edge
  list into the current Tri. If it is not present, add one or more
  points to guarentee that the next iteration with this pt will
  include the current edge. This is slower because it performs the
  entire Del Tri to find out if it is constrained.

  Note that if, for some reason, the triangulation fails, then all
  inserted points are kept.
  
*/


/* create a Del Triangle in the local format */
void cre8cirtri(int p1,int p2,int p3, 
	   double xc,double yc, double rad,
	   CIRTRI tri[CTI_NTRI],
	   int *ntri_p)
{
  int ntri = *ntri_p;
  // unused var// int a;

  tri[ntri].p[0] = p1;
  tri[ntri].p[1] = p2;
  tri[ntri].p[2] = p3;
  tri[ntri].xc = xc;
  tri[ntri].yc = yc;
  tri[ntri].rad = rad;
  
  (*ntri_p)++;
}

/* remove a Del Triangle */
void rmcirtri(int trrm, CIRTRI tri[CTI_NTRI], int *ntri_p)
{
    /* slide the array from trrm to the end down 1 entry */
    memcpy((void*)(& tri[trrm]), (void*)(& tri[trrm+1]),
	  sizeof (CIRTRI) * (--(*ntri_p) - trrm));
}



/*
 * "XOR" an edge s..e to the spec'd edge list. Adding edge E1 which is
 * *already* in the list as E2 will delete E2 from the list and ignore
 * the request to add E1.  Ret: 1 or -1 if con edge was deleted (1 for
 * edge s..e, -1 for edge e..s), 0 for ok
 */
int addedge(int s, int e,	/* start & end point indices */
	    EDGE edge[], int *nedge_p /* the edge list & count */
	    )
{
    int e1;
    PRINTF("edges[%d]: ",*nedge_p);
    for (e1 = 0; e1 < *nedge_p; e1++)
    {
	if (   (edge[e1].s == s && edge[e1].e == e)
	    || (edge[e1].s == e && edge[e1].e == s) )
	{
	    edge[e1].deleted = 1;
	    PRINTF("DUP!\n");

	    if (s >= 0 && e >= 0)
	    {
		if (nextpt[s] == e)
		    return 1;	/* deleted a con edge */
		if (prevpt[s] == e)
		    return -1;	/* deleted a reverse con edge */
	    }
	    return 0;	/* duplicate */
	}
	else PRINTF("%d/%d ",edge[e1].s, edge[e1].e);
    }
    edge[*nedge_p].s = s;
    edge[*nedge_p].e = e;
    edge[*nedge_p].deleted = 0;
    (*nedge_p)++;
    PRINTF("edge added: %d/%d\n",s,e);
    return 0;
}


/* add the edges in the triangle to the edge list. See addedge.
 * Ret: err: index of pt whose con edge was deleted
 *      ok: -1
 */
int addedges(int tr, CIRTRI tri[CTI_NTRI], EDGE edge[], int *nedge_p)
{
    int pt, ret, s,e;
    for (pt = 0; pt < 3; pt++)
    {
	s = tri[tr].p[pt];
	e = tri[tr].p[(pt+1)%3];
	if ((ret = addedge(s, e, edge, nedge_p)) != 0)
	{
	    if (ret == 1)
		return s;
	    else if (ret == -1)
		return e;
	}
    }
    return -1;
}

/* select one of the 2 sets of pre-calc'd frame points, and make a
 * 2-triangle frame from them. */
void cre8frametris(CIRTRI tri[CTI_NTRI], int *ntri_p)
{
    double xc,yc, xd,yd, rad;	/* center, distance, radius */
    int fr = framestart;

    if (turn45)			/* if frame is tilted... */
	fr  += 4;		/* ...use the second frame */
    turn45 = 1 - turn45;	/* prepare for next slice */

    /* points 0 and 2 are assumed to be opposite each other */
    xd = (savept[fr+2].x-savept[fr].x)/2.0; /* distance */
    yd = (savept[fr+2].y-savept[fr].y)/2.0;
    rad = sqrt(xd*xd + yd*yd);
    xc = (savept[fr+2].x+savept[fr].x)/2.0; /* midpoint */
    yc = (savept[fr+2].y+savept[fr].y)/2.0;
    cre8cirtri(fr  ,
	       fr+1,
	       fr+2,
	       xc,yc, rad, tri,ntri_p);

    /* the frame is assumed square, so (rad, xc, yc) are the same */
    cre8cirtri(fr  ,
	       fr+2,
	       fr+3,
	       xc,yc, rad, tri,ntri_p);
}

/*
 * convert the format of a tri from the local format to the returned
 * format; do not allow a triangle to be converted twice
 */
void cre8tri(int s, int e, int t,
	double xc, double yc,
	CTI_TRI tri[CTI_NTRI],
	int *ntri_p
	)
{
    int tr, i, p[3], ntri = *ntri_p;

    p[0] = s;
    p[1] = e;
    p[2] = t;
    sortpts(p);

    /* this could use a binary search (both lists are ordered) */
    for (tr=0; tr < ntri; tr++)
    {
	if ((tri[tr].p[0] == p[0]) &&
	    (tri[tr].p[1] == p[1]) &&
	    (tri[tr].p[2] == p[2]))
	    return;
    }

    for (i=0; i<3; i++) {
	tri[ntri].p[i] = p[i];
	tri[ntri].npts[i] = 0;
	tri[ntri].t[i] = -1;
    }
    tri[ntri].xc = xc;
    tri[ntri].yc = yc;
    tri[ntri].outside = 0;
    (*ntri_p)++;
}


/*
 * convert ncirtri tri's, even those containing frame points, from the
 * local format to the returned format
 */
void cvtcir2tri(CIRTRI cirtri[CTI_NTRI],
	   int ncirtri,
	   CTI_TRI ctitri[CTI_NTRI],
	   int *nctitri_p
	   )
{
    int tr1;

    for (tr1 = 0; tr1 < ncirtri; tr1++)
	cre8tri(cirtri[tr1].p[0], cirtri[tr1].p[1], cirtri[tr1].p[2],
		cirtri[tr1].xc, cirtri[tr1].yc, ctitri, nctitri_p);
}




/* compute del tri of edges in a slice; overwrite tri[] with new triangles */
int del(int s,			/* slice to triangulate */
	CTI_TRI tri[CTI_NTRI],	/* RET: Del Tri */
	int *ntri_p,		/* RET: # of tris in tri[] */
	float maxper		/* maximum pertubation distance */
	)
{
    int pt, badpt;	/* pt index; -1 means "invalid" */
    // unused var//int detourpt;
    DETOUR *detour_hd, *detour_ele; /* record detoured points */
    double rad,dist, xc,yc, x1,y1, x2,y2, x3,y3; 
    // unused vars// double x4,y4;
    int con, startpt;		/* curr contour & first pt on con */
    int ntries;

    /* LATER: all of these arrays should be made dynamic */
    CIRTRI nowtri[CTI_NTRI];	/* current incremental triangulation */
    int   nnowtri;		/* # of tris in nowtri */
    CIRTRI lasttri[CTI_NTRI];	/* last incremental triangulation */
    int   nlasttri;		/* # of tri's in lasttri */

    EDGE edge[4*CTI_NTRI];	/* tri edges */
    int nedge;

    int trn, ed;		/* tri index, edge index */
    double xd,yd;		/* delta */
    float newx,newy;		/* location of a newly inserted pt */
    float tx,ty;		/* temp perturbed postion */

    int contour_edge_in;	/* if a contour edge is a side of a triangle */
    int num_obtuse_check=0;

 restart:

    nnowtri = nlasttri = 0;
    cre8frametris(nowtri,&nnowtri);

#ifdef BREAK_LENGTH
    /* break long contour segments into smaller ones so that double points may 
     * not appear.
     */
    break_long_segments(s);
#endif

    for(con=0; con<num_contours[s]; con++)
    {
	startpt = pt = begin_contour[s][con];
	detour_hd = NULL;
	do
	{
	    nedge = trn = 0;
	    x3 = savept[pt].x; y3 = savept[pt].y;

	    /* create tri save state in case pt is ON circumcenter or
	     * a middle point is added.
	     */
	    memcpy(lasttri, nowtri, sizeof (CIRTRI) * nnowtri);
	    nlasttri = nnowtri;

	    while (trn < nnowtri)
	    {
		xd = x3 - nowtri[trn].xc;
		yd = y3 - nowtri[trn].yc;
		dist = sqrt(xd*xd + yd*yd); /* LATER: take out sqrt() */
		
		/* ON circumcircle? */
		if (FUZZY_ZERO(dist - nowtri[trn].rad, COCIR_EPS)) 
		{
		    /* perturb and restart */
		    ntries = 0;
		    do
		    {
			tx = savept[pt].x + 
			     maxper * ((float)(rand() % 32675) / 32676.);
			ty = savept[pt].y + 
			     maxper * ((float)(rand() % 32675) / 32676.);
			ntries++;
		    }
		    while (unsafe_perturb(pt, tx, ty) && ntries < 1000);
		    if (ntries == 1000)
		    {
			printf(" ERR: PERTURB required, but couldn't: pt = %d\n",
			       pt);		    
			return 6;
		    }

		    PRINTF("  PERTURB: move pt%d from (%f,%f) to (%f,%f)\n",
			   pt,x3,y3,tx,ty);

		    savept[pt].x = tx;
		    savept[pt].y = ty;
		    x3 = savept[pt].x; y3 = savept[pt].y;

		    /* revert to save state */
		    nedge = trn = 0; /* ignore previous edges; reset index */
		    memcpy(nowtri, lasttri, sizeof (CIRTRI) * nlasttri);
		    nnowtri = nlasttri;
		    continue;
		}
		else if (dist < nowtri[trn].rad) /* INSIDE circumcircle? */
		{
		    badpt = addedges(trn, nowtri, edge, &nedge);
		    if (badpt == -1) /* Del Tri is still constrained... */
		    {
			rmcirtri(trn, nowtri, &nnowtri);
		    }
		    else	/* Del Tri is NO LONGER constrained... */
		    {
			newx = (savept[badpt].x + savept[nextpt[badpt]].x)/2;
			newy = (savept[badpt].y + savept[nextpt[badpt]].y)/2;
			PRINTF("  INSERT: midpt %d (%f,%f) after badpt %d\n",
			       numpts,newx,newy,badpt);
			insertpt(newx, newy, s,con,badpt);
			x3 = newx; y3 = newy;

			/* revert to save state */
			nedge = trn = 0; /* ignore previous edges; reset index */
			memcpy(nowtri, lasttri, sizeof (CIRTRI) * nlasttri);
			nnowtri = nlasttri;
		    
			detour_ele = (DETOUR *)(malloc(sizeof(DETOUR)));
			detour_ele->pt = pt;
			detour_ele->next = detour_hd;
			detour_hd = detour_ele;

			pt = nextpt[badpt];
			continue;
		    }
		}
		else		/* pt is outside circumcircle, so ignore tri */
		    trn++;
	    }

	    /* if the contour edge prevpt[pt]-pt is not added as a 
	     * triangle edge this time, it is very likely it will be
	     * missed after trianglation is done. Add a middle point
	     * before pt will solve this problem. But this may cause
	     * the current slab to abort because of too many points
	     * added if the contours of this slab are complicated. 
	     * Thus a threshhold is given here.
	     */
	    contour_edge_in=0;
	    for (ed=0; ed<nedge; ed++)
	    {
		if (edge[ed].s==prevpt[pt] || edge[ed].e==prevpt[pt])
		{
		    contour_edge_in = 1;
		    break;
		}
	    }

	    if ( !contour_edge_in && pt != startpt )
	    {
		if (fabs(savept[prevpt[pt]].x - x3) > ADD_EPS ||
		    fabs(savept[prevpt[pt]].y - y3) > ADD_EPS)
		{
		    newx = (savept[prevpt[pt]].x + x3)/2.0;
		    newy = (savept[prevpt[pt]].y + y3)/2.0;
		    
		    PRINTF("  ADD: midpt %d (%f,%f) before pt %d\n", 
			   numpts,newx,newy,pt);
		    insertpt(newx,newy,s,con,prevpt[pt]);

		    /* revert to save state */
		    nedge = trn = 0; /* ignore previous edges; reset index */
		    memcpy(nowtri, lasttri, sizeof (CIRTRI) * nlasttri);
		    nnowtri = nlasttri;
		    
		    detour_ele = (DETOUR *)(malloc(sizeof(DETOUR)));
		    detour_ele->pt = pt;
		    detour_ele->next = detour_hd;
		    detour_hd = detour_ele;

		    pt = prevpt[pt]; /* now pt is the new added point */
		    continue;
		}
		else
		{
		    printf("\n");
		    printf("  WARNING: Triangulation is not constrained!\n");
		    printf("           A contour edge is too short to break.\n");
		}
	    }

	    /* for each non-deleted edge ed in the list, create a tri
	     * containing ed and this pt */
	    for (ed=0; ed<nedge; ed++)
	    {
		if (edge[ed].deleted)
		    continue;

		x1 = savept[edge[ed].s].x;
		y1 = savept[edge[ed].s].y;
		x2 = savept[edge[ed].e].x;
		y2 = savept[edge[ed].e].y;

		if (! circenter(x1,y1, x2,y2, x3,y3, &xc,&yc))
		{
		    printf("  ERR: circenter doesn't exist!\n");
		    continue;
		}

		xd = x2-xc;
		yd = y2-yc;
		rad = sqrt(xd*xd + yd*yd);
		cre8cirtri(edge[ed].s,edge[ed].e,pt, xc,yc, rad, nowtri,&nnowtri);
	    }
		
	    if (detour_hd != NULL) /* continue to a detour point if any */
	    {
		DETOUR *temp;
		temp=detour_hd;
		pt = detour_hd->pt;
		detour_hd = detour_hd->next;
		PRINTF("  DETOUR: %d\n", pt);
		free(temp);
	    }
	    else
	    {
		pt = nextpt[pt]; /* success! continue to next point */
	    }
	} while (pt != startpt);
    }

#ifdef MIN_OBTUSE_ANGLE
    /* usually one time of eliminating obtuse triangles is enough. */
    if (num_obtuse_check<1 && find_obtuse(s,nowtri,nnowtri))
    {
	num_obtuse_check++;
	PRINTF("  Found obtuse triangles. Retry triangulating.\n");
	goto restart;
    }
#endif

    /* convert tri's from local format to returned format. */
    *ntri_p = 0;
    cvtcir2tri(nowtri, nnowtri, tri, ntri_p);

    return 0;		/* aok */
}



/* this routine suggests that moving the spec'd pt to location (x,y)
 * *may* corrupt the geometry of the contours on the slice. It checks
 * that the triangle swept out as the contour segments on which this
 * pt exists move from it's original position to the newly suggested
 * position does not intercept any other contour edges.  LATER:
 * Bernhard suggests that this could go faster if we precompute a list
 * (once per pt, based on it's original position) of con edges which
 * are within the max perturb distance from pt, and only check those
 * edges in this test.
 *
 * Ret: 0 if safe
 *      1,2,3: unsafe because of cross
 *      4,6: unsafe because perturbance is along con edge
 *      5: pt is a frame pt -- cannot be perturbed
*/

int unsafe_perturb(int pt,	/* point to be perturbed */
		   float x,	/* new position of pt */
		   float y
		   )
{
    int prev = prevpt[pt];
    int next = nextpt[pt];

    int slice, con;		/* of pt */
    int thispt, thisnext;	/* index of pt's on con edge */
    int startpt;		/* first pt on con */
    float x1,y1, x2,y2;		/* con edge to test */
    float xi,yi;		/* intersept */

    /* calc intersection of each of 2 triangles and all contour lines
     * in the slice (that's 2 * numpts intersections): one tri goes
     * from prevpt[pt] to old position of pt to new position of pt,
     * and the other uses a different first point: nextpt[pt] instead
     * of prevpt[pt]. No intercepts means the point can (probably)
     * safely be moved to the new position without crossing a contour
     * segment. Note that only 3 of the edges in the 2 triangles need
     * to be tested, (sorry -- I can't remember why...). 2 of the
     * contour edges, obviously, don't count as true intersections:
     * the 2 contained in the triangle edges. */

    /* touches a neighboring contour edge */
    if (coli(savept[pt].x, savept[pt].y,
	     savept[next].x, savept[next].y,
	     x, y))
	return 6;
    if (coli(savept[pt].x, savept[pt].y,
	     savept[prev].x, savept[prev].y,
	     x, y))
	return 4;

    if (-1 == (slice = slicenum[pt]))
	return 5;

    /* testing for failure */
    for(con=0; con < num_contours[slice]; con++)
    {
	thispt = startpt = begin_contour[slice][con];
	do
	{
	    thisnext = nextpt[thispt];
	    x1 = savept[thispt  ].x;
	    y1 = savept[thispt  ].y;
	    x2 = savept[thisnext].x;
	    y2 = savept[thisnext].y;

	    if (cross_exclusive (savept[pt].x, savept[pt].y, x, y,
				 x1, y1, x2, y2,
				 &xi, &yi))
	    {
		return 1;
	    }

	    if ( (thispt != prev) && (thispt != pt) )  /**/
	    {
		if (cross_exclusive (savept[next].x, savept[next].y, x,y,
				     x1, y1, x2, y2,
				     &xi, &yi))
		{
		    return 2;
		}
	    }

	    if ( (thispt != prev) && (thispt != pt) ) /**/
	    {
		if (cross_exclusive (savept[prev].x, savept[prev].y, x,y,
				     x1, y1, x2, y2,
				     &xi, &yi))
		{
		    return 3;
		}
	    }
	    
	    thispt = thisnext;
	} while (thispt != startpt);
    }

    return 0;
}

#ifdef BREAK_LENGTH
/* simple solution to solve double point problem. doesn't work in all cases and 
 * it usually adds too many points thus slows down the processing */
break_long_segments(slice)
int slice;
{
    int con, pt, startpt, ppt;
    float dx,dy,dist,newx,newy;
    float maxdist = 0.0;

    for(con=0; con<num_contours[slice]; con++)
    {
	startpt = pt = begin_contour[slice][con];
	do
	{
	    ppt = prevpt[pt];
	    dx = savept[pt].x - savept[ppt].x;
	    dy = savept[pt].y - savept[ppt].y;
	    dist = dx*dx + dy*dy;
	    if (dist > maxdist)
	    {
		maxdist = dist;
	    }

	    if (dist > BREAK_LENGTH)
	    {
		newx = savept[pt].x - dx/2.0;
		newy = savept[pt].y - dy/2.0;
		insertpt(newx,newy,slice,con,ppt);
		PRINTF("   Break long segment.\n");
	    }
	    pt = nextpt[pt];
	} while (pt != startpt);
    }
    PRINTF("The max length of contour segment of slice %d is %f\n",
	   slice, (float)sqrt((double)maxdist));
}
#endif

#ifdef MIN_OBTUSE_ANGLE
/* A better solution to eliminate double points. It fails in some rare cases.
 * It searches all obtuse triangles and breaks the opposite contour segment
 * at the perpendicular root point.
 */
int find_obtuse(int slice,CIRTRI *tri,int ntri)
{
    int tr, pt,pt0,pt1,pt2;
    float xp,yp;
    int find = 0;

    for (tr=0; tr<ntri; tr++)
    {
	for (pt=0; pt<3; pt++)
	{
	    pt0 = tri[tr].p[pt];
	    pt1 = tri[tr].p[(pt+1)%3];
	    pt2 = tri[tr].p[(pt+2)%3];
	    if ((tri[tr].rad>MIN_OBTUSE_TRI_RAD) &&
		(pt1==nextpt[pt2] || pt1==prevpt[pt2]) &&
		obtuse(savept[pt0].x,savept[pt0].y,
		       savept[pt1].x,savept[pt1].y, savept[pt2].x,savept[pt2].y,
		       xp,yp))
	    {
		insertpt(xp,yp,slice,savept[pt1].con,(pt1==nextpt[pt2])?pt2:pt1);
		//return 1; /**/
		find = 1;
		break;
	    }
	}
    }

    return(find);
}

void get_perpendicular_root_pt(float x0,float y0,float x1,float y1,float x2,float y2,float &xp,float &yp)
{
    double a,b,c, u,v;

    a = y1-y2;
    b = x2-x1;
    c = (double)x1*(double)y2-(double)x2*(double)y1;
    u = (b*(double)x0-a*(double)y0)/(a*a+b*b); /* a,b can't both be zero */
    v = c/(a*a+b*b);

    xp = u*b - v*a;
    yp = -u*a - v*b;
}

int obtuse(float x0,float y0,float x1,float y1,float x2,float y2,float &xp,float &yp)
{
    double dx0,dy0,dx1,dy1,dx2,dy2,a2,b2,c2;

    dx0 = x1-x2;
    dy0 = y1-y2;
    dx1 = x1-x0;
    dy1 = y1-y0;
    dx2 = x2-x0;
    dy2 = y2-y0;
    a2  = dx1*dx1 + dy1*dy1;	/* square length of one angle edge */
    b2  = dx2*dx2 + dy2*dy2;	/* square length of another angle edge */
    c2  = dx0*dx0 + dy0*dy0;	/* square length of opposite edge of the angle */

    /* use formular c^2 = a^2 + b^2 - 2ab*cos(A) */
    if ((a2+b2-c2)/(2.0*sqrt(a2)*sqrt(b2)) < MIN_OBTUSE_ANGLE)
    {
	get_perpendicular_root_pt(x0,y0,x1,y1,x2,y2,xp,yp);
	/*
	drawline(x0,y0,x1,y1,1,1);
	drawline(x1,y1,x2,y2,1,1);
	drawline(x2,y2,x0,y0,1,1);
	drawline(x0,y0,*xp,*yp,2,1);
	try_fiddle();
	*/
	return 1;
    }
    else return 0;
}
#endif
