/* mer.c - merge 2 planar Delaunay Triangulation's & their
 * corresponding Voronoi diagrams into a single 3-d triangulation
 */

/* Author: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#include <stdlib.h>
#include "math.h"
#include "errno.h"		/* for errno */
#include "cti_globals.h"
#include "cross_extern.h"
#include "insertpt_extern.h"
#include "del_extern.h"
#include "vor_extern.h"
#include <cfloat>

#define EPS (.00001)		/* threshhold for intersection point and end point */
#define MAX_SHIFT (.001)	/* maximum shift distance for V-diagram */
#define NEAREST_EPS (.001)	/* threshhold for two nearest points */

void shift_vor(int &slice,CTI_TRI tri[],int &ntri,float maxshift);

/* - - - - - - - - - - -   o u t s i d e c o n t o u  r   - - - - - - - - - - - */

static int boi_outside(int i,int j, int n,int p) /* edge from pt i to pt j */
{
    /* Is the segment from i to j inside or outside? */
    /* This assumes points are "close enough" as determined */
    /* by checker */

    /*
     * From Boissonnat's paper, "Shape Reconstruction from Planar
     * Cross Sections", p13:
     * If the vectors ip,ij,in are in counter-clockwise order, then ij
     * is an internal edge. If they are in clockwise order, then ij
     * is an external edge.
     *
     * The paper doesn't mention edges that are coplanar. I'm assuming
     * that this can't happen: (1) the farther points would not have
     * an edge in the Delaunay triangulation because the closer one
     * would (probably) get the edge, and (2) they can't be contour
     * points because the contours cannot have overlapping edges.
     */

    /*
     * These are components of vectors:
     *   ip goes from i to prevpt[i]
     *   in goes from i to nextpt[i]
     *   ij goes from i to j
     */
    double ipx,ipy,inx,iny,ijx,ijy;

    double pa,ja,na;  /* polar angles formed by each point with
		       * point i as the initial point in the vector */

    if (j == p || j == n )   /* is ij "ON" contour? */
	return(FALSE);

    ipx = (double)(savept[p].x - savept[i].x);
    ipy = (double)(savept[p].y - savept[i].y);
    inx = (double)(savept[n].x - savept[i].x);
    iny = (double)(savept[n].y - savept[i].y);
    ijx = (double)(savept[j].x - savept[i].x);
    ijy = (double)(savept[j].y - savept[i].y);

    /* test for failure: p->j->n in CCW order */
    /* Note: the values for p&n can be cached in the POINT struct */
    /* Another way to do this faster is to use cross products of the */
    /* vectors, but speed is not important here, because this routine */
    /* is not called often (i think). */

    /* Note; if atan2 is passed (0,0), bad things will happen (like
     * DOMAIN ERROR) if this occurs, it is because 2 of the points
     * are the same (which is not allowed).
     */

    /* convert to "angle" (polar coords) */
    errno = 0;
    pa = atan2(ipy,ipx);
    if (errno)
    {
	fprintf(stderr,"  ERR: atan2(ip=%f,%f) errno = %d\n", ipy,ipx, errno);
	errno = 0;
    }
    ja = atan2(ijy,ijx);
    if (errno)
    {
	fprintf(stderr,"  ERR: atan2(ij=%f,%f) errno = %d\n", ijy,ijx, errno);
	errno = 0;
    }
    na = atan2(iny,inx);
    if (errno)
    {
	fprintf(stderr,"  ERR: atan2(in=%f,%f) errno = %d\n", iny,inx, errno);
	errno = 0;
    }

    if ( (na==pa) || (pa==ja) || (na==ja) )
    {
	PRINTF("mer.c:boi_outside: ERR: unexpected coplanar edges.\n");
	PRINTF("      : edge: [i-j] = [%d-%d] (with prev=%d next=%d)\n", i,j,p,n);
	PRINTF("      : i(%7.3f,%7.3f) j(%7.3f,%7.3f) p(%7.3f,%7.3f) n(%7.3f,%7.3f)\n",
		savept[i].x, savept[i].y,
		savept[j].x, savept[j].y,
		savept[p].x, savept[p].y,
		savept[n].x, savept[n].y);
	PRINTF("      : pa=%f na=%f ja=%f\n", pa,na,ja);
	PRINTF("      : ja>pa:%d + pa>na:%d + na>ja:%d  sum:%d\n", ja>pa,pa>na,na>ja,
		((ja>pa) + (pa>na) + (na>ja)));
    }

    /* hey! we're adding logicals:
     *   if 2 of the following conditions occur,
     *     then the points are in CCW order
     */
    if ((size_t(ja>pa) + size_t(pa>na) + size_t(na>ja)) == 3)
	PRINTF("mer.c:boi_outside: ERR: unexpected ((ja>pa + pa>na + na>ja) == 3)\n");
    if ( ((ja>pa) + (pa>na) + (na>ja)) == 2)
	return FALSE;		/* boi inside */
    else
	return TRUE;		/* boi outside */
}



int outsidecontour(int i, int j,  /* edge from pt i to pt j */
		   char cw
		   )
{
    int n,p, boi, coni,conj;
    
    n = nextpt[i];
    p = prevpt[i];


    /* these tests are sorted by anticipated frequency of occurence */
    if (n == -1 || p == -1)
	return TRUE;		/* edge contains a frame pt */

    if ( (j==p) || (j==n) )
	return FALSE;		/* edge ON contour */

    if (savept[i].z != savept[j].z)
	return TRUE;		/* edge not contained in a single slice. */

    /* edge connecting two same depth but different contours should be outside */
    coni = savept[i].con;
    conj = savept[j].con;
    if (coni != conj && depth[coni] == depth[conj])
	return TRUE;

    /* if two contours have depths difference 1, then the edge must be internal */
    if (abs(coni-conj) == 1) return FALSE;

    /* if we get here, the edge is either INSIDE or OUTSIDE the
     * contour (not ON). The contours can be either in CCW order or
     * (Boissonnat's) CW order, so reverse the sense of Boissonnat's
     * test routine for CCW contours.  */

    boi = boi_outside(i,j, n,p);
    if (cw)
	return boi;
    return !boi;
}


static void cre8tetra(int p0,int p1,int p2,int p3,
		 CTI_TETRA tetra[],
		 int *ntetra_p, int tettype)
{
    tetra[*ntetra_p].p[0] = p0;
    tetra[*ntetra_p].p[1] = p1;
    tetra[*ntetra_p].p[2] = p2;
    tetra[*ntetra_p].p[3] = p3;
    tetra[*ntetra_p].type = tettype;
    tetra[*ntetra_p].searchnum = 0;
    (*ntetra_p)++;
}



/* Intersect each V-edge in plane a with each V-edge in plane b.
 * Create a T12 tet for each intersection. Perturb tri centers around
 * so that endpoints of the V-edges are not involved in an
 * intersection, which means that either 5 or 6 points in 2 planes are
 * cospherical.
 *
 * NOTE: a side effect of using outsidecontour in this routine is that
 * segments on the frame will be considered outside, and therefore
 * will never be included in any tet.
 *
 * Ret: 0 aok
 *      4 err: all solutions attempted had 4 pts cospherical, so they
 *        were all rejected
 *      5 err: intersection of v-edge was too close to end of v-edge,
 *        which probably means that the solution is bogus
 *      9 err: v-edge reversed; I don't understand why this err
 *        happens, but I think that it is safe to say that it should
 *        result in the slab being aborted.
 */

static int mer_find_in_vor(CTI_TRI tria[],
			   int ntria,
			   CTI_TRI trib[],
			   int ntrib,
			   CTI_TETRA tetra[],
			   int *ntetra_p,
			   float maxper,
			   char cw)
{
    int stra,etra,strb,etrb;	/* Del tri's whose dual V-edges we are intersecting */
    int rota,rotb;		/* which edge of tra or trb we are considering */
    int sa,ea, sb,eb;		/* dual Del edge being considered */
    // unused vas//int t,t1,t2;		/* tri across tra's V-edge of rota */
    float xi,yi;		/* intercept of 2 V-edges */
    float xsa,ysa, xea,yea;	/* start/end of first V-edge */
    float xsb,ysb, xeb,yeb;	/* start/end of second V-edge */
    //unused var//int retval = 1;		/* overall success or failure of the call*/
    //unused var//int rota1, rota2;

    /* for each triangle, store up the T12 tetra's that it's v-edges
     * would produce until it is certain that all of the intersections
     * are valid (that is, that no V-center is too close to another
     * v_center or v_edge in the other plane), and then allow this
     * tri's T12 tet's admission into the master list. This iterative
     * method allows us to displace a v_center and backtrack when
     * invalid intersections are found.  */

    stra=0;

    while (stra < ntria)
    {
	xsa = tria[stra].xc; ysa = tria[stra].yc;
	for (rota=0; rota<3; rota++)
	{
	    /* process in forward direction only */ 
	    if (tria[stra].npts[rota] != 3) continue; 

	    etra = tria[stra].t[rota];
	    if (etra == -1)	/* real edge? (is this test necessary?) */
		continue;
	    
	    xea = tria[etra].xc; yea = tria[etra].yc;
	    sa = tria[stra].p[(rota+1)%3];
	    ea = tria[stra].p[(rota+2)%3];

	    if ( (xsa==xea) && (ysa==yea) )
	    {
		printf("  WARNING: skipping 0 length v-edge @ TRIa ");
		printf("%d =%d-%d-%d, edge %d-%d.\n",
		       stra,
		       tria[stra].p[0], tria[stra].p[1], tria[stra].p[2],
		       tria[stra].p[(rota+1)%3], tria[stra].p[(rota+2)%3]);
		tria[stra].npts[rota] = 0;  /* ignore in future passes */
		//drawpt(xsa,ysa, merwindow, 1);/**/
		continue;
	    }
	    
	    for (strb=0; strb < ntrib; strb++)
	    {
		xsb = trib[strb].xc; ysb = trib[strb].yc;
		for (rotb=0; rotb<3; rotb++)
		{
		    if (trib[strb].npts[rotb] != 3) continue;

		    etrb = trib[strb].t[rotb];
		    if (etrb == -1)	/* real edge? (is this test necessary?) */
			continue;
	    
		    xeb = trib[etrb].xc; yeb = trib[etrb].yc;
		    sb = trib[strb].p[(rotb+1)%3];
		    eb = trib[strb].p[(rotb+2)%3];

		    if ( (xsb==xeb) && (ysb==yeb) )
		    {
			printf("  WARNING: skipping 0 length v-edge @ TRIb ");
			printf("%d =%d-%d-%d, edge %d-%d.\n",
			       strb,
			       trib[strb].p[0], trib[strb].p[1], trib[strb].p[2],
			       trib[strb].p[(rotb+1)%3], trib[strb].p[(rotb+2)%3]);
			trib[strb].npts[rotb] = 0;  /* ignore in future passes */
			//drawpt(xsb,ysb, merwindow, 1);/**/
			continue;
		    }
		    
		    PRINTF("working on %d-%d & %d-%d\n", sa,ea, sb,eb);
		    if (cross_exact(xsa,ysa, xea,yea, xsb,ysb, xeb,yeb, &xi,&yi))
		    {
		      //drawpt(xi,yi,3,5); /**/

			/*
			if (outsidecontour(sa,ea, cw) || outsidecontour(sb,eb, cw))
			{
			    PRINTF("tet OUTSIDE - eliminated!");
			    continue;
			}
			*/

			/* see if the endpt of either V-edge is too close to the
			 * intercept pt, using a poor man's distance formula
			 * (the Manhattan distance, L1, is a square instead of a
			 * circle -- at this tiny distance, the difference is
			 * not important) */

			if ((FUZZY_ZERO(xsa-xi, EPS) && FUZZY_ZERO(ysa-yi, EPS)) ||
			    (FUZZY_ZERO(xea-xi, EPS) && FUZZY_ZERO(yea-yi, EPS)) ||
			    (FUZZY_ZERO(xsb-xi, EPS) && FUZZY_ZERO(ysb-yi, EPS)) ||
			    (FUZZY_ZERO(xeb-xi, EPS) && FUZZY_ZERO(yeb-yi, EPS)))
			{
			    PRINTF("  Too close to one of end points of v-edges.\n");
			    return 5;
			}
			/*
			drawpt(xi,yi, 1, 14);
			char text[100];
			sprintf(text, "%d-%d %d-%d", sa,ea,sb,eb);
			drawtext(xi,yi, text, 1);
			*/

			cre8tetra(sa, ea, sb, eb, tetra,ntetra_p, CTI_T12);

		    } 				    /* if (cross(xsa,...)) */
		}				    /* for (rotb=0;...) */
	    }					    /* for (strb=0;..). */
	}					    /* for (rota=0;...) */
	stra++;
    }						    /* while (stra < ntria) */
    return 0;
}

/* replace pt_closest by a list of indices of points which are closest
 * *and* equidistance to pt (x,y), among all contour points in this plane */

static void nearest_pts(float x,float y,
			int slice,		/* slice to search */
			int pt_closest[],	/* RET */
			int *npts_p		/* RET: # of pts in pt_closest */
			)
{
    int firstpt, npts, t2, con;
    float dist, dist_closest, x1, y1;
    // !!! added by bcd, is this right ????
    // !!! dist_closest must be initialized, but to what
    // !!! im guessing something big
    // dist_closest = FLT_MAX; replaced by (npts == 0)

    npts=0;
    for (con=0; con < num_contours[slice]; con++)
    {
      t2 = firstpt = begin_contour[slice][con];
      do
	{
	  x1 = savept[t2].x; y1 = savept[t2].y;
	  dist = (x1-x)*(x1-x) + (y1-y)*(y1-y);
	  if (npts == 0)	// list is empty so add pt to list
	    {
	      dist_closest = dist;
	      pt_closest[npts++] = t2;
	    }
	  else
	    {
	      // t2's dist is equidistant, closer, or farther. Since
	      // equidistant is a fuzzy test, it includes part of
	      // "closer" and "farther" and is more important than them.
	      
	      if ( FUZZY_ZERO(dist - dist_closest, NEAREST_EPS) ) // equidistant?
		{
		  // this test needs to come before dist_closest test because the tests overlap.
		  // this point is nearly as close as current closest, add it to the list.
		  // do not update dist_closest - set it by the first point only.
		  pt_closest[npts++] = t2;
		}
	      else if ( dist < dist_closest )	      // closer?
		{
		  // we found a closer point, start a new list and set a new distance
		  npts = 0;
		  pt_closest[npts++] = t2;
		  dist_closest = dist;
		}
	      else
		{
		  // otherwise this point is too far away, forget about it
		}
	    }
	  
	  t2 = nextpt[t2];
	} while (t2 != firstpt);
    }

#ifdef NOT_DEFINED
    //THIS NEEDS SOME ANALYSIS AS TO 1) WHAT IT'S DOING and 2) IS IT DOING IT RIGHT?
    // bcd & gregg 02July03
    for (int framept = framestart + (turn45*4);
	 framept < framestart + (turn45*4) + 4;
	 framept++)
    {
	x1 = savept[framept].x; y1 = savept[framept].y; 
	dist = (x1-x)*(x1-x) + (y1-y)*(y1-y);
	if ( (dist <= dist_closest) || (npts == 0) )
	{
	    if ( FUZZY_ZERO(dist - dist_closest, .001) && (npts == 0) ) 
	    {
		/* add pt to list */
		pt_closest[npts++] = t2;
	    }
	    else
	    {
		/* empty the list; add pt to list */
		dist_closest = dist;
		npts = 0;
		pt_closest[npts++] = t2;
	    }
	}
    }
#endif

    *npts_p = npts;
}



/*
Find the closest contour pt to each Del Tri's circumcenter. There are
2 ways to do this:

1) for each tri:
     find the closest pt

2) for each pt:
     find the closest tri


I used (1). BUT, (2) can run faster, because the tri's know their
neighbors, and can always advance towards a neighbor which is closer
to the target. LATER: see comment about (2) at end of file
*/

static int
    mer_find_in_del(CTI_TRI tria[],
		    int ntria,	       /* # of tri's in slice a */
		    int sliceb,	       /* slice number of slice b */
		    CTI_TETRA tetra[], /* RET: additional T1 or T2 tet's */
		    int *ntetra_p,     /* # of tets in tetra[] */
		    float maxper,
		    char cw,	       /* relevant? */
		    int tettype	       /* tells which type of tet's we should */
				       /* create: CTI_T1 or CTI_T2 */
		    )
{
    int s,e,t,pt_closest[100], tra, done, npts;
    float xc,yc;
    
    for (tra=0; tra < ntria; tra++)
    {
	if (tria[tra].outside) continue;

	xc = tria[tra].xc;   yc = tria[tra].yc;

	/* if 5 or more pts (3 on this plane and 2 or more on the
	 * other plane) pts are cospherical, perturb the circumcenter
	 * of this Del Tri. */

	done = 0;
	while (! done)
	{
	    nearest_pts(xc, yc, sliceb, pt_closest, & npts);
	    done = 1;

	    if (npts == 1)
	    {
		done = 1;
	    }
	    else
	    {
		printf("  WARNING: Moving center of tri %d. No safety check!\n",tra);
		xc = tria[tra].xc + maxper * ((float)(rand() % 32675) / 32676.);
		yc = tria[tra].yc + maxper * ((float)(rand() % 32675) / 32676.);
	    }
	    /**/
	}
	tria[tra].xc = xc;
	tria[tra].yc = yc;

	PRINTF("TRI @ s=%d e=%d t=%d  center=(%f %f)\n",s,e,t,xc,yc);

	/* if there is more than 1 closest point, then we have a
	 * number of cospherical points, so we need to move some of
	 * them. */

	PRINTF("closest pt = %d  (dist=unknown)\n", pt_closest[0]);

	s = tria[tra].p[0];
	e = tria[tra].p[1];
	t = tria[tra].p[2];
	    
	/*
	if (outsidecontour(s,e, cw) ||
	    outsidecontour(e,t, cw) ||
	    outsidecontour(t,s, cw) )
	{
	    PRINTF("tet OUTSIDEa - eliminated!\n");
	    continue;
	}
	*/
	//drawpt(xc,yc, merwindow, 10+tra);/**/

	cre8tetra(s, e, t, pt_closest[0], tetra,ntetra_p, tettype);
    }

    return 1;
}


int mer(CTI_TRI tria[],
	int ntria,
	int slicea,
	CTI_TRI trib[],
	int ntrib,
	int sliceb,
	CTI_TETRA tetra[],
	int *ntetra_p,
	float maxper,
	char cw  /* TRUE: contour points are considered in clockwise order, else CCW */
	)
{
    int retval;

    /* the v-edge intersection test should probably be done before the
     * other 2 tests, because the former test might perturb some
     * v-centers to arrive at a unique solution, which might affect
     * the latter test.
     */

 restart:
    *ntetra_p = 0;

    PRINTF("  INTERSECT V-EDGES...");
    retval = mer_find_in_vor(tria,ntria,
			     trib,ntrib,
			     tetra,ntetra_p,
			     maxper, cw);
    PRINTF("%d tet's produced so far.\n",*ntetra_p);

    if (retval == 5) {
#ifdef MAX_SHIFT
	PRINTF("  Retry intersecting ! ...\n");
	shift_vor(slicea,tria,ntria,float(MAX_SHIFT));
	goto restart;
#else
    return retval;
#endif
    }

    PRINTF("  CREATE T1 Tetras in SLICE 1...");
    retval = mer_find_in_del(tria,ntria, sliceb, tetra,ntetra_p, maxper, cw, CTI_T1);
 
    if (retval != 1) return 1;
    PRINTF("%d tet's produced so far.\n",*ntetra_p);

    PRINTF("  CREATE T2 Tetras in SLICE 2...");
    retval = mer_find_in_del(trib,ntrib, slicea, tetra,ntetra_p, maxper, cw, CTI_T2);

    if (retval != 1) return 2;
    PRINTF("%d tet's produced so far.\n",*ntetra_p);

    if (*ntetra_p >= CTI_NTRI)
    {
	printf("  ERR: tetra[] overflow (%d) \n",*ntetra_p);
	retval = -1;
    }
    return 0;
}

#ifdef MAX_SHIFT
void shift_vor(int &slice,CTI_TRI tri[],int &ntri,float maxshift)
{
    int con, pt, firstpt, t;
    float delta_x, delta_y;

    delta_x = maxshift * ((float)(rand() % 997) / 997.);
    delta_y = maxshift * ((float)(rand() % 997) / 997.);
    /* shift points */
    for (con=0; con < num_contours[slice]; con++)
    {
	pt = firstpt = begin_contour[slice][con];
	do
	{
	    savept[pt].x += delta_x;
	    savept[pt].y += delta_y;
	    pt = nextpt[pt];
	} while (pt != firstpt); 
    }
    /* shift centers of triangles */
    for (t=0; t<ntri; t++) {
	tri[t].xc += delta_x;
	tri[t].yc += delta_y;
    }
}
#endif

/*
Rewrite mer_find_in_del like this:

1) start at any pt;
2) traverse adj graph (of del tri);
3) for each pt connected to this pt:
   - if dist from other plane's pt is less than running min,
     replace min by this pt`s dist
     count = 0
   - else if dist is equal,
     if "equal" min == this min,
     count++
4) switch (count)
   0: move to new pt
   1: on line: move pt
*/
