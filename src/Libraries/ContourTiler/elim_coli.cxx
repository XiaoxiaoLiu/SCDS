/* elim_coli.c - eliminate some pts, so that resultant contour has no
 * consecutive colinear pts */

/* Author: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#include "cti_globals.h"
#include "coli_extern.h"

/* macro to delete a point */
#define DELPT(PT)                \
	    nextpt[PT] = -1;     \
	    prevpt[PT] = -1;     \
	    savept[PT].con = -1; \
	    slicenum[PT] = -1;

/* delete the pt from the contour; if there are less than 3 pts left
 * in the con after deletion, delete the entire con. CAREFUL: can
 * change all globals concerning contours, such as: savept,
 * begin_contours, num_contours, num_points, slicenum, nextpt,
 * prevpt...  RET: 1 if con was deleted, else 0 */

int deletept(int s, int c, int pt)
{
    int next = nextpt[pt];
    int prev = prevpt[pt];
    int con, firstpt;
    int allgone = 0;		/* all pts were deleted from con */
    int delpt, nextdelpt;

    if (next == -1 || prev == -1)
    {
	fprintf(stderr, "ERR: deletept: pt %d is not on any con\n", pt);
	return 0;		/* pt is no good */
    }

    num_points[s][c]--;
    if (num_points[s][c] < 3)	   /* are all pts on this contour gone? */
    {
	fprintf(stderr,
		"deletept: ERR: too few points ON THIS CONTOUR; slice %d, con %d\n",
		s,c);
	fprintf(stderr,"num_points = %d, pt=%d, next=%d\n",
		num_points[s][c],pt,next);
	allgone = 1;
    }

    if (allgone)
    {
	/* del all pts on contour (so much code -- so little function) */
	delpt = firstpt = begin_contour[s][c];
	do
	{
	    nextdelpt = nextpt[delpt];
	    DELPT(delpt);
	    delpt = nextdelpt;
	} while (delpt != firstpt); 

	/* delete contour from slice */
	for (con = c; con < ((num_contours[s] - c) - 1); con++)
	{
	    /* move remaining cons down 1 place */
	    num_points[s][con] = num_points[s][con+1];
	    begin_contour[s][con] = begin_contour[s][con+1];
	}

	num_contours[s]--;
    }
    else
    {
	/* remove the pt from the contour */
	DELPT(pt);
	nextpt[prev] = next;	/* reroute con around this pt */
	prevpt[next] = prev;
	if (begin_contour[s][c] == pt)
	    begin_contour[s][c] = next;
    }

    return allgone;
}

/* elim coli pts on the spec'd slice */
void elim_coli(int slice	/* slice to eliminate colinear pts on */
	       )
{
    int pt, npt, nnpt;		/* pt, next pt, next next pt */
    int con;			/* first pt in contour */
    int allgone;		/* 1: all pts in a con were colinear */

    /* CAREFUL: num_contours, nextpt & begin_contour may be changed by
     * deletept(), so don't cache their values during loops.
     * In each iteration of the loop, either con will be incremented
     * or num_contours[] will be decremented */

    con = 0;
    while (con < num_contours[slice])
    {
	/* if pt does not exist, con has been deleted */
	if ((pt = begin_contour[slice][con]) != -1)
	{
	    con++; continue;
	}

	do
	{
	    npt = nextpt[pt];
	    nnpt = nextpt[npt];
	    allgone = 0;
	    if (coli(savept[  pt].x, savept[  pt].y,
		     savept[ npt].x, savept[ npt].y,
		     savept[nnpt].x, savept[nnpt].y))
	    {
		printf("elim_coli: deleting pt %d (slice %d, con %d)\n",
		       npt, slice, con);
		allgone = deletept(slice,con,npt); /* delete npt */
	    }

	    if (! allgone)
		pt = nextpt[pt];   /* goto next pt, skipping del'd pt */
	    else
		break;		   /* con deleted -- goto next con */
	} while (pt != begin_contour[slice][con]); 
	if (! allgone)
	    con++;
    }
}

/* elim coli pts on all slices */
void elim_coli_all()
{
    int slice;
    for (slice=0; slice < numslices; slice++)
    {
	if (num_contours[slice] != 0)
	{
	    elim_coli(slice);
	}
    }
}
