/* sculpt.c - sculpt from the Delauney triangulation (really,
 *   tetrahedrulation) to leave the maximal (largest) solid object.
 *
 * The procedure called "correct" mentioned in the Boissonnat
 * has not yet been implemented. It may improve the results.
 *
 * 01-Dec-89 gst original implementation of Boissonnat's algorithm
 * 04-Oct-90 gst integrate sculpting & tetrahedrulation routines
 *               into a single program
 */

/* Author: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#include <stdlib.h>
#include "cti_globals.h"
#include <stdexcept>

static int oppdir(int dir)
{
    switch (dir)
    {
    case 0: return 1;
    case 1: return 0;
    case 2: return 3;
    case 3: return 2;
    default: 
      {
	throw std::runtime_error("invalid case");
      }
    }
}


/* - - - - - - - - -  s o l i d  - - - - - - - - - -  */

/* find solidity of this T12 tet (Rule 2)
 * Ret: 0 not solid
 *      1 solid
 *      2 cycle err
 */
static int solid(int nt,	/* starting tet */
		 int dir	/* dir of neigh tet to search */
		 )
{
  //unused var//int diff;
    int newdir;			/* dir along a line of T12 tet's */
				/* relative to a specific tet  */
    int neigh;			/* nt's neighboring tet */
    static int searchnum = 1;	/* unique id of this search */

    PRINTF(" classifying solidity of tet #%d in dir %d\n", nt, dir);

    searchnum++;		/* start a new search */

    if (tetra[nt].type & (CTI_ELIM))
    {
	PRINTF("  #%d is ELIM\n",nt);
	return 0;   /* can't possibly be in solid */
    }

    /*
     * look at my neighbors along the spec'd dir to determine my solidity;
     * LATER: recode this in better style.
     */
    while (1)
    {

	/* If I'm solid, return good */
#define THE_RIGHT_WAY
#ifdef THE_RIGHT_WAY
	if ( (dir == 0 || dir == 1)
	    && (tetra[nt].type & CTI_T2))
	{
	    PRINTF("  SOLID: tet %d is T2\n", nt);
	    return 1;
	}
	else if ( (dir == 2 || dir == 3)
	    && (tetra[nt].type & CTI_T1))
	{
	    PRINTF("  SOLID: tet %d is T1\n", nt);
	    return 1;
	}
#else
	if (tetra[nt].type & (CTI_T1 | CTI_T2))
	{
	    PRINTF("  SOLID: tet %d is T2/T1\n", nt);
	    return 1;
	}
#endif

	/* if I have no neighbor, or my neigh has been elim'd, return bad */
	neigh = tetra[nt].t[dir];
	if ((neigh == -1) || (tetra[neigh].type & CTI_ELIM))
	{
	    PRINTF("  NOT solid: tet %d %s %d (neigh:%d)\n",
		   nt,
		   (neigh == -1) ?
		     "has no neigh in dir" :
		     "has ELIM'd neigh in dir",
		   dir, neigh);
	    return 0;
	}

	/* if my neighbor has been searched in THIS call of solid(), */
	/* then there is a loop, so return err */
	if (tetra[neigh].searchnum == searchnum)
	{
	    return 2;
	}
	else
	    tetra[neigh].searchnum = searchnum;	/* mark as searched */

	/* find the current direction in terms of my neighbor's points */
	if ( dir == 0 || dir == 1 )
	{
	    if (tetra[neigh].p[0] == tetra[nt].p[oppdir(dir)])
		newdir = 0;
	    else
		newdir = 1;
	}
	else if ( dir == 2 || dir == 3 )
	{
	    if (tetra[neigh].p[2] == tetra[nt].p[oppdir(dir)])
		newdir = 2;
	    else
		newdir = 3;
	}

	/* move to my neighbor, and use his reference model */
	PRINTF("  moving to tet %d, dir %d\n", neigh, newdir);
	nt = neigh;
	dir = newdir;
    }
}

/* search & destroy - search along the specified dir's, eliminating all
 * T12 tet's until no more T12's are found. */
static int reject(int nt,
		  int dir1,	/* dir of neigh tet to search */
		  int dir2	/* other dir of neigh tet to search */
		  )
{
  // unused var // int diff;
    int newdir;			/* dir along a line of T12 tet's */
				/* relative to a specific tet  */
    int ntsave = nt;

    PRINTF("  rejecting tets connected to #%d in dir %d-%d\n", nt,dir1,dir2);

    /*
     * look at my neighbors along the spec'd dir to determine my solidity;
     * LATER: recode this in better style.
     */
    while (1)
    {
	if (tetra[nt].type & (CTI_T1 | CTI_T2))
	    break;
	
	tetra[nt].type |= CTI_ELIM;

	if (tetra[nt].t[dir1] == -1)
	    break;

	if ( dir1 == 0 || dir1 == 1 )
	{
	    if (tetra[tetra[nt].t[dir1]].p[0] == tetra[nt].p[oppdir(dir1)])
		newdir = 0;
	    else
		newdir = 1;
	}
	else if ( dir1 == 2 || dir1 == 3 )
	{
	    if (tetra[tetra[nt].t[dir1]].p[2] == tetra[nt].p[oppdir(dir1)])
		newdir = 2;
	    else
		newdir = 3;
	}

	/* move to my neighbor, and his dir */
	nt = tetra[nt].t[dir1];
	dir1 = newdir;
    }

    nt = ntsave;		/* reset to original tet */
    while (1)
    {
	if (tetra[nt].type & (CTI_T1 | CTI_T2))
	    break;
	
	tetra[nt].type |= CTI_ELIM;

    //It is not clear what should be returned in this case.
	if (tetra[nt].t[dir2] == -1)
	    return -2;

	if ( dir2 == 0 || dir2 == 1 )
	{
	    if (tetra[tetra[nt].t[dir2]].p[0] == tetra[nt].p[oppdir(dir2)])
		newdir = 0;
	    else
		newdir = 1;
	}
	else if ( dir2 == 2 || dir2 == 3 )
	{
	    if (tetra[tetra[nt].t[dir2]].p[2] == tetra[nt].p[oppdir(dir2)])
		newdir = 2;
	    else
		newdir = 3;
	}

	/* move to my neighbor, and his dir */
	nt = tetra[nt].t[dir2];
	dir2 = newdir;
    }
}


/* - - - - - - - - -  c l a s s i f y T 1 2 s e t  - - - - - - - - - -  */
/* Ret: 0 err
 *      1 aok
*/
static int classifyT12set()
{
    int nt;
    int got_one = 1;		/* eliminated a tet during this pass */
    int s1,s2;			/* retval's from solid() */

    /* continue searching and rejecting tets until we get an entire
     * pass in which no tets are rejected. */
    while (got_one)
    {
	got_one = 0;
	for (nt=0; nt<ntetra; nt++)
	{
	    if (! (tetra[nt].type & (CTI_T1 | CTI_T2 | CTI_ELIM)))
	    {
		PRINTF("classifyT12set: starting on tet #%d {neigh= #%d #%d #%d #%d}\n",
		       nt,
		       tetra[nt].t[0],
		       tetra[nt].t[1],
		       tetra[nt].t[2],
		       tetra[nt].t[3]);

		/* elim tet's which are NOT solid in either direction
		 * and are colinear (?) */
		if (! (s1=solid(nt,0)) && ! (s2=solid(nt,1)))
		{
		    if (s1 == 2 || s2 == 2)
			return 0; /* err */
		    PRINTF("  NOT solid in dir 0-1\n");
		    got_one = 1;
		    /*reject(nt,0,1);*/
		    PRINTF("  REJECT (plane 1): %d\n", nt);
		    tetra[nt].type |= CTI_ELIM;
		}
		//
		// bcd replaced the following line
		// so that s1 and s2 would get initialized
		// 
		// its not clear that this is the correct thing to do
		// 
		// if (! solid(nt,2) && ! solid(nt,3))
		//
		if (! (s1=solid(nt,2)) && ! (s2=solid(nt,3)))
		{
		    if (s1 == 2 || s2 == 2)
			return 0; /* err */
		    PRINTF("  NOT solid in dir 2-3\n");
		    got_one = 1;
		    /*reject(nt,2,3);*/
		    PRINTF("  REJECT (plane 2): %d\n", nt);
		    tetra[nt].type |= CTI_ELIM;
		}
	    }
	}
    }
    return 1;			/* success */
}



#ifdef NOT_USED
/* find the face which points to me in my neighbor's tet */
static int findtet(int nt, int neit)
{
    int face;
    for (face=0; face<4; face++)
    {
	neit = tetra[nt].t[face];
	if (neit != -1 && tetra[neit].t[face] == nt)
	    return face;
    }
    printf("findtet: err\n");
}
#endif


/* - - - - - - - - -  c l a s s i f y T 1 T 2 - - - - - - - - - -  */
/* find solidity of this T1 or T2 tet */
static void classifyT1T2()
{
#define NTETLIST 300		/* alloc chunk size for tet list */
    int *tetlist;		/* list of T1/T2 tets to be searched */
    int ntetlist;		/* # of alloc`d entries in tet list */

    int start, stop;		/* indexes into tetlist: start is */
				/* the current tet, and stop */
				/* is the end of the list */
    int neigh, neighindex;
    int keep;			/* set: found a valid tet */
    int tetmark;
    int tet;			/* seed tet */

    tetlist = (int*) malloc(NTETLIST * sizeof(int));
    ntetlist = NTETLIST;
    
    for (tet=0; tet < ntetra; tet++)
    {
	/* Rule 3 is for non-eliminated T1/T2, so ignore elim'd tet's
	 * and ignore T12 tet's
	 */
	if ((tetra[tet].type & (CTI_ELIM | CTI_MARK)) ||
	    (!(tetra[tet].type & (CTI_T1 | CTI_T2))) )
	    continue;

	/* set up the search table to have a single entry -- tet */
	start = 0;
	stop = 1;
	tetlist[0] = tet;

	PRINTF("START @ tet %d: found: ", tet);

	tetra[tet].type |= CTI_MARK;
	keep = FALSE;		/* assume false, search for true */

	/* while there are more tet's to search, apply Rule 3.
	 *   Mark a tet as searched if the result of this search will
	 *   apply to it.
	 */
	do
	{
	    /* for MOST neigh's of this tet [Note that only [0,1,2] are
	     * used ... 4 will not occur since a t1/t2 tet does not have a
	     * neigh at this position.] */
	    for (neighindex=0; neighindex<3; neighindex++)
	    {
		neigh = tetra[tetlist[start]].t[neighindex]; /* get neigh tet */
		
		/* if a neigh tet exists, and has been neither elim'd nor searched,
		 *    if neigh is T12,
		 *       keep the whole list of tet's
		 *    else
		 *       mark tet as searched
		 *       add to end of list
		 */
		if (   (neigh != -1)
		    && !(tetra[neigh].type & (CTI_ELIM | CTI_MARK)) )
		{
		    if (!(tetra[neigh].type & (CTI_T1 | CTI_T2)))
			keep = TRUE; /*do not reject*/
		    else
		    {
			/* new T1 or T2 (fresh meat!) -- record for future search*/
			tetra[neigh].type |= CTI_MARK;
			if (stop == ntetlist)
			{
			    ntetlist += NTETLIST;
			    tetlist = (int*)realloc(tetlist,ntetlist * sizeof(int));
			}
			tetlist[stop] = neigh;
			stop++;
		    }
		}
	    }
	    start++;
	} while(start < stop /* &&!keep */);
	
	/* if no valid neigh of the group was found,
	 *    mark every tet in the group as "searched and elim'd"
	 */
	if (! keep)
	{
	    for (tetmark=0; tetmark < stop; tetmark++)
	    {
		PRINTF("%d ", tetlist[tetmark]);
		tetra[tetlist[tetmark]].type |= (CTI_MARK | CTI_ELIM);
	    }
	    PRINTF(" (ELIM'd)\n");
	}
	else
	    PRINTF("(ALL KEPT)\n");
    }

    free(tetlist);
}



/* - - - - - - - - -   s c u l p t  - - - - - - - - - - - - - -  */
/* Ret: 0 err
 *      1 aok
*/
int sculpt()
{
    cti_debug = 0;

    if (!classifyT12set())
	return 0;
    classifyT1T2();
    return 1;
}


