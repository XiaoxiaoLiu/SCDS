/* From Mike Holman, St. Louis  -gst (10-Sep-90) */

/***************************************************
 *                                                 *
 * Program: massage                                *
 * Version: 3.0                                    *
 * Purpose: To "correct" various errors in the     *
 *          contour data before the tiling routine *
 *          processes it.  This code will:         *
 *                                                 *
 *               - eliminate duplicate nodes       *
 *                 in each contour.                *
 *               - eliminate crossed contour loops.*
 *               - enforce contour direction       *
 *                                                 *
 * Author:   Michael Holman  (original version)    *
 * Location: Radiation Oncology Center             *
 * Comments: Here is the contour coord. sys. used. *
 *                                                 *
 *                     y | / z                     *
 *                       |/                        *
 *                   ----+---- x                   *
 *                      /|                         *
 *                       |                         *
 *                                                 *
 ***************************************************/

/*
 * NOTES:
 * (1) Some nomenclature:
 *   pt == point == node
 *   contour == slice == scan
 *   anastruct == set of related contours == volume == solid
 *
 * History:
 *   31-Aug-93 gst @ UNC:
 *             - ANSIfy (prototypes, reorder routines, etc)
 *   21-Oct-92 jc @ UNC:
 *             - modification history moved to RELEASE notes
 *   24-Feb-92 gst @ UNC:
 *             - total rewrite; move all I/O outside this file,
 *               because the foundation I/O broke the mold; move
 *               surface, nodeCnt and ncon to globals, and dont pass
 *               them into routines, because passing a huge array
 *               might not work in gcc.
 *             - ANSIfy: add prototypes; routines now have ANSI
 *               parameter lists.  Output routines removed: output to
 *               file no longer done.
 *   06-Aug-91 gst @ UNC:
 *             - moved line segment intersection code to ../util/cross.c 
 *               and ../include/cross.h, because it is general purpose.
 *   18-Dec-90 gst @ UNC:
 *             - fixed bug (twice) in which wrong parameter was passed
 *               to Toggle_Contour_Orient
 *   23-Oct-90 gst @ UNC:
 *             - added STLA (St. Louis ASCII) format support, with ifdef
 *   08-Oct-90 gst @ UNC:
 *             - massive rename of everything to make this production code
 *             - added concept of read/write packets: an extensible
 *               format-dependent I/O data structure
 *             - added GRATIS format support, with ifdef
 *             - retained IMEX format support, with ifdef
 *             - changed intersection routines to test for intersections
 *               BETWEEN contours on the same plane as well as self-crossing
 *               in a single contour
 *             - rewrote intersection routine to handle divide by zero cases
 *             - changed Elim_The_Loops to repeat until ALL loops have
 *               been removed, even those FORMED by deleting pts
 *   18-Sep-90 gst @ UNC:
 *             - fixed bug in Force_Contour_CCW which misinterpreted a CCW contour
 *             - reorganized Force_Contour_CCW to call prev_index &
 *               next-index
 *             - rewrote Input_Surface & Output_Surface to read imex format
 *             - commented out include of gl.h
 *             - created stubs & ifdef mechanism (not code) for sorting
 *               out nested contours
 * 
 * To Do:
 *   - allocation of arrays should be dynamic, so that any size array can
 *     be handled
 */
#include <stdlib.h>
#include "math.h"
#include "cti_globals.h"
#include "cross.h"
#include "cross_extern.h"
#include <stdexcept>

/* defining NESTED_CONTOURS will force the orientations of nested
 * contours to alternate, otherwise, they will all be set to the same
 * orientation (either situation is controlled by parameter want_cw in
 * massage()); when alternating, the outermost con is set to want_cw's
 * value. */

/* define to alternate con orientations, else, all orients are the same  */
#define NESTED_CONTOURS  /**/
#define ONLINE_EPS 0.000001       

      /* * * * * *   l o c a l   r o u t i n e s   * * * * * */

int num_dups, num_loops;

   /* * * * * *   e x t e r n a l   r o u t i n e s   * * * * * */



     /* * * * * *   R e m o v e _ D u p _ N o d e   * * * * * */

static void Remove_Dup_Node(int cNbr,
                     int dupNbr
                     )
{
    int pt;
   
    for (pt=dupNbr; pt<nodeCnt[cNbr]-1; pt++)
    {
        surface[cNbr][pt][0] = surface[cNbr][pt+1][0];
        surface[cNbr][pt][1] = surface[cNbr][pt+1][1];
    }
}

      /* * * * * *   E l i m _ D u p _ N o d e s   * * * * * */

static void Elim_Dup_Nodes()
{
    int con;                  /* Current contour index              */
    int pt1;                  /* Moving base index for each contour */
    int pt2;                  /* Node index                         */

    for (con=0; con<ncon; con++)
        for (pt1=0; pt1<nodeCnt[con]-1; pt1++)
            for (pt2=pt1+1; pt2<nodeCnt[con]; )
            {
                if ((surface[con][pt1][0] == surface[con][pt2][0]) &&
                    (surface[con][pt1][1] == surface[con][pt2][1]))
                {
                    FPRINTF(stderr,
                            "INFO: duplicate node (%f %f) removed from con %d\n",
                            surface[con][pt1][0], surface[con][pt1][1], con);   
                    Remove_Dup_Node(con, pt2);
                    nodeCnt[con] -= 1;
		    num_dups++;
                }
                else
                    pt2++;    /* only incr to next pt when not a dup */
            }
}


           /* * * * * *   I n t e r s e c t   * * * * * */

/* Tests for intersection between 2 separate lines (that is, it is NOT
 * expected that the lines share an endpoint -- if they do, that is a
 * cross)
 *
 * Ret: 1 if the lines cross; else 0
 */

static int Intersect(float from1x, float from1y, /* endpts of the 2 lines */
              float to1x,   float to1y,
              float from2x, float from2y,
              float to2x,   float to2y
              )
{
    int cross_case;

    cross_case = Classify_Intersection(from1x,from1y, to1x,to1y, from2x,from2y, to2x,to2y);

    switch (cross_case)
    {
    case PARALLEL:
    case COLINEAR:
    case NO_CROSS:
        return 0;
    case CROSS:
    case ENDPTCROSS:
    case COINCIDE:
        return 1;
    default:
      throw std::runtime_error("unknown case");
    }
}

      /* * * * * *   S e l f _ I n t e r s e c t   * * * * * */

static int Self_Intersect(int cNbr,     /* Contour NUMber */
                   int firstNode1,
                   int firstNode2,
                   int secondNode1,
                   int secondNode2
                   )
{
    float x1 = surface[cNbr][firstNode1][0],      /* endpts of line 1 */
          y1 = surface[cNbr][firstNode1][1],
          x2 = surface[cNbr][firstNode2][0],
          y2 = surface[cNbr][firstNode2][1],
          x3 = surface[cNbr][secondNode1][0],     /* endpts of line 2 */
          y3 = surface[cNbr][secondNode1][1],
          x4 = surface[cNbr][secondNode2][0],
          y4 = surface[cNbr][secondNode2][1];

    return(Intersect(x1,y1, x2,y2, x3,y3, x4,y4));
}
    /* * * * * *   R e m o v e _ T h e _ L o o p   * * * * * */

static int Remove_The_Loop(int nbrOfNodes, /* Nodes in contour */
                    int cNbr,
                    int baseNode,
                    int crossNode
                    )
{
    int pt1;
    int pt2;
    int nodesToDelete = crossNode - baseNode - 1;

    if (nodesToDelete <= nbrOfNodes - nodesToDelete)
    {
        for (pt1=baseNode+2; pt1<nbrOfNodes-nodesToDelete; pt1++)
        {
            surface[cNbr][pt1][0] = surface[cNbr][pt1+nodesToDelete][0];
            surface[cNbr][pt1][1] = surface[cNbr][pt1+nodesToDelete][1];
        }
    }
    else
    {
        for (pt1=0,pt2=baseNode+2; pt2<=nodesToDelete; pt1++,pt2++)
        {
            surface[cNbr][pt1][0] = surface[cNbr][pt2][0];
            surface[cNbr][pt1][1] = surface[cNbr][pt2][1];
        }
        nodesToDelete = nbrOfNodes - nodesToDelete;
    }
    return(nodesToDelete);
}

      /* * * * * *   E l i m _ T h e _ L o o p s   * * * * * */

static void Elim_The_Loops()
{
    int con,pt1,pt2, loop_found, dn, i;
    
    for (con=0; con<ncon; )
    {
        if (nodeCnt[con] < 3)
        {
            con++;
            continue;         /* con has too few pts -- go on to next con */
        }

        loop_found = FALSE;

        /* look for loops formed by 2 edges which are NOT connected */
        for (pt1=0; pt1<nodeCnt[con]-2; pt1++)
        {
            for (pt2=pt1+2; pt2<nodeCnt[con]-1; )
            {
                if (Self_Intersect(con, pt1, pt1+1, pt2, pt2+1))
                {
                    FPRINTF(stderr,"INFO: loop removed from con %d between edge", 
			    con);
		    FPRINTF(stderr,"%d (%f %f) & edge %d (%f %f)\n",
                            pt1, surface[con][pt1][0], surface[con][pt1][1],
                            pt2, surface[con][pt2][0], surface[con][pt2][1]);
                    nodeCnt[con] -= Remove_The_Loop(nodeCnt[con],con,pt1,pt2);
                    loop_found = TRUE;
		    num_loops++;
                }
                else
                    pt2++;    /* only incr to next pt when not a dup */
            }
        }

        /* special case of above */
        for (pt2=1; pt2<nodeCnt[con]-2;   )
        {
            if (Self_Intersect(con, nodeCnt[con]-1, 0, pt2, pt2+1))
            {
                FPRINTF(stderr, "INFO: loop removed from con %d between edge", con);
		FPRINTF(stderr, "%d (%f %f) & edge %d (%f %f)\n",
                        nodeCnt[con]-1,
                        surface[con][nodeCnt[con]-1][0],
                        surface[con][nodeCnt[con]-1][1],
                        pt2,
                        surface[con][pt2][0],
                        surface[con][pt2][1]);

		dn = pt2;	/* number of nodes to be deleted */
		if (dn <= nodeCnt[con] - dn)
		{
		    for (i=1; i<nodeCnt[con]-1; i++)
		    {
			surface[con][i][0] = surface[con][i+dn][0];
			surface[con][i][1] = surface[con][i+dn][1];
		    }
		}
		nodeCnt[con] -= dn;

                loop_found = TRUE;
		num_loops++;
            }
            else
                pt2++;        /* only incr to next pt when not a dup */

        }

        /* look for a contour which doubles back on itself. */

        /* too much of a good thing... :) */
        if (nodeCnt[con] < 3)
        {
            FPRINTF(stderr, "WARNING: removing loop on con %d leaves insufficient", con);
	    FPRINTF(stderr, " points (%d)\n", nodeCnt[con]);
        }

        /* repeat until clean */
        if (! loop_found)
            con++;
    }
}

/* Tests for intersection between 2 lines which share an endpoint
 *
 * Ret: 1 if the lines cross (other than at the shared endpoint); else 0
 */

static int Intersect_touch(float from1x, float from1y, /* endpts of the 2 lines */
                    float to1x,   float to1y,
                    float from2x, float from2y,
                    float to2x,   float to2y
                    )
{
    int cross_case;

    cross_case = Classify_Intersection(from1x,from1y,
                                       to1x,to1y,
                                       from2x,from2y,
                                       to2x,to2y);

    switch (cross_case)
    {
    case PARALLEL:
    case COLINEAR:
    case NO_CROSS:
    case ENDPTCROSS:
        return 0;
    case CROSS:
    case COINCIDE:
        return 1;
    default:
      throw std::runtime_error("unknown case");
    }
}
 
static int Self_Intersect_touch(int cNbr, /* Contour NUMber */
                         int firstNode1,
                         int firstNode2,
                         int secondNode1,
                         int secondNode2
                         )
{
    float x1 = surface[cNbr][firstNode1][0],      /* endpts of line 1 */
          y1 = surface[cNbr][firstNode1][1],
          x2 = surface[cNbr][firstNode2][0],
          y2 = surface[cNbr][firstNode2][1],
          x3 = surface[cNbr][secondNode1][0],     /* endpts of line 2 */
          y3 = surface[cNbr][secondNode1][1],
          x4 = surface[cNbr][secondNode2][0],
          y4 = surface[cNbr][secondNode2][1];
    // unused vars// float jx, jy;                                         /* junk (ignored) */

    return (Intersect_touch(x1,y1, x2,y2, x3,y3, x4,y4));
}
          
 
   /* * * * * *   F o r c e _ A l l _ C o n t o u r _ O r i e n t  * * * * * */

static int prev_index(int current,
               int nodeCnt)
{
    if (!current)
        return nodeCnt-1;
    return current-1;
}

static int next_index(int current,
               int nodeCnt)
{
    if (current == nodeCnt-1)
        return 0;
    return current+1;
}

/* * * * * *   T o g g l e _ C o n t o u r _ O r i e n t  * * * * * */

static void Toggle_Contour_Orient(int nbrOfNodes,
				  int cNbr
				  )
{
    float temp;
    int pt;

    PRINTF("Toggle orientation of contour %d\n", cNbr); /**/
   
    for (pt=0; pt<nbrOfNodes/2; pt++)
    {
        temp = surface[cNbr][pt][0];
        surface[cNbr][pt][0] = surface[cNbr][nbrOfNodes-pt-1][0];
        surface[cNbr][nbrOfNodes-pt-1][0] = temp;
        temp = surface[cNbr][pt][1];
        surface[cNbr][pt][1] = surface[cNbr][nbrOfNodes-pt-1][1];
        surface[cNbr][nbrOfNodes-pt-1][1] = temp;
    } 
}


#ifndef NESTED_CONTOURS
/* test if a contour is ClockWise: if ccw is true,
 * reverse the order of the pts in the contour
 */
static void Force_All_Contour_Orient(int ccw /* SET: con should be */
				             /* counter-clockwise  */)
{
    int con;

    for (con=0; con<ncon; con++)
    {
        if (nodeCnt[con] < 3) /* ignore malformed contours */
            continue;

        if (Detect_Contour_Orient(con) != ccw)
        {
            FPRINTF(stderr,"INFO: contour #%d was changed to %s\n",
		    con, ccw ? "CCW" : "CW");
            Toggle_Contour_Orient(nodeCnt[con], con);
        }
    }
}
#endif

/* return TRUE: con is cw.   FALSE: con is ccw or not determinable */
static int Detect_Contour_Orient(int con)
{
    float maxY;
    int indexToMaxY, origindexToMaxY;
    int maxYLess1;
    int maxYPlus1;
    int pt;
    float xl,yl,xc,yc,xp,yp;                      /* shorthands: coords for "less",
                                                   * "center" and "more" points.
                                                   */
    float dp;                                     /* dot product */
    
    maxY = surface[con][0][1];
    indexToMaxY = 0;
    for (pt=0; pt<nodeCnt[con]; pt++)
    {
        if (surface[con][pt][1] > maxY)
        {
            indexToMaxY = pt;
            maxY = surface[con][pt][1];
        }
    }
    origindexToMaxY = indexToMaxY;
    maxYLess1 = prev_index(indexToMaxY, nodeCnt[con]);
    maxYPlus1 = next_index(indexToMaxY, nodeCnt[con]);

    xl = surface[con][maxYLess1][0];
    yl = surface[con][maxYLess1][1];
    xc = surface[con][indexToMaxY][0];
    yc = surface[con][indexToMaxY][1];
    xp = surface[con][maxYPlus1][0];
    yp = surface[con][maxYPlus1][1];
    
    /* rotate second vector by 90 degrees, and take dot product */
        
#define DOT_PRODUCT(x1,x2,y1,y2) ((x1)*(x2) + (y1)*(y2))
    /* The dot product is the projection of one 2d vector onto another
     * vector; it results in a scalar which is:
     *   > 0   angle between vectors is acute
     *   == 0  vectors are perpendicular
     *   < 0   angle between vectors is obtuse
     */
    while ((dp = DOT_PRODUCT(yc-yl,xp-xc, xl-xc,yp-yc)) == 0)
    {
        /*
         * while there are 2 co-linear vectors at the Max which touch, we
         * need to "slide" down 1 point to make the judgement as to
         * whether it is a left- or right-hand turn.
         */
        FPRINTF(stderr,
                "INFO: SLIDDING! decision point on contour %d (0-rel) is #%d:\n",
                con, indexToMaxY);/**/
        FPRINTF(stderr,"   dp = %f = %f*%f  +   %f*%f\n",
                dp, yc-yl,xp-xc, xl-xc,yp-yc);/**/
        indexToMaxY = prev_index(indexToMaxY, nodeCnt[con]);
        maxYLess1   = prev_index(indexToMaxY, nodeCnt[con]);
        maxYPlus1   = next_index(indexToMaxY, nodeCnt[con]);
        
        /*
         * test for non-computable contour, to avoid looping
         * forever, looking for the answer.
         */
        if (indexToMaxY == origindexToMaxY)
        {
            fprintf(stderr,
		    "WARNING: Contour %d could not be identified as\n", con);
            fprintf(stderr,
		    "         CW or CCW. Passed through untouched\n");
            return FALSE;
        }
        
        xl = surface[con][maxYLess1][0];
        yl = surface[con][maxYLess1][1];
        xc = surface[con][indexToMaxY][0];
        yc = surface[con][indexToMaxY][1];
        xp = surface[con][maxYPlus1][0];
        yp = surface[con][maxYPlus1][1];
    }

    return (dp > 0);
}


#ifdef NESTED_CONTOURS
 /* * * * * *   C a l c _ C o n t o u r _ O r i e n t   * * * * * */


/* calc the proper orientation of each contour (whether the contour is
 * nested or not) and ensure that the contour actually has this
 * orientation.  Proper orientations are ccw for outermost contours,
 * cw for the contours inside of the outermost, ccw for the contours
 * inside of this second set, and so on, alternating until the
 * innermost contour is reached.  */


/* return TRUE iff contour P surrounds contour C
 * (contours are known NOT to cross)
 */

static int isparent(int conP, int conC)
{
    /* From Preparata & Shamos, "Introduction to Computational
     * Geometry", p42:
     *  L := 0;
     *  l := a horizontal line from some pt S (which is on contour C)
     *       to the left
     *  FOR i := 1 UNTIL numedges[P] DO
     *    IF (edge[i] is not horizontal) THEN
     *      IF (edge[i] intersects l to the left of S at any of l's points
     *         except its lower extreme) THEN L++;
     *  IF (L is odd) THEN z is internal
     *                ELSE z is external
     */

    int L = 0;
    float Sx = surface[conC][0][0],      /* coords of point S */
          Sy = surface[conC][0][1],
          ex1,ey1, ex2,ey2;              /* edge (on P) endpoints */
    double xi;                           /* x-dimension value of intersection pt */
    double a, b, c;			 /* coefficients of test segment */
    int ptP;                             /* starting pt of test contour segment */
    int nextP;                           /* end pt of test contour segment */

    for (ptP=0; ptP < nodeCnt[conP]; ptP++)
    {
        nextP = (ptP+1) % nodeCnt[conP]; /* wrap around at last pt */
        ex1 = surface[conP][  ptP][0];
        ey1 = surface[conP][  ptP][1];
        ex2 = surface[conP][nextP][0];
        ey2 = surface[conP][nextP][1];
        if (( ey1 != ey2 ) &&		
	    ( ey1<Sy && Sy<=ey2 || ey2<Sy && Sy<=ey1 ))
	{
	    /* the segment is not horiz and a cross exists (but the
	     * cross is NOT at the lower extreme)
	     */
	    if (ex1<=Sx && ex2<=Sx)
		L++;			 /* cross is on left side of S */
	    else if (Sx<=ex1 && Sx<=ex2)
		continue;		 /* cross is on right side of S */
	    else
	    {
		/* unfortunately we still have to calculate the
		 * intersection point but we don't have to call
		 * Intersect() because this is not the general case.
		 * [There is a faster way to do this with similiar
		 * triangles] 
		 */
		a = ey1 - ey2;
		b = ex2 - ex1;
		c = (double)ex1*(double)ey2 - (double)ex2*(double)ey1;
		xi = -(b*(double)Sy + c)/a;
		if (FUZZY_ZERO(xi-(double)Sx,ONLINE_EPS))
		    return 1;		 /* (Sx,Sy) is on contour */
		if (xi < (double)Sx)
		    L++;		 /* cross is at left side */
	    }
	}
    }

    return L % 2;			 /* return TRUE if z is internal */
}


void Print_Contour_For_Xgraph(int con)
{
    int pt;
    printf("\"contour %d\n", con);
    for (pt=0; pt<nodeCnt[con]-1; pt++)
    {
	printf("%.1f %.1f\n",surface[con][pt][0],surface[con][pt][1]);
    }
    printf("%.1f %.1f\n\n",surface[con][0][0],surface[con][0][1]);
}

void Print_All_Contours_For_Xgraph()
{
    int con;
    for (con=0; con<ncon; con++)
    {
	Print_Contour_For_Xgraph(con);
    }
}

static void Calc_Contour_Orient(int ccw, float minx[], float maxx[], 
				         float miny[], float maxy[])
{
   int cons_on_this_slice[MAS_NCONS];             /* this secondary */
                                                  /* array keeps track */
                                                  /* of the cons on a */
                                                  /* specific slice */
   int ncons_on_this_slice = 0;                   /* # of contours on the cur slice */
   int ncons_left = ncon;                         /* # left to process */
   int P,C;                                       /* parent & child contours */
   float slicenum[MAS_NCONS];                     /* for now, a slice is */
                                                  /* identified by its */
                                                  /* nominal z position */
   int used[MAS_NCONS];                           /* SET: con has been proc'd */
   int con1, con2, firstcon;
   float slice;
   int Pindex, Cindex;
   int orient;
   // unused var // int pt;
   
   /* init; collect nominal z positions */
   for (con1=0; con1 < ncon; con1++)
   {
       slicenum[con1] = surface[con1][0][2];
       cons_on_this_slice[con1] = 0;
       depth[con1] = 0;
       used[con1] = 0;
   }

   /* each pass of the loop checks all contours which are on the same
    * slice */

   while (ncons_left)
   {
       /* find the first unprocessed contour (on any slice) */
       for (con1=0; con1 < ncon; con1++)
       {
           if (! used[con1])
           {
               slice = slicenum[con1];            /* id for slice */
               firstcon = con1;
               break;
           }
       }

       /* collect (& mark) all contours on this slice */
       ncons_on_this_slice = 0;
       for (con2 = firstcon; con2 < ncon; con2++)
       {
           if (! used[con2] && slicenum[con2] == slice)
           {
               used[con2] = 1;
               cons_on_this_slice[ncons_on_this_slice++] = con2;
           }
       }

       ncons_left -= ncons_on_this_slice;

       /* for each child con: count how many other contours are outside */
       for  (Cindex=0; Cindex < ncons_on_this_slice; Cindex++)
       {
           for (Pindex=0; Pindex < ncons_on_this_slice; Pindex++)
           {
	       P = cons_on_this_slice[Pindex];
	       C = cons_on_this_slice[Cindex];
	       if (P == C || minx[C] > maxx[P] || minx[P] > maxx[C] ||
		             miny[C] > maxy[P] || miny[P] > maxy[C])
		   continue;

	       if (isparent(P,C))
	       {
		   depth[C]++;
		   /*
		     Print_Contour_For_Xgraph(P);
		     Print_Contour_For_Xgraph(C);
		   */
	       }
	   }
       }
       
       /*
       printf("Depths: slice z = %f\n", slice);
       for (con1=0; con1<ncons_on_this_slice; con1++)
       {
	   printf("%d ", depth[con1]);
       }
       printf("\n\n");
       */

       /* if the detected orientations is not the same as the expected */
       /* orientation, toggle the orientation */
       for (Cindex=0; Cindex < ncons_on_this_slice; Cindex++)
       {
           C = cons_on_this_slice[Cindex];
	   orient = Detect_Contour_Orient(C); /* 1 means clockwise */
	   
	   /*
	   printf("CCW? ");
	   orient?printf("yes\n"):printf("no\n");
	   Print_Contour_For_Xgraph(C);
	   getchar();
	   */

	   /* if we consider clockwise as 1 for orient, then
	    * for ccw == 0, nested contours' orientation should be :
	    * depth  0   1   2   3   ...
	    * clock  cw  ccw cw  ccw ...
	    * orient 1   0   1   0   ...
	    *
	    * for ccw == 1, nested contours' orientation should be :
	    * depth  0   1   2   3   ...
	    * clock  ccw cw  ccw cw  ...
	    * orient 0   1   0   1   ...
 	    */
	   if ((!ccw && orient == depth[C] % 2) ||
	       ( ccw && orient != depth[C] % 2))
               Toggle_Contour_Orient(nodeCnt[C], C);
       }
   }

}

/* if two contours intersect, return 1 */
int Contour_Intersect(int con1, int con2)
{
    int pt1, np1, pt2, np2;

    for (pt1=0; pt1<nodeCnt[con1]; pt1++)
    {
	np1 = (pt1+1) % nodeCnt[con1];
	for (pt2=0; pt2<nodeCnt[con2]; pt2++)
	{
	    np2 = (pt2+1) % nodeCnt[con2];
	    if (Intersect(surface[con1][pt1][0], surface[con1][pt1][1],
			  surface[con1][np1][0], surface[con1][np1][1],
			  surface[con2][pt2][0], surface[con2][pt2][1],
			  surface[con2][np2][0], surface[con2][np2][1]))
		return 1;
	}
    }
    return 0;
}
	
/* check if there are two contours intersecting */
int Detect_Contours_Intersection (float minx[], float maxx[], 
				  float miny[], float maxy[])
{
   int cons_on_this_slice[MAS_NCONS];             /* this secondary */
                                                  /* array keeps track */
                                                  /* of the cons on a */
                                                  /* specific slice */
   int ncons_on_this_slice = 0;                   /* # of contours on the cur slice */
   int ncons_left = ncon;                         /* # left to process */
   int P,C;                                       /* parent & child contours */
   float slicenum[MAS_NCONS];                     /* for now, a slice is */
                                                  /* identified by its */
                                                  /* nominal z position */
   int used[MAS_NCONS];                           /* SET: con has been proc'd */
   int con1, con2, firstcon;
   float slice;
   int Pindex, Cindex;
   // unused var // int pt;

   PRINTF("Detecting intersection of contours ... \n"); /**/

   /* init; collect nominal z positions */
   for (con1=0; con1 < ncon; con1++)
   {
       slicenum[con1] = surface[con1][0][2];
       cons_on_this_slice[con1] = 0;
       used[con1] = 0;
   }

   /* each pass of the loop checks all contours which are on the same
    * slice */

   while (ncons_left)
   {
       /* find the first unprocessed contour (on any slice) */
       for (con1=0; con1 < ncon; con1++)
       {
           if (! used[con1])
           {
               slice = slicenum[con1];            /* id for slice */
               firstcon = con1;
               break;
           }
       }

       /* collect (& mark) all contours on this slice */
       ncons_on_this_slice = 0;
       for (con2 = firstcon; con2 < ncon; con2++)
       {
           if (! used[con2] && slicenum[con2] == slice)
           {
               used[con2] = 1;
               cons_on_this_slice[ncons_on_this_slice++] = con2;
           }
       }

       ncons_left -= ncons_on_this_slice;

       /* for each child con: count how many other contours are outside */
       for  (Cindex=0; Cindex < ncons_on_this_slice; Cindex++)
       {
           for (Pindex=Cindex+1; Pindex < ncons_on_this_slice; Pindex++)
           {
	       P = cons_on_this_slice[Pindex];
	       C = cons_on_this_slice[Cindex];
	       if (minx[C] > maxx[P] || minx[P] > maxx[C] ||
		   miny[C] > maxy[P] || miny[P] > maxy[C])
		   continue;
	       if (Contour_Intersect(C,P)) 
	       {
		   //printf("   ERR: Two contours %d and %d intersects.",C,P);
		   //printf(" z = %.1f, %.1f\n", slicenum[C], slicenum[P]);
		   
		   /* throw away the contour with less nodes */
		   //printf("        The contour with less nodes is discarded.\n");
		   if (nodeCnt[P] < nodeCnt[C]) nodeCnt[P] = 0;
		   else 			nodeCnt[C] = 0;

// 		   Print_Contour_For_Xgraph(C);
// 		   Print_Contour_For_Xgraph(P);
// 		   exit(1);	/* that's enough for now... */

	       }
           }
       }
   }
   return 1;
}

/* calculate the bounding box of each contour */
void Get_Contour_Bounding_Box(float minx[], float maxx[], float miny[], float maxy[])
{
    int con,p;
    float x,y;
    
    PRINTF("Calculating bounding box ...\n"); /**/
    for (con=0; con<ncon; con++)
    {
	minx[con] = miny[con] = 10000.0;
	maxx[con] = maxy[con] = -10000.0;
	for (p=0; p<nodeCnt[con]; p++)
	{
	    x = surface[con][p][0];
	    y = surface[con][p][1];
	    if (minx[con] > x)      minx[con] = x;
	    else if (maxx[con] < x) maxx[con] = x;
	    if (miny[con] > y)      miny[con] = y;
	    else if (maxy[con] < y) maxy[con] = y;
	}
    }
}

#endif


                /* * * * * *   m a s s a g e   * * * * * */

void massage(int want_cw	/* SET: want orientation of outermost
                                 * cons to be clockwise; CLR: want
                                 * counter-clockwise */

             )
{
    float minx[MAS_NCONS], maxx[MAS_NCONS], miny[MAS_NCONS], maxy[MAS_NCONS];
    // unused vars//int con,pt;

    /* do the work (using Phase 1 format) */
    num_dups = num_loops = 0;
    Elim_Dup_Nodes();
    Elim_The_Loops();

#ifdef NESTED_CONTOURS
    Get_Contour_Bounding_Box(minx,maxx,miny,maxy);
    Detect_Contours_Intersection(minx,maxx,miny,maxy);
    Calc_Contour_Orient(!want_cw,minx,maxx,miny,maxy);
#else
    Force_All_Contour_Orient(!want_cw);
#endif
    
    printf("%d duplicated points are deleted.\n", num_dups);
    printf("%d contour loops are deleted.\n", num_loops);

    //Print_All_Contours_For_Xgraph(); /**/
}
