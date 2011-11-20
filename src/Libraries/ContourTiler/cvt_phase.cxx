/* cvt_phase - a total hack to meld the 2 phases together; LATER,
 * I should remove all Phase 1 formats from all routines and then the
 * need for this routine will dissolve.
 */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#include "cti_globals.h"

/* Nomenclature: scan == slice */

int cvt_phase()
{
   int conin,			/* index of con to try to convert */
       conout=0,		/* index of con sucessfully converted */
       connext;			/* index of next signif con (>2 pts), */
				/* if any */
   int pt;			/* index of pt being copied */
   float z;			/* z of current con */
   int con_on_this_scan = 0;	/* index of current, rel to this slice */

   numslices = 1;
   for (conin=0; conin<ncon; conin++)
   {
       if (nodeCnt[conin] < 3)
	   continue;		/* ignore trivial cons */

       begin_contour[numslices-1][con_on_this_scan] = numpts;
       num_points[numslices-1][con_on_this_scan] = nodeCnt[conin];
       depth[conout] = depth[conin];

       z = surface[conin][0][2]; /* z of first pt */
       for (pt=0; pt<nodeCnt[conin]; pt++)
       {
	   savept[numpts].x = savept[numpts].xo = surface[conin][pt][0];
	   savept[numpts].y = savept[numpts].yo = surface[conin][pt][1];
	   savept[numpts].z = z;
	   savept[numpts].con = conout;
	   slicenum[numpts] = numslices-1;

	   if (pt == 0)			/* first pt wraps to last pt */
	       prevpt[numpts] = numpts + nodeCnt[conin] - 1;
	   else
	       prevpt[numpts] = numpts - 1;

	   if (pt == nodeCnt[conin]-1)	/* last pt wraps to first pt */
	       nextpt[numpts] = begin_contour[numslices-1][con_on_this_scan];
	   else
	       nextpt[numpts] = numpts + 1;
	       
	   numpts++;
       }

       /*
	* assuming that cons are sorted in z, ask if the next
        * non-trivial con is on the same slice (same z) as this con;
        * if it is in a diff slice, update the slice info
        */

       /* find next non-trivial con, if any */
       connext = conin+1;
       while (connext < ncon && nodeCnt[connext] < 3)
	   connext++;

       if (   (connext >= ncon)	/* no signif cons left? */
	   || (surface[connext][0][2] != z) /* next signif con on diff slice? */
	   )
       {
	   num_contours[numslices-1] = con_on_this_scan +1;
	   numslices++;
	   con_on_this_scan = 0;
       }
       else
       {
	   con_on_this_scan++;
       }

       conout++;		/* count of non-trivial cons */
   }
   return 0;
}
