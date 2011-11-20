/* insert pt (x,y) after pt @ index k (on contour c of slice s); phase
 * 2 data
 */

/* Author: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#include "cti_globals.h"
#include <stdlib.h>

void insertpt (float x, float y,
	       int s, int c, int k)
{
    int next = nextpt[k];

FPRINTF(stderr,
	"insertpt: pt %d (%7.3f,%7.3f) after pt %d on contour %d of slice %d\n",
	numpts,x,y,k,c,s);

    /* make sure there is room */
    if (numpts >= CTI_MAX_PTS)
    {
	fprintf(stderr,"insertpt: ERR: too many points IN ALL\n");
	//return; /**/
	exit(1);
    }
    if (num_points[s][c] >= CTI_MAX_ENDS)
    {
	fprintf(stderr,
		"insertpt: ERR: too many points ON THIS CONTOUR; slice %d, con %d\n", s,c);
	fprintf(stderr,"num_points = %d, k=%d, next=%d\n",
		num_points[s][c],k,next);
	// return; /**/
	exit(1);
    }

    savept[numpts].x = savept[numpts].xo = x;
    savept[numpts].y = savept[numpts].yo = y;
    savept[numpts].z = savept[k].z;
    savept[numpts].con = savept[k].con;
    

    slicenum[numpts] = s;
    num_points[s][c]++;

    nextpt[k]= numpts;
    nextpt[numpts] = next;
    prevpt[next] = numpts;
    prevpt[numpts] = k;
    numpts++;
}
