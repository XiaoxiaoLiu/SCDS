/* extractslab.c - extract the pts from a particular slab */

/* Author: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#include "cti_globals.h"
#include <stdlib.h>

void extract_slabpt(int slice, int *count_p, int slabpt[]);

/* alloc & init slabpt array; tell how many were alloc'd;
   this array can be used to subset the 2 particular contours in use. */
void extractslab(int loslice, int hislice,
	    int **pslabpt, int *pnumslabpts)
{
    int count,num, con;
    int *slabpt;

    /*cti_debug = 1; */

    num = 0;
    for(con=0; con<num_contours[loslice]; con++)
	num += num_points[loslice][con];
    for(con=0; con<num_contours[hislice]; con++)
	num += num_points[hislice][con];

    if (*pslabpt)
	free(*pslabpt);

    slabpt= (int *) malloc(sizeof(int) * num);
    *pslabpt = slabpt;
    if (slabpt== NULL)
    {
	fprintf(stderr,"extractslab: out of memory (slabpt)\n");
	exit(1);
    }

    /* concatenate slabs together */
    count = 0;
    extract_slabpt(loslice, &count, slabpt);
    extract_slabpt(hislice, &count, slabpt);
    
    if (count != num)
	fprintf(stderr,"extractslab: ERR: count(%d) != num(%d)\n",count,num);
    
    *pnumslabpts = num;
}

void extract_slabpt(int slice, int *count_p, int slabpt[])
{
    int con, pt, startpt;

    PRINTF("extract: starting on slice %d at count=%d\n", slice,
	   *count_p);

    for(con=0; con < num_contours[slice]; con++)
    {
	pt = startpt = begin_contour[slice][con];
	do
	{
	    slabpt[(*count_p)++] = pt;
	    PRINTF("extract:  slabpt[%d] = %d\n",*count_p -1,pt);
	    pt = nextpt[pt];
	} while (pt != startpt);
    }
}
