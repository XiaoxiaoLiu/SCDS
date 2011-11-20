/* read_contours.c - open and read all contours from a single Volume
 * or Solid into Phase 1 globals. LATER: integrate 2 phases into 1.
 * This routine isolates OS-dependent file mechanisms and
 * site-specific contour formats from the caller.  (gst 24-feb-92)
 *   
 * There are multiple readers here which all read different contour
 * formats, selected at compile-time with define directives (gst
 * 22-oct-92) 
 */

/* Author: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */



/* ****************************************************************
 * this is the UNC-specific contour format (called imex or .con
 * format)
 * **************************************************************** */

#include "cti_globals.h"
#include "stdio.h"

#include <iostream>
using namespace std;

//#ifdef UNC
int read_contours(char* volume_data_name)
{
    FILE *cstrm;
    int trash, x,y,z, pt;
    // unused var // int con;

    cstrm = fopen (volume_data_name, "r");
    if (cstrm == NULL)
    {
	fprintf (stderr, "ERR: Can't open %s\n", volume_data_name);
	perror ("contour file: ");
	return 0;
    }

    if (fscanf (cstrm, "factor %d\n", &trash) != 1)
    {
	fprintf (stderr, "ERR: read_contours: expecting \"factor %%d\"\n");
	return 0;
    }
    if (fscanf (cstrm, "numslices %d\n", &trash) != 1)
    {
	fprintf (stderr, "ERR: read_contours: expecting \"numslices %%d\"\n");
	return 0;
    }

    ncon=0;
    while (fscanf (cstrm, "slice %d\n", &z) == 1)
    {
    	nodeCnt[ncon] = 0;
	do
	{
	    if (fscanf (cstrm, "%d %d\n", &x, &y) != 2)
	    {
		fprintf(stderr, "ERR: contour point format error. Goodbye\n");
		return 0;
	    }
	    if ((x != -1) && (y != -1))
	    {
		/* surface & nodeCnt are GLOBAL */
		pt = nodeCnt[ncon]++;
		surface[ncon][pt][0] = x;
		surface[ncon][pt][1] = y;
		surface[ncon][pt][2] = z;

		if (pt >= MAS_NNODES)
		{
		    fprintf(stderr, "ERR: too many points (%d) on con %d;\n", pt, ncon);
		    fprintf(stderr, "     recompile with higher MAS_NNODES (%d). Goodbye.\n",
			   MAS_NNODES);
		    return 0;
		}
	    }
	} while ((x != -1) && (y != -1));
	ncon++;			/* GLOBAL */
    }

    if (ncon >= MAS_NCONS)
    {
	fprintf(stderr, "ERR: too many contours;\n");
	fprintf(stderr, "     recompile with higher MAS_NCONS (%d). Goodbye.\n",
	       MAS_NCONS);
	return 0;
    }

    fclose (cstrm);
    return 1;			/* aok */
}

#ifdef UNC
//#else

/* ****************************************************************
 * this is the FOUNDATION-specific contour reader
 * **************************************************************** */
/*
#include "cti_globals.h"
#include "rtpt.h"
#include "ids.h"
#include "fou_support_extern.h"

#include "string.h"
#include "stdlib.h"

char *cls[] =
{
    RTPT_VOLUME,
    RTPT_SOLID,
};

*/
int read_contours(char* volume_data_name)
{
    // to read a Foundation Volume or Solid 
    int clindex;
    rtpt_set_t setin[NITEMS(cls)];
    rtpt_set_t setout;
    enum rtpt_status_e ret_status;
    rtpt_string_t volume_data_id = volume_data_name;

    // to find the correct fou obj in the set 
    rtpt_solid_t *sol_p;
    rtpt_volume_t *vol_p;
    rtpt_set_t *node_p;
    int fail = 0;

    // to convert from fou format into Phase 1 format 
    rtpt_contour_t *con_p;
    int con, pt, npt;
    rtpt_set_t *connode_p;
    rtpt_polyline_t *poly_p;
    float z;



    // build set of class names to retrieve 
    for (clindex = 0; clindex < NITEMS(cls); clindex++) {
	    setin[clindex].obj_id.type = RTPT_SET;
	    setin[clindex].set_element_p = malloc(sizeof(rtpt_id_t));
	    ((rtpt_id_t*)(setin[clindex].set_element_p))->type = cls[clindex];
	    setin[clindex].next_p = & setin[clindex+1];
    }
    setin[NITEMS(cls)-1].next_p = 0; // last element has no forward link 

    setout.set_element_p = setout.next_p = 0;

    // open & read data in site-specific format 

    if (! volume_data_id)
	fail = 1;

    if (! fail)
	ret_status = fou_fetchall(volume_data_id, setin, & setout);

    if (   fail
	|| (ret_status != RTPT_ST_OK))
    {
	printf("ERR: fetch returns err (%d) while opening Volume dataset <%s>. Goodbye.\n",
	       ret_status,
	       volume_data_id ? volume_data_id : "nil");
	return 0;
    }

    // find either the first Solid or the first Volume; if neither, 
    // ret failure  
    if (node_p = findnode_by_type(& setout, RTPT_SOLID))
    {
	sol_p = (rtpt_solid_t*)(node_p->set_element_p);
	vol_p = & sol_p->volume;
    }
    else if (node_p = findnode_by_type(& setout, RTPT_VOLUME))
    {
	vol_p = (rtpt_volume_t*)(node_p->set_element_p);
    }
    else
    {
	printf("read_contours: Volume or Solid not found in dataset. Goodbye.\n");
	return 0;
    }


    // convert into Phase 1 format, observing limits on fixed-size arrays 
    
    connode_p = vol_p->contours_p;
    ncon = count_nodes(connode_p); // ncon is GLOBAL 

    if (ncon >= MAS_NCONS)
    {
	printf("ERR: too many contours;\n");
	printf("     recompile with higher MAS_NCONS (%d). Goodbye.\n",
	       MAS_NCONS);
	return 0;
    }

    // copy each pt of each con, observing limits on fixed-size arrays 
    for (con = 0; con < ncon; con++, connode_p = connode_p->next_p)
    {
	con_p = (rtpt_contour_t *) connode_p->set_element_p;
	poly_p = & con_p->poly;
	npt = nodeCnt[con] = poly_p->nvertices;	// nodeCnt is GLOBAL 
	z = poly_p->trans.transform[2][3];

	if (npt >= MAS_NNODES)
	{
	    printf("ERR: too many points (%d) on con %d;\n",
		   npt, con);
	    printf("     recompile with higher MAS_NNODES (%d). Goodbye.\n",
		   MAS_NNODES);
	    return 0;
	}

	for (pt = 0; pt < npt; pt++)
	{
	    // surface is GLOBAL 
	    surface[con][pt][0] = poly_p->vertices_p[pt].x;
	    surface[con][pt][1] = poly_p->vertices_p[pt].y;
	    surface[con][pt][2] = z;
	}
    }
    return 1;
}

#endif
