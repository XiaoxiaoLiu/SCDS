/* cti_globals.c -- this file contains global storage for the massage
 * phase (generated from cti_globals.h) and tiling phases (listed in
 * this file) of the Contour Tiler
 */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#define CTI_GLOBAL_STORAGETYPE
#include "cti_globals.h"

int cti_debug = 0;

/* Phase 1 globals */
//float surface[MAS_NCONS][MAS_NNODES][MAS_NDIMS]; /* Array holding surface */
//int   nodeCnt[MAS_NCONS];	/* Nodes in contour */
//int   ncon = 0;			/* # of Contours */


/* Phase 2 globals */
//CTI_TETRA tetra[CTI_MAX_TETRAS];
//int ntetra=0;           /* # of entries in tetra array */

//CTI_TRI tri1[CTI_NTRI], tri2[CTI_NTRI];
  /* how much is enough?  make dynamic */
//int ntri1, ntri2;  /* # of entries used in tri1, tri2 */

//CTI_VEDGE vedge1[6*CTI_MAX_PTS_PER_SLICE],
          //vedge2[6*CTI_MAX_PTS_PER_SLICE];
//int nvedge1, nvedge2;   /* # of entries used in vedge1, vedge2 */

//char edge1[CTI_MAX_PTS_PER_SLICE][CTI_MAX_PTS_PER_SLICE],
     //edge2[CTI_MAX_PTS_PER_SLICE][CTI_MAX_PTS_PER_SLICE];

