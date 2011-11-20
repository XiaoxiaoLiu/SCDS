#ifndef _CTI_GLOBAL_H_
#define _CTI_GLOBAL_H_

/* cti_globals.h - globals used in the Contour Tiler; note that there
 * are two styles of globals here, due to the 2 programs that were
 * merged into this program (massage and per).  One style is used
 * during the massage phase of the program, and the other during the
 * tiling phase. LATER, make these into a single style and dynamically
 * allocate them.
 */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

/*
 * The default storage type of these globals is "extern";
 * If the module you are writing actally contains the storage
 * for these globals, do this:
 *   #define CTI_GLOBAL_STORAGETYPE
 *   #include "cti_globals.h"
 */
#ifndef CTI_GLOBAL_STORAGETYPE
#define CTI_GLOBAL_STORAGETYPE extern
#endif

/* **************************************************************** */

#include "cti_gen.h"
#include "cti_contours.h"
#include "cti_defines.h"
#include <vector>

/* Phase 1 globals */
CTI_GLOBAL_STORAGETYPE
    float surface[MAS_NCONS][MAS_NNODES][MAS_NDIMS]; /* Array holding surface */
CTI_GLOBAL_STORAGETYPE
    int   nodeCnt[MAS_NCONS];	/* Nodes in contour */
CTI_GLOBAL_STORAGETYPE
    int   depth[MAS_NCONS];	/* Depth of contour in its slice, used in ph2 too */
CTI_GLOBAL_STORAGETYPE
    int   ncon;			/* # of Contours */

/* **************************************************************** */

/* Phase 2 globals */

CTI_GLOBAL_STORAGETYPE
CTI_POINT savept[CTI_MAX_PTS];  /* all of the points in this */
				/* anatomical structure. */
CTI_GLOBAL_STORAGETYPE
    int numpts,	    /* # of entries used in savept */
        numslices;  /* # of slices (scans) represented by the contours */

/*
 * These arrays are per point:
 *   keep up with its slice and adjacent pts on the same contour
 */
CTI_GLOBAL_STORAGETYPE
    int	slicenum[CTI_MAX_PTS],  
	prevpt[CTI_MAX_PTS],
	nextpt[CTI_MAX_PTS];

/*
 * These arrays are per slice (scan):
 *   keep up with the number of contours on the same slice, and the
 *   index to the point at which each one begins, along with its
 *   length in points.
 */
CTI_GLOBAL_STORAGETYPE
    int num_contours[CTI_MAX_SLICES],
	begin_contour[CTI_MAX_SLICES][CTI_MAX_CONTOURS],
	num_points[CTI_MAX_SLICES][CTI_MAX_CONTOURS];
/*
 * slice numbers might not be consecutive, so keep the correspondence
 *   to slabs, too
 */
CTI_GLOBAL_STORAGETYPE
    int numslabs,  /* half of the # of entries used in slabslice */
	slabslice[2*CTI_MAX_SLICES];

CTI_GLOBAL_STORAGETYPE
    float con_factor,
          zdistance;

CTI_GLOBAL_STORAGETYPE
    char turn45;		/* 0: use orthogonal frame
			         * 1: turn frame 45 degrees
				 */
CTI_GLOBAL_STORAGETYPE
    int framestart;		/* start of 8 frame points: the first */
				/* 4 are orthongonal to the x and y */
				/* axes and the next 4 are tilted 45 */
				/* degrees to the axes. */

/* **************************************************************** */

#include "cti_vedge.h"
#include "cti_tri.h"
#include "cti_tetra.h"


CTI_GLOBAL_STORAGETYPE
    CTI_TETRA tetra[CTI_MAX_TETRAS];
CTI_GLOBAL_STORAGETYPE
    int ntetra;           /* # of entries in tetra array */

#ifndef CTI_NTRI
#define CTI_NTRI 3000
#endif
CTI_GLOBAL_STORAGETYPE
    CTI_TRI tri1[CTI_NTRI], tri2[CTI_NTRI];
  /* how much is enough?  make dynamic */
CTI_GLOBAL_STORAGETYPE
    int ntri1, ntri2;  /* # of entries used in tri1, tri2 */

CTI_GLOBAL_STORAGETYPE
    CTI_VEDGE vedge1[6*CTI_MAX_PTS_PER_SLICE],
              vedge2[6*CTI_MAX_PTS_PER_SLICE];
CTI_GLOBAL_STORAGETYPE
    int nvedge1, nvedge2;   /* # of entries used in vedge1, vedge2 */

CTI_GLOBAL_STORAGETYPE
    char edge1[CTI_MAX_PTS_PER_SLICE][CTI_MAX_PTS_PER_SLICE],
         edge2[CTI_MAX_PTS_PER_SLICE][CTI_MAX_PTS_PER_SLICE];

/* **************************************************************** */

#include "cti_cap.h"

/* is this contour bottom, top, or not capped?
   Indices are [slice][con]
*/
CTI_GLOBAL_STORAGETYPE
    int caps[MAS_NCONS];

#endif // #ifndef _CTI_GLOBAL_H_
