/* ctiler.c - the Contour TIler (CTI); from a set of cross-sections
 * (contours) of a single volume, this will produce a set of
 * triangular tiles representing the surface of the volume.
 *
 * Note: main should call ctiler() [see bottom of this file]
 */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

#include "cti_globals.h"
#include "rtpt.h"
#include "cross_extern.h"
#include "cvt_phase_extern.h"
#include "del_extern.h"
#include "mer_extern.h"
#include "elim_coli_extern.h"
#include "addrand_extern.h"
#include "sculpt_extern.h"
#include "surfacetotile_extern.h"
#include "slabtocaptile_extern.h"
#include "vor_extern.h"
#include "massage_extern.h"

#include <iostream>
#include <string.h>

#include <stdlib.h>                               /* for memcpy */
#include <stdio.h>
#ifndef WIN32
#include <unistd.h>
#else
#include <io.h>
#endif
using namespace std;

#define NTRIES 10                                 /* # of tries before we abort the slab */
//#define SOLID_VOLUME 	/* define if volumes are considered one piece */

#define EXTRACT_SLAB                              /**/
#define DO_ONE_SLICE                              /**/

extern void adj(CTI_TETRA tetra[], int ntetra);
int read_contours(char *);

/* create a pt which is not on any contour; you may test that a pt was
 * created by this routine by checking if (nextpt[index] == -1);
 * RET: index of created pt */
static int cre8pt(float x,float y)
{
/* make sure there is room */
    if (numpts >= CTI_MAX_PTS)
    {
        fprintf(stderr,"  ERR(cre8pt): too many points IN ALL\n");
	// return; /**/
        exit(1);
    }

    savept[numpts].x = savept[numpts].xo = x;
    savept[numpts].y = savept[numpts].yo = y;
    savept[numpts].z = 99999;                     /* no z */
    savept[numpts].con = 99999;                   /* no con */

    slicenum[numpts] = 99999;

    nextpt[numpts] = -1;                          /* no next */
    prevpt[numpts] = -1;                          /* no prev */
    return numpts++;
}


/* even pts are across the frame from one another; same for odd. That
 * is, the 2 even pts are opposite corners of the bounding box -- same
 * for odd. */

static int cre8framepts(float bsize, float bminx, float bminy)
{
    int retframestart;                            /* first pt in the 8-pt frame */
    float xc,yc, side_len;

/* this is a hack which guarantees that all V-edges will be inside
 * the frame's V-edges (Bernhard figured this fudge factor
 * empirically).  */
    bminx -= 1.5*bsize;                           /* expand to 4 times normal size */
    bminy -= 1.5*bsize;
    bsize *= 4.0;

/* LATER: the order of these pts probably depends on cw */

    retframestart =
        cre8pt(bminx        , bminy);
    cre8pt(bminx        , bminy + bsize);
    cre8pt(bminx + bsize, bminy + bsize);
    cre8pt(bminx + bsize, bminy        );

    side_len = 1.4*bsize/2;                       /* approx of hypotenuse of a 45-90-45 triangle */

    xc = bminx + (bsize/2);                       /* center of the bounding box */
    yc = bminy + (bsize/2);

    cre8pt(xc           , yc - side_len);
    cre8pt(xc - side_len, yc           );
    cre8pt(xc           , yc + side_len);
    cre8pt(xc + side_len, yc           );

    return retframestart;
}


static int create_slab_frame_points(float bsize, float bminx, float bminy)
{
    int retframestart;                            /* first pt in the 8-pt frame */
    float xc,yc, side_len;

/* this is a hack which guarantees that all V-edges will be inside
 * the frame's V-edges (Bernhard figured this fudge factor
 * empirically).  */
    bminx -= 1.5*bsize;                           /* expand to 4 times normal size */
    bminy -= 1.5*bsize;
    bsize *= 4.0;

/* LATER: the order of these pts probably depends on cw */

    retframestart =
        cre8pt(bminx        , bminy);
    cre8pt(bminx        , bminy + bsize);
    cre8pt(bminx + bsize, bminy + bsize);
    cre8pt(bminx + bsize, bminy        );

    side_len = 1.4*bsize/2;                       /* approx of hypotenuse of a 45-90-45 triangle */

    xc = bminx + (bsize/2);                       /* center of the bounding box */
    yc = bminy + (bsize/2);

    cre8pt(xc           , yc - side_len);
    cre8pt(xc - side_len, yc           );
    cre8pt(xc           , yc + side_len);
    cre8pt(xc + side_len, yc           );

    return retframestart;
}


/* for debuggery */
static void print_tet()
{
    int i;

    if (cti_debug)                                /**/
        for (i=0; i < ntetra; i++)
            printf("tet %4d: [%s%s%s] %s %s pts %2d %2d %2d %2d  nei %2d %2d %2d %2d\n",
                i,
                (tetra[i].type & CTI_T1) ? "T1 " : "",
        (tetra[i].type & CTI_T2) ? "T2 " : "",
        (!(tetra[i].type & (CTI_T1|CTI_T2))) ? "T12" : "",
        (tetra[i].type & CTI_ELIM) ? "ELIM" : "    ",
        (tetra[i].type & CTI_MARK) ? "MARK" : "    ",
        tetra[i].p[0],
        tetra[i].p[1],
        tetra[i].p[2],
        tetra[i].p[3],
        tetra[i].t[0],
        tetra[i].t[1],
        tetra[i].t[2],
        tetra[i].t[3]);
}


/* slice numbers might not be consecutive. After findslabs() is
 * called, slabslice[i] will contain the lower slice number for slab
 * i; slabslice[CTI_MAX_SLICES + i] will contain the upper. Note that
 * slabslice[CTI_MAX_SLICES + numslabs - 1] is the highest slice
 */

static void findslabs()
{
    int slice;

    numslabs = 0;

#ifdef SOLID_VOLUME
    int savenumslabs;
    for (slice=0; slice < numslices; slice++)
    {
        if (num_contours[slice] != 0)
        {
            slabslice[numslabs++] = slice;
            numslabs++;
        }
    }

    numslabs--;                                   /* there are 1 fewer slabs than slices */
    savenumslabs = numslabs;
    numslabs--;                                   /* convert count to index */

/* count backwards to generate the upper slices */
    for (slice = numslices - 1; slice >= 0 && numslabs >= 0; slice--)
    {
        if (num_contours[slice] == 0)
            continue;
        slabslice[numslabs-- + CTI_MAX_SLICES] = slice;
    }

    numslabs = savenumslabs;

#else

/* scan all possible slices, looking for 2 IN A ROW which have contours */
    for (slice=0; slice < CTI_MAX_SLICES-1; slice++)
    {
/* see if there is anything in missing in 2 consecutive slices */
        if ( (num_contours[slice] == 0) ||  (num_contours[slice+1] == 0) )
            continue;
        slabslice[numslabs                 ] = slice;
        slabslice[numslabs + CTI_MAX_SLICES] = slice+1;
        numslabs++;
    }
#endif
}


static void init_globals()
{
    int pt1, pt2, tet, e;
    setbuf(stderr, NULL);
    setbuf(stdout, NULL);

    numpts = numslices = numslabs = ntetra = 0;
    con_factor = zdistance = 1.0;

    for (tet=0; tet < CTI_MAX_TETRAS; tet++)
    {
        tetra[tet].p[0] =
            tetra[tet].p[1] =
            tetra[tet].p[2] =
            tetra[tet].p[3] = -1;

        tetra[tet].t[0] =
            tetra[tet].t[1] =
            tetra[tet].t[2] =
            tetra[tet].t[3] = -1;

        tetra[tet].type = tetra[tet].slabnum = tetra[tet].searchnum = -1;
    }
    for (e=0; e < 6*CTI_MAX_PTS_PER_SLICE; e++)
    {
        vedge1[e].xs   =
            vedge1[e].ys   =
            vedge1[e].xe   =
            vedge1[e].ye   = 0;
        vedge1[e].npts =
            vedge1[e].s    =
            vedge1[e].e    =
            vedge1[e].ts   =
            vedge1[e].te   = 0;

        vedge2[e].xs   =
            vedge2[e].ys   =
            vedge2[e].xe   =
            vedge2[e].ye   = 0;
        vedge2[e].npts =
            vedge2[e].s    =
            vedge2[e].e    =
            vedge2[e].ts   =
            vedge2[e].te   = 0;
    }
    for (pt1=0; pt1 < CTI_MAX_PTS_PER_SLICE; pt1++)
        for (pt2=0; pt2 < CTI_MAX_PTS_PER_SLICE; pt2++)
            edge1[pt1][pt2] = edge2[pt1][pt2] = 0;
}


/* find the bounding box around the contours */
/* RET: length of a side, and position of lower left corner */
static void find_range(float *bsize_p, float *bminx_p, float *bminy_p)
{
    int con, pt, startpt;
    float bminx, bminy;
    float bmaxx, bmaxy;
    float brangex, brangey;
    int slice;
    int virgin = 1;                               /* 1: no contours seen yet */

    for (slice=0; slice < CTI_MAX_SLICES-1; slice++)
    {
        if (num_contours[slice] == 0)
            continue;

        if (virgin)
        {
            bminx = bmaxx = savept[begin_contour[slice][0]].x;
            bminy = bmaxy = savept[begin_contour[slice][0]].y;
            virgin = 0;
        }

        for(con=0; con < num_contours[slice]; con++)
        {
            pt = startpt = begin_contour[slice][con];
            do
            {
                if (bminx > savept[pt].x)         /* min x & y */
                    bminx = savept[pt].x;
                if (bminy > savept[pt].y)
                    bminy = savept[pt].y;

                if (bmaxx < savept[pt].x)         /* max x & y */
                    bmaxx = savept[pt].x;
                if (bmaxy < savept[pt].y)
                    bmaxy = savept[pt].y;

                pt = nextpt[pt];
            } while (pt != startpt);
        }

        brangex = bmaxx - bminx;
        brangey = bmaxy - bminy;
    }

    if (virgin)
        printf("  ERR: (find_range) no contours found\n");
    else
    {
/* copy to RET values */
        *bminx_p = bminx;
        *bminy_p = bminy;
        *bsize_p = MAX(brangex, brangey);         /* square bounding box */

        PRINTF("bsize=%f, bminx=%f, bminy=%f\n",
            *bsize_p, *bminx_p, *bminy_p);

    }
}


/* find the bounding box around the contours */
/* RET: length of a side, and position of lower left corner */
static void find_slab_range(int slab, float *bsize_p, float *bminx_p,
float *bminy_p)
{
    int i, con, pt, startpt;
    float bminx, bminy;
    float bmaxx, bmaxy;
    float brangex, brangey;
    int slice;
    int virgin = 1;                               /* 1: no contours seen yet */

    for (i=0; i<2; i++)
    {
        if (i==0) slice=slabslice[slab];
        else      slice=slabslice[slab+CTI_MAX_SLICES];

        if (num_contours[slice] == 0)
            continue;

        if (virgin)
        {
            bminx = bmaxx = savept[begin_contour[slice][0]].x;
            bminy = bmaxy = savept[begin_contour[slice][0]].y;
            virgin = 0;
        }

        for(con=0; con < num_contours[slice]; con++)
        {
            pt = startpt = begin_contour[slice][con];
            do
            {
                if (bminx > savept[pt].x)         /* min x & y */
                    bminx = savept[pt].x;
                if (bminy > savept[pt].y)
                    bminy = savept[pt].y;

                if (bmaxx < savept[pt].x)         /* max x & y */
                    bmaxx = savept[pt].x;
                if (bmaxy < savept[pt].y)
                    bmaxy = savept[pt].y;

                pt = nextpt[pt];
            } while (pt != startpt);
        }

        brangex = bmaxx - bminx;
        brangey = bmaxy - bminy;
    }

    if (virgin)
        printf("  ERR: (find_range) no contours found\n");
    else
    {
/* copy to RET values */
        *bminx_p = bminx;
        *bminy_p = bminy;
        *bsize_p = MAX(brangex, brangey);         /* square bounding box */

        PRINTF("bsize=%f, bminx=%f, bminy=%f\n",
            *bsize_p, *bminx_p, *bminy_p);

    }
}


/* * * * * * * *   c t i l e r   * * * * * * * */

/* Ret: 0 if err */
int ctiler(char *volume_data_id,                  /* fou data id to read volume */
char *tile_filename,                              /* filename to write tile set */
cti_flags_t flags,                                /* user specified parameters */
float the_perturb_dist                            /* (cm) max perturbation distance, if */
/* flags & CTI_PERT_DIST */
)
{
  // unused vars // int npts1, npts2, con;
    int con_index;
    int pt;
    FILE *fptile;
    int loslice, hislice, slab;
    // unused var // int slice;
    float bsize, bminx, bminy;
    // unused vars // float dx,dy, xi,yi;
    // unused var // int col;                       /* DEBUG: color */

/* record of multiple attempts */
    int retval;
    int ntries;
    int nslabfails = 0;
    int nslabaok   = 0;
    int nslabsaved = 0;

    int retries[150][NTRIES];
    int *ip;                                      /* int ptr */
    int attempt;
    int old_hislice;
    // unused var // char fail = 0;        /* failure getting data or files */
    char cw, allow_horiz, allow_caps, draw_normals;
    float perturb_dist = DEFAULT_PERTURB_DIST;

/*
 * options from caller
 */

    cw = flags & CTI_CW;
    allow_horiz = flags & CTI_ALLOW_HORIZ;
    if (flags & CTI_PERT_DIST)                    /* override default? */
        perturb_dist = the_perturb_dist;
    allow_caps = flags & CTI_CAPS;
    draw_normals = flags & CTI_DRAW_NORMALS;

    printf("Options: %scaps %scw p%f %snormals %s %s\n",
        allow_caps ? "" : "no",
        cw         ? "" : "c",
        perturb_dist,
        draw_normals ? "" : "no",
        volume_data_id,
        tile_filename);

/*
 * read contours for a single volume; open tile data file
 */

    ncon = 0;

    if (! read_contours(volume_data_id))          /* read into Phase 1 globals */
        return 0;

    if (! (fptile = fopen(tile_filename, "w")))   /* overwrite file contents */
    {
        printf("  ERR: tile file <%s> could not be opened for writing. Goodbye\n",
            tile_filename);
        return 0;
    }

/*
  Phase 1 - massage contours
*/
/* pre-process contours & pts:
 *   - eliminate duplicate pts in each contour
 *   - eliminate crossed contour loops
 *   - enforce contour direction (CW/CCW)
 */

    massage(cw);

/*
  Phase 2 - tile
*/

/* pre-process contours & pts for tiling:
 *   - convert con format from phase 1 format
 *   - eliminate colinear vertices
 *   - connect slices into slabs,
 *   - find the bounding box (range) of all pts in all cons
 *   - create an initial Del Tri frame, per Bernhard's
 *     instructions
 *   - init # of retries for each slab to -none- (-1)
 */

    init_globals();                               /* (yuk!) guarantee that phase 2 globals */
/* are initialized */
    cvt_phase();                                  /* read in-memory data from last */
/* phase */
    elim_coli_all();                              /* dump consecutive colinear points */
    findslabs();                                  /* group slices into slabs */
    for (pt=0; pt < numpts; pt++)
    {
        savept[pt].xo = savept[pt].x;             /* save unperturbed copies */
        savept[pt].yo = savept[pt].y;
    }
    //addrand_all(perturb_dist);/**/

#ifndef EXTRACT_SLAB
    find_range(&bsize, &bminx, &bminy);           /* how big are the contours? */
    turn45 = 0;
    framestart = cre8framepts(bsize, bminx, bminy);
/* frame the contour area */
#endif

/* init record of retries */
    unsigned int i;
    for (i=0, ip = &retries[0][0];
        i < sizeof(retries)/sizeof(retries[0][0]);
        i++, ip++)
    *ip = -1;

/*
 * Tile. For each slab (pair of slices):
 * - isolate this slab's points (points on loslice and hislice)
 * - [UNIMPLEMENTED] calc the buckets for each edge in the slab
 * - reattempt for each slab until sucessful or out of patience:
 *   - calc Del Tri, Vor diagram (add pts during Del Tri calc to
 *     guarantee that the contour pts are sampled finely enough)
 *   - create the tetrahedrilation (Merge step)
 *   - sculpt the tetrahedrilation
 *   - if successful:
 *     - convert tet's to tile-set of tet's surface
 *     - write the tiles (surfacetotile())
 *   - else:
 *     - change some data heuristically and attempt again...
 */

    slab = 0;
    ntries = 0;
    old_hislice = -1;
    printf("\n");

    while (slab < numslabs)
    {
#ifdef EXTRACT_SLAB
/* debuggery */
        find_slab_range(slab,&bsize, &bminx, &bminy);
        framestart = create_slab_frame_points(bsize, bminx, bminy);
        srand(1);
        turn45 = 0;
#endif

        loslice = slabslice[slab];
        hislice = slabslice[slab + CTI_MAX_SLICES];

        printf("Slab %d/%d ... ", slab, numslabs);
        PRINTF("slab %d/%d ( slice %d-%d, \"z1=%.1f z2=%.1f\" )\n",
            slab, numslabs, loslice, hislice,
            savept[begin_contour[loslice][0]].z,
            savept[begin_contour[hislice][0]].z);

/* assume caps & remove cap status if con causes tiles */
        for (con_index=0; con_index < MAS_NCONS; con_index++)
            caps[con_index] = 1;

#ifdef DO_ONE_SLICE
        if (loslice == old_hislice)
        {
            memcpy((char*)tri1,(char*)tri2, sizeof(CTI_TRI)*ntri2);
            ntri1 = ntri2;
        }
        else
        {
            retval = del(loslice, tri1,&ntri1, perturb_dist);
            if (retval != 0)
                goto attemptagain;                /* ugly, but simplifies code */
        }
        old_hislice = hislice;
#else
        retval = del(loslice, tri1,&ntri1, perturb_dist);
        if (retval != 0)
            goto attemptagain;                    /* ugly, but simplifies code */
#endif

        retval = del(hislice, tri2,&ntri2, perturb_dist);
        if (retval != 0)
            goto attemptagain;                    /* ugly, but simplifies code */

        vor(tri1,ntri1,char(cw));
        vor(tri2,ntri2,cw);

        if (ntries & 1)                           /* odd number of tries */
        {
            retval = mer(tri1,ntri1, loslice,
                tri2,ntri2, hislice,
                tetra, &ntetra,
                (float)(ntries+1)/NTRIES * perturb_dist, cw);
        }
        else                                      /* even number of tries and first merge */
        {
            retval = mer(tri2,ntri2, hislice,
                tri1,ntri1, loslice,
                tetra, &ntetra,
                (float)(ntries+1)/NTRIES * perturb_dist, cw);
        }
        if (ntries & 1)                           /* odd number of tries */

            attemptagain:
            cout << "i" << endl;
        if (ntries & 1)                           /* odd number of tries */
            switch (retval)
            {
                case 0:                           /* aok */
                    cout << "0" << endl;
/* drawvor(tri1,ntri1,1,1); drawdel(tri1,ntri1,2,1);*/
/* drawvor(tri2,ntri2,3,2); drawdel(tri2,ntri2,4,2);*/
/* drawvor(tri1,ntri1,1,3);*/
/* drawvor(tri2,ntri2,3,3);*/
/* attempt_fiddle();*/
                    PRINTF("  INFO:  SCULPTING %d tet's...", ntetra); fflush(stdout);

                adj(tetra, ntetra);

                if (sculpt())
                {
                    PRINTF("OK\n");
                    printf("COMPLETED\n\n");

                    surfacetotile(fptile, tetra,ntetra, allow_horiz,
                        allow_caps, draw_normals);/**/

/* make_tiles is used when debugging. Set TETFLAG to
 * 999 to generate the normal result. a value of 0
 * exports the tiles of sculpted tetras. See
 * make_tiles for other options.
 */
// make_tiles(fptile, tetra, ntetra, allow_horiz,0,0); /**/

/* draw caps in tile-less cons; force end-caps on
 * first & last slices
 */
                    if (allow_caps)
                    {
                        draw_caps(fptile, loslice, (slab == 0),
                            CTI_T2, tetra,ntetra, draw_normals);
                        draw_caps(fptile, hislice, (slab == (numslabs-1)),
                            CTI_T1, tetra,ntetra, draw_normals);
                    }

                    slab++;                       /* go to next slab */
                    if (ntries) nslabsaved++;
                    nslabaok++;
                    ntries = 0;
                    break;
            }
            else
            {
                printf("  ERR: sculpter reports a cycle (internal err)\n");
                retries[slab][ntries++] = 8;
/* no break! */
            }
            case 1:                               /* find_in_del - slice 1 */
            case 2:                               /* find_in_del - slice 2 */
            case 3:                               /* */
            case 4:                               /* 30 tries in this plane */
            case 5:                               /* other plane err */
            case 6:                               /* del perturb err */
            case 7:                               /* addrand_slice err */
            case 8:                               /* cycle found by sculpter */
            case 9:                               /* find_in_vor - edge reversed */
                cout << "9" << endl;
                if (ntries == NTRIES)
                {
                    abortion:
                    printf("ABORTED: slab %d (slices %d & %d)\n\n",
                        slab, loslice, hislice);
                    retries[slab][ntries] = retval;
                    nslabfails++;
                    slab++;                       /* abort slab */
                    ntries = 0;
                }
                else if (ntries == 5)
                {
                    retries[slab][ntries] = retval;
                    addval_slice(hislice, .7 * perturb_dist, -perturb_dist);
                    printf("ROLL BACK: offset pts in hi slice (attempt #%d of %d)\n",
                        ntries+1, NTRIES);
                    ntries++;
                }
                else if (ntries == 6)
                {
                    retries[slab][ntries] = retval;
                    addval_slice(loslice, -.7 * perturb_dist, perturb_dist);
                    printf("ROLL BACK: offset pts in lo slice (attempt #%d of %d)\n",
                        ntries+1, NTRIES);
                    ntries++;
                }
                else if (ntries == 7)
                {
                    retries[slab][ntries] = retval;
                    addval_slice(hislice, 2, .5 * perturb_dist);
                    printf("ROLL BACK: offset pts in hi slice (attempt #%d of %d)\n",
                        ntries+1, NTRIES);
                    ntries++;
                }
                else
                {
                    retries[slab][ntries] = retval;
                    ntries++;
                    printf("ROLL BACK: perturb pts in slab (attempt #%d of %d)\n",
                        ntries+1, NTRIES);

/* a mini reattempt loop which tries to mix up points a */
/* few times */
                    while (   (ntries < 5 || ntries > 7)
                        && (   ! addrand_slice(loslice,
                        (float)(ntries+1)/NTRIES * perturb_dist)
                        || ! addrand_slice(hislice,
                        (float)(ntries+1)/NTRIES * perturb_dist)))
                    {
                        retries[slab][ntries] = 7;/* record err */
                        ntries++;
                        printf("  ERR: perturb may corrupt.");
                        printf("(attempt #%d of %d) w/ perturb slack=%f\n",
                            ntries, NTRIES, (float)(ntries+1)/NTRIES * perturb_dist);

                        if (ntries >= NTRIES)
                            goto abortion;        /* last chance failed - induce abortion */
                    }
		    //turn45 = 1 - turn45;/**/
                }
                break;
            case -1:
                fprintf(stderr, "FATAL ERR. ABORT.\n");
                return 0;
        }

        print_tet();
    }
    cout << "l" << endl;

/* the only real measure of success is the status msg */
    printf("\nSummary of all slabs:\n");
    printf("%3d aok\n",                         nslabaok);
    printf("    %3d no problem\n",              nslabaok - nslabsaved);
    printf("    %3d saved by reattempt facility\n", nslabsaved);
    printf("%3d gave up\n",                     nslabfails);
    printf("%3d total\n\n",                     numslabs);

/* slab-by-slab summary */
    printf("\nRecap of problem slabs:\n");
    for (slab = 0; slab < numslabs; slab++)
    {
        for (attempt = 0; attempt < NTRIES; attempt++)
        {
            if (attempt == 0 && retries[slab][0] != -1)
                printf("slab %d: ", slab);
            if (retries[slab][attempt] != -1)
                printf("%d..", retries[slab][attempt]);
        }
        if (retries[slab][0] != -1)
            printf("\n");
    }
    cout << "m" << endl;

    fclose(fptile);
    return 1;                                     /* aok */
}
