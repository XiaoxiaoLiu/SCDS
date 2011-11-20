/* cti_tetra.h  - CTIler tetrahedra definitions for Boissonnat's algor.*/

/* Later, make allocation of tetra's and faces a dynamic process */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

/* History:
   22-Jul-91 gst  Combine tetratype & solidtype into type (cf Bernhard's code)
   17-DEC-89 gst  added tetratype, solidtype & assoc'd defines
   13-APR-90 gst  added slabnum
    1-OCT-90 gst  added CTI_ to all capitalized names, to avoid
                  name conflicts.
*/
 
#define CTI_MAX_TETRAS     (CTI_MAX_PTS)  /* fudge * (total # of pts) */

/* these values encode the type field of CTI_TETRA, which identifies a
 * tet as:
 *   T1, T2, or T12,
 *   eliminated or not eliminated,
 *   marked or not marked (as seen).
 * These bits are meant to be individually tested: */
#define CTI_T12    0x10		/* this tet is a T12 tet */
#define CTI_T1     0x1		/* this tet is a T1 tet (mutually */
				/* exclusive of T2) */ 
#define CTI_T2     0x2		/* this tet is a T2 tet (mutually */
				/* exclusive of T1) */
#define CTI_ELIM   0x4		/* this tet has been eliminated */
#define CTI_MARK   0x8		/* this tet has been marked as seen */

typedef struct tetra_t
{
    int p[4];			/* indices of four vertices in increasing order */
    int t[4];			/* index of tetra opposite each vertex */
    int type;			/* T1, T2, T12, elim, mark */
    int slabnum;
    int searchnum;		/* last search id that found this tet */
} CTI_TETRA;

