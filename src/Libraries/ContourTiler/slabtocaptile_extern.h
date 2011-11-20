/* slabtocaptile_extern.h - prototypes for slabtocaptile.c */

/* Author: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

extern void draw_caps(FILE *fptile,
		      int slice,
		      int all_cons, /* 1: cap all cons in slice */
				    /* 0: cap cons having (caps[con]==1) */
		      int tettype,  /* get cap from a tet that points in this dir */
		      CTI_TETRA tetra[],
		      int ntetra,
		      int draw_normals);
