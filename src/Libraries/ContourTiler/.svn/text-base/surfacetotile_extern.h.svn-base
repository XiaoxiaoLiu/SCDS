/* surfacetotile_extern.h - prototypes for surfacetotile.c */

/* Author: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

extern void surfacetotile(FILE *fptile,
			  CTI_TETRA tetra[],
			  int ntetra,
			  int allow_horiz, /* 1; allow horizontal tiles on surface */
			  int collect_caps,/* 1: update cap status */
			  int draw_normals);

extern int draw_tile(FILE *fptile,
		     CTI_TETRA *tetra_p,
		     int not_included, /* [0..3] the vertex who is NOT in this tile */
		     int allow_horiz,
		     int flip_normal,
		     int draw_normals
		     );

