extern int mer(CTI_TRI tria[],
	       int ntria,
	       int slicea,
	       CTI_TRI trib[],
	       int ntrib,
	       int sliceb,
	       CTI_TETRA tetra[],
	       int *ntetra_p,
	       float maxper,
	       char cw 	       /* TRUE: contour points are considered in */
	       		       /* clockwise order, else CCW */
	       );
