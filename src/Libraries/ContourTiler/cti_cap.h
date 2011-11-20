/* cti_cap.h - a set of tiles in the plane of a slice is a cap,
 *   which can be allowed or disallowed, per contour.
 */

/* Authors: Gregg Tracton & Jun Chen, UNC Dept of Radiation Oncology */

/* these values are meant to be OR'd */
#define     CTI_NOCAP 0x0
#define CTI_BOTTOMCAP 0x1
#define    CTI_TOPCAP 0x2
#define   CTI_BOTHCAP (CTI_BOTTOMCAP | CTI_TOPCAP)

