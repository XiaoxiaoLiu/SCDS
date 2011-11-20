/* rtpt.h  --- sample for C binding document --- J. Jacky 3-Mar-1992  */

#ifndef _RTPT_H
#define _RTPT_H

/*
revision 1.14:
gst@UNC: change blocks to blocks_p per Aug 92 TR-92-1 spec; untabify

revision 1.13:
gst@UNC: add #ifndef _RTPT_H to protect file from multiple inclusions

revision 1.12:
massive changes: Jon Jacky's version of this file is now in use.
Remove image-specific typedefs (ct_pixel_t, ...), enums are positive
for valid values, negative for errors (only listed once in unk_abs_e);
RTPT_ST_OUTOFMEM added; collimator dir enum's collapsed into a single
enum; arrays handled uniformly -- max size at site is array size, and
# of elements in use is always specified somewhere; several struct
members were set's -- these were changed to ptr's to sets; "name"
removed from grid values (for now -- consider adding it back in when a
tool actually requires it); pre-defined scanner names removed; lots
of spellings changed for consistency; lots of large_floats changed to
small_floats; all calls now return only status -- additional info
is now an updated parameter.

gst@UNC added gratis_id to rtpt_id_t for identifying gratis objs; -1
means "invalid" id.

gst@UNC added name to grid values (for UNC Fou code ONLY) -- it is not
set by tools, so the Fou will ignore requests to store grid values.

*/



/* All strings must be declared as rtpt_string_t */
typedef char *        rtpt_string_t;

/* Strings */

/* Special strings for attribute values in RTPT classes */ 
#define RTPT_UNKNOWN    "UNKNOWN"
#define RTPT_ABSENT     "ABSENT"
#define RTPT_DRP        "DRP"
#define RTPT_UNITS_HOUNSFIELD "hounsfield number"

/* Must assign this value to obj_id member of struct that represents a set */
#define RTPT_SET "set"

/* C binding of all class names from TR-91-1 */
#define RTPT_POLYLINE                         "polyline"
#define RTPT_CONTOUR                          "contour"
#define RTPT_SLAB                             "slab"
#define RTPT_VOLUME                           "volume"
#define RTPT_SOLID                            "solid"
#define RTPT_POINT                            "point"
#define RTPT_POINT_VALUE                      "point value"
#define RTPT_GRID_GEOM                        "grid geometry"
#define RTPT_WEDGE_DESCR                      "wedge description"
#define RTPT_BEAM_DELIV                       "beam delivery"
#define RTPT_S_BEAM_DELIV                     "symmetric beam delivery"
#define RTPT_C_BEAM_DELIV                     "combination beam delivery"
#define RTPT_VJ_BEAM_DELIV                    "variable jaw beam delivery"
#define RTPT_M_BEAM_DELIV                     "multileaf beam delivery"
#define RTPT_BLOCK                            "block"
#define RTPT_BEAM                             "beam"
#define RTPT_S_BEAM                           "symmetric beam"
#define RTPT_C_BEAM                           "combination beam"
#define RTPT_VJ_BEAM                          "variable jaw beam"
#define RTPT_M_BEAM                           "multileaf beam"
#define RTPT_GRID_VALUES                      "grid values"
#define RTPT_IMAGE                            "image"
#define RTPT_IMAGE_2D                         "image 2d"
#define RTPT_IMAGE_3D                         "image 3d"


/*
 * Scanner names are for FOU code ONLY (not for tools)
 */
#define RTPT_SCANNER_GE_9800     "GE9800"
#define RTPT_SCANNER_SOMATOM_DRX "SomatomDrx"


/* Enums */

/* Negative values so won't conflict with values of other enums */
enum unk_abs_e { RTPT_UNKNOWN_E = -1, RTPT_ABSENT_E = -2 };

/* Foundation operation status return values */
enum  rtpt_status_e
{
    RTPT_ST_OK,                 /* success */
    RTPT_ST_FAIL,               /* generic failure */
    RTPT_ST_OUTOFMEM            /* insufficient dynamic memory failure */
};

/* Beam modalities */
enum  beam_modality_e
{
    RTPT_BEAM_PHOTON,
    RTPT_BEAM_ELECTRON
};

/* Collimator motion direction for leaf collims and one-way-only var jaw */
enum collim_dir_e
{
    RTPT_X,
    RTPT_Y
};

/* Image modalities */
enum image_type_e
{
    RTPT_IMTYPE_CT,             /* Computed Tomography */
    RTPT_IMTYPE_SIM,            /* Rad Therapy Simulation Film */
    RTPT_IMTYPE_PORT,           /* Rad Therapy Port Film */
    RTPT_IMTYPE_DRR             /* Digitally Recomputed Radiograph */
}; 

/* Image pixel type */
enum pixel_type_e
{
    RTPT_PIXTYPE_SMALL_UINT,
    RTPT_PIXTYPE_MEDIUM_UINT,
    RTPT_PIXTYPE_MEDIUM_INT,
    RTPT_PIXTYPE_LARGE_INT
};

/* Numeric types.  Choose types to provide at least the specified range */

/*      compiler types     RTPT types      Approx low..Approx high */

typedef unsigned char      small_uint;    /*    0       255  */
typedef short int          medium_int;    /*  -32K       32K */
typedef unsigned short int medium_uint;   /*    0        64K */
typedef int                large_int;     /*   -2B        2B */

typedef float              small_float;   /* see body of document for range */
typedef double             large_float;

/* Special numbers */
#define RTPT_MAXN_LEAFPAIRS 25
#define RTPT_MAXN_WEDGEPOS  4

/* Contour vertices */
typedef struct 
{
   small_float  x;
   small_float  y;
} rtpt_vertex_t;

/* Transformation matrices */
typedef struct rtpt_transform
{
   large_float transform[4][4]; /* transformation matrix */
}  rtpt_transform_t;            /* [column:0..3][row:0..3] */

/* Identifiers - first member of every struct that is a set or RTPT class */
typedef struct rtpt_id
{
    rtpt_string_t type;         /* object's class type from #define's above */
    int gratis_id;              /* added to id obj, gratis-wise */
} rtpt_id_t;

/* Set */
typedef struct rtpt_set
{
    rtpt_id_t         obj_id;           /* has value RTPT_SET */
    struct rtpt_set  *next_p;           /* ptr to next node */
    void             *set_element_p;    /* element can be any type */
} rtpt_set_t; 


/* Foundation classes from TR-91-1.  

Where the C struct member name is spelled  exactly as in TR-91-1, no comment
is necessary.  When the C struct member name is not spelled exactly as in
TR-91-1, the TR-91-1 spelling  appears in the comment.  Additional comments
explain C binding machinery. */

/* Polyline */
typedef struct rtpt_polyline
{
    rtpt_id_t         obj_id;    /* has value RTPT_POLYLINE */
    rtpt_vertex_t    *vertices_p;/* vertices -- array of rtpt_vertex_t */
    medium_uint       nvertices; /* n of elements in vertices */
    rtpt_transform_t  trans;     /* transform */  
} rtpt_polyline_t;

/* Contour */
typedef struct rtpt_contour
{
    rtpt_id_t       obj_id;     /* has value RTPT_CONTOUR */
    rtpt_polyline_t poly;       /* ancestry of contour */
} rtpt_contour_t;


/* Slab */
typedef struct rtpt_slab
{
    rtpt_id_t        obj_id;    /* has value RTPT_SLAB */
    rtpt_contour_t   contour;   /* ancestry of slab */ 
    small_float      z_plus;    /* z plus */ 
    small_float      z_minus;   /* z minus */ 
} rtpt_slab_t;


/* Volume */
typedef struct rtpt_volume
{
    rtpt_id_t         obj_id;           /* has value RTPT_VOLUME */
    rtpt_string_t     name;             /* name of volume */
    rtpt_set_t       *contours_p;       /* contours - set of rtpt_contour_t */
    small_float       z_plus;           /* z plus */ 
    small_float       z_minus;          /* z minus */ 
} rtpt_volume_t;

/* Solid */
typedef struct rtpt_solid
{
    rtpt_id_t       obj_id;     /* has value RTPT_SOLID */
    rtpt_volume_t   volume;     /* ancestry of solid */
    small_float     density;    
} rtpt_solid_t;

/* Point */
typedef struct rtpt_point
{
    rtpt_id_t         obj_id;           /* has value RTPT_POINT */
    rtpt_string_t     name;             
    small_float       x;                /* x coord */
    small_float       y;                /* y coord */
    small_float       z;                /* z coord */
} rtpt_point_t;

/* Point With Value */
typedef struct rtpt_point_value
{
    rtpt_id_t         obj_id;       /* has value RTPT_POINT_VALUE */
    rtpt_point_t      point;        /* ancestry */
    small_float       value;        
} rtpt_point_value_t;

/* Grid Geometry */
typedef struct rtpt_grid_geom
{
    rtpt_id_t            obj_id;        /* has value RTPT_GRID_GEOM */
    rtpt_string_t        name;          
    small_float          x_size;        /* x size */ 
    small_float          y_size;        /* y size */ 
    small_float          z_size;        /* z size */ 
    medium_uint          n_x;           /* n x samples */
    medium_uint          n_y;           /* n y samples */
    medium_uint          n_z;           /* n z samples */
    rtpt_transform_t     trans;         /* transform */
} rtpt_grid_geom_t;

/* Wedge Description */
typedef struct rtpt_wedge_desc
{
    rtpt_id_t       obj_id;             /* has value RTPT_WEDGE_DESCR  */
    rtpt_string_t   name;               
    small_float     angle;              
    small_float     retraction;         
    small_float     orientations[RTPT_MAXN_WEDGEPOS]; 
    small_uint      num_orient;         /* n of elements in orientations */
} rtpt_wedge_desc_t;

/* Beam Delivery System */
typedef struct rtpt_beam_deliv
{
    rtpt_id_t             obj_id;            /* has value RTPT_BEAM_DELIV*/
    rtpt_string_t         unit_name;         /* unit name */ 
    rtpt_string_t         mfgr_model;        /* mfgr model */ 
    enum beam_modality_e  modality;          
    small_float           energy;            
    small_float           source_size;       /* source size */
    small_float           source_axis;       /* source to axis distance */   
    rtpt_set_t           *wedges_p;          /* wedges --- 
                                                 set of rtpt_wedge_desc_t */ 
} rtpt_beam_delivery_t;

/* Symmetric Beam Delivery System */
typedef struct rtpt_s_beam_deliv
{
    rtpt_id_t               obj_id;     /* has value RTPT_S_BEAM_DELIV*/
    rtpt_beam_delivery_t    beam_deliv; /* ancestry of symmetric beam deliv. */
    small_float             coll_x_lim; /* collimator x limit */
    small_float             coll_y_lim; /* collimator y limit */
} rtpt_s_beam_delivery_t;

/* Combination Beam Delivery System */
typedef struct rtpt_c_beam_deliv
{
    rtpt_id_t              obj_id;      /* has value RTPT_C_BEAM_DELIV*/
    rtpt_beam_delivery_t   beam_deliv;  /* ancestry of combin. beam deliv.*/
    enum collim_dir_e      direction;   /* asymmetry direction */
    small_float            coll_retraction;  /* collimator retraction */
    small_float            coll_extension; /* collimator extension */
    small_float            coll_limit;  /* collimator limit */ 
} rtpt_c_beam_delivery_t;

/* Variable Jaw Beam Delivery System */
typedef struct rtpt_vj_beam_deliv
{
    rtpt_id_t               obj_id;     /* has value RTPT_VJ_BEAM_DELIV*/
    rtpt_beam_delivery_t    beam_deliv;/* ancestry of vari. jaw beam deliv. */
    small_float             coll_retraction;/* collimator retraction  */
    small_float             coll_extension; /* collimator extension   */
} rtpt_vj_beam_delivery_t;

/* Multileaf Beam Delivery System */
typedef struct rtpt_m_beam_deliv
{
    rtpt_id_t               obj_id;     /* has value RTPT_M_BEAM_DELIV*/
    rtpt_beam_delivery_t    beam_deliv; /* ancestry of multileaf beam deliv.*/
    enum collim_dir_e       direction;  /* leaf direction */ 
    small_uint              num_leaf;   /* n leaf pairs */
    small_float             retraction; /* leaf retraction */
    small_float             extension;  /* leaf extension */
    small_float             thickness;  /* leaf thickness */
    small_float             positions[RTPT_MAXN_LEAFPAIRS];/* leaf positions */
    small_uint              num_positions; /* n of elements in positions */
} rtpt_m_beam_delivery_t;

/* Block */
typedef struct rtpt_block
{
    rtpt_id_t         obj_id;           /* has value RTPT_BLOCK*/
    rtpt_contour_t    contour;          /* ancestry of block */
} rtpt_block_t;

/* Beam */
typedef struct rtpt_beam
{
    rtpt_id_t          obj_id;         /* has value RTPT_BEAM*/
    rtpt_string_t      name;           
    rtpt_string_t      delivery;       /* delivery system */
    small_float        units;          /* monitor units */
    small_float        lateral;        /* couch lateral */
    small_float        longit;         /* couch longitudinal */
    small_float        height;         /* couch height */
    small_float        c_angle;        /* couch angle */
    small_float        g_angle;        /* gantry angle */
    small_float        coll_angle;     /* collimator angle */
    rtpt_string_t      wedge;          
    small_float        orientation;    /* wedge orientation */
    rtpt_set_t        *blocks_p;       /* set of rtpt_block_t */
} rtpt_beam_t;

/* Symmetric Beam */
typedef struct rtpt_s_beam
{
    rtpt_id_t     obj_id;          /* has value RTPT_S_BEAM */
    rtpt_beam_t   beam;            /* ancestry of symmetric beam */
    small_float   coll_x;          /* collimator x */
    small_float   coll_y;          /* collimator y */
} rtpt_s_beam_t;

/* Combination Beam */
typedef struct rtpt_c_beam
{
    rtpt_id_t     obj_id;       /* has value RTPT_C_BEAM */
    rtpt_beam_t   beam;         /* ancestry of combination beam */
    small_float   coll_sup;     /* collimator sup */ 
    small_float   coll_inf;     /* collimator inf */ 
    small_float   coll_width;   /* collimator */
} rtpt_c_beam_t;

/* Variable Jaw Beam */
typedef struct rtpt_vj_beam
{
    rtpt_id_t     obj_id;       /* has value RTPT_VJ_BEAM*/
    rtpt_beam_t   beam;         /* ancestry of variable jaw beam */
    small_float   coll_supx;    /* collimator x sup */ 
    small_float   coll_infx;    /* collimator x inf */
    small_float   coll_supy;    /* collimator y sup */ 
    small_float   coll_infy;    /* collimator y inf */ 
} rtpt_vj_beam_t;

/* Multileaf Beam */
typedef struct rtpt_m_beam
{
    rtpt_id_t     obj_id;       /* has value RTPT_M_BEAM*/
    rtpt_beam_t   beam;         /* ancestry of multileaf beam */
    small_float   leaf_sups[RTPT_MAXN_LEAFPAIRS]; /* leaf sup settings */
    small_float   leaf_infs[RTPT_MAXN_LEAFPAIRS]; /* leaf inf settings */
                                /* n of leaves actually used is in */
                                /* corresponding rtpt_m_beam_deliv_t */
} rtpt_m_beam_t;

/* Grid Values */
typedef struct rtpt_grid_values
{
    rtpt_id_t          obj_id;         /* has value RTPT_GRID_VALUE */
    rtpt_string_t      grid;           
    rtpt_string_t      name;           /* for UNC fou code only (not tools) */
    small_float      * values;         /* 3-d matrix:
                                        * [z:0..nzsamples-1]
                                        * [y:0..nysamples-1]
                                        * [x:0..nxsamples-1]
                                        * Order is plane, raster, pix
                                        */
} rtpt_grid_values_t;

/* Image */
typedef struct rtpt_image
{
    rtpt_id_t            obj_id;        /* has value RTPT_IMAGE*/
    rtpt_string_t        id;            
    rtpt_string_t        descrip;       /* description */
    rtpt_string_t        date;          /* acquisition date */
    rtpt_string_t        time;          /* acquisition time */
    enum image_type_e    type;          /* image type */ 
    rtpt_grid_geom_t     grid;          
    rtpt_string_t        scanner;       /* scanner type */ 
    large_int            range_low;     /* range low */ 
    large_int            range_high;    /* range high */ 
                                        /* cast range_low and range_high to */
                                        /* integer type of pixel before use */
    rtpt_string_t        units;         
    enum pixel_type_e    pixtype;       /* indicates integer type of pixel */
} rtpt_image_t;

/* Image 2D */
typedef struct rtpt_image_2d
{
    rtpt_id_t      obj_id;              /* has value RTPT_IMAGE_2D*/
    rtpt_image_t   image;               /* ancestry of 2d image */
    small_float    thickness;           
    void         * pixels;              /* 2d array of image; cast to 
                                         * integer type of pixel before use;
                                         * [y:0..nysamples-1]
                                         * [x:0..nxsamples-1]
                                         * Order is raster, pix
                                         */
} rtpt_image_2d_t;

/* Image 3D */
typedef struct rtpt_image_3d
{
    rtpt_id_t        obj_id;            /* has value RTPT_IMAGE_3D*/
    rtpt_image_t     image;             /* ancestry of 3d image */
    void           * voxels;            /* 3d array of image; cast to 
                                         * integer type of pixel before use;
                                         * [z:0..nzsamples-1]
                                         * [y:0..nysamples-1]
                                         * [x:0..nxsamples-1]
                                         * Order is plane, raster, pix
                                         */
} rtpt_image_3d_t;

/* Foundation operations */

/* Fetch All */
enum rtpt_status_e fou_fetchall(        /* output: status */
    rtpt_string_t      data_id,         /* input: data set identifer */
    rtpt_set_t       * classes_p,        /* input: classes */
    rtpt_set_t       * objects_p        /* output: objects */
    );

/* Store All */
enum rtpt_status_e fou_storeall(        /* output: status */
    rtpt_string_t      data_id,         /* input: data set identifier */
    rtpt_set_t       * objects_p        /* input: objects */
    );

/* Compute Dose Matrix */
enum rtpt_status_e fou_compute_dose_matrix( /* output: status */
        rtpt_set_t       * ana_set_p,    /* input: anatomy-solids */
        rtpt_image_3d_t  * im3d_p,       /* input: anatomy-image */
        void             * field_p,      /* input: field */
        rtpt_point_t     * drp_p,        /* input: delivery reference point */
        rtpt_grid_geom_t * dose_geom_p,  /* input: dose grid */
        rtpt_grid_values_t  * grid_values /* output: dose rate matrix */
     );

/* Compute Dose Points */
enum rtpt_status_e  fou_compute_dose_points( /* output: status */
        rtpt_set_t       * ana_set_p,    /* input: anatomy-solids */
        rtpt_image_3d_t  * im3d_p,       /* input: anatomy-image */
        void             * field_p,      /* input: field */
        rtpt_point_t     * drp_p,        /* input: delivery reference point */
        rtpt_set_t       * pt_values_p   /* updated: point dose rates */
        );


#endif
