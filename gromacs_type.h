!!gromacs/src/gromacs/legacyheaders/types

!!simple.h

#define XX      0           /* Defines for indexing in	*/
#define YY      1           /* vectors			*/
#define ZZ      2
#define DIM     3           /* Dimension of vectors		*/
#define XXXX    0           /* defines to index matrices */
#define XXYY    1
#define XXZZ    2
#define YYXX    3
#define YYYY    4
#define YYZZ    5
#define ZZXX    6
#define ZZYY    7
#define ZZZZ    8

typedef int gmx_bool;

typedef int         atom_id;      /* To indicate an atoms id         */
#define NO_ATID     (atom_id)(~0) /* Use this to indicate invalid atid */

typedef real            rvec[DIM];
typedef double          dvec[DIM];
typedef real            matrix[DIM][DIM];
typedef real            tensor[DIM][DIM];
typedef int             ivec[DIM];
typedef int             imatrix[DIM][DIM];

typedef long long int gmx_large_int_t;  //64bit

// others for compiler and check later

!! atoms.h
#inlude "simple.h"

enum {
    eptAtom, eptNucleus, eptShell, eptBond, eptVSite, eptNR
};
/* The particle type */

typedef struct {
    real           m, q;        /* Mass and charge                      */
    real           mB, qB;      /* Mass and charge for Free Energy calc */
    unsigned short type;        /* Atom type                            */
    unsigned short typeB;       /* Atom type for Free Energy calc       */
    int            ptype;       /* Particle type                        */
    int            resind;      /* Index into resinfo (in t_atoms)      */
    int            atomnumber;  /* Atomic Number or NOTSET              */
    char           elem[4];     /* Element name                         */
} t_atom;

typedef struct {
    char          **name;       /* Pointer to the residue name          */
    int             nr;         /* Residue number                       */
    unsigned char   ic;         /* Code for insertion of residues       */
    int             chainnum;   /* Iincremented at TER or new chain id  */
    char            chainid;    /* Chain identifier written/read to pdb */
    char          **rtp;        /* rtp building block name (optional)   */
} t_resinfo;

typedef struct {
    int      type;              /* PDB record name                      */
    int      atomnr;            /* PDB atom number                      */
    char     altloc;            /* Alternate location indicator         */
    char     atomnm[6];         /* True atom name including leading spaces */
    real     occup;             /* Occupancy                            */
    real     bfac;              /* B-factor                             */
    gmx_bool bAnisotropic;      /* (an)isotropic switch                 */
    int      uij[6];            /* Anisotropic B-factor                 */
} t_pdbinfo;

typedef struct {
    int   nr;                   /* Number of different groups           */
    int  *nm_ind;               /* Index in the group names             */
} t_grps;

typedef struct {
    int            nr;          /* Nr of atoms                          */
    t_atom        *atom;        /* Array of atoms (dim: nr)             */
                                /* The following entries will not       */
                                /* always be used (nres==0)             */
    char          ***atomname;  /* Array of pointers to atom name       */
                                /* use: (*(atomname[i]))                */
    char          ***atomtype;  /* Array of pointers to atom types      */
                                /* use: (*(atomtype[i]))                */
    char          ***atomtypeB; /* Array of pointers to B atom types    */
                                /* use: (*(atomtypeB[i]))               */
    int              nres;      /* The number of resinfo entries        */
    t_resinfo       *resinfo;   /* Array of residue names and numbers   */
    t_pdbinfo       *pdbinfo;   /* PDB Information, such as aniso. Bfac */
} t_atoms;

typedef struct {
    int           nr;           /* number of atomtypes                          */
    real         *radius;       /* GBSA radius for each atomtype                */
    real         *vol;          /* GBSA efective volume for each atomtype       */
    real         *surftens;     /* implicit solvent surftens for each atomtype  */
    real         *gb_radius;    /* GB radius for each atom type                 */
    real         *S_hct;        /* Overlap factors for HCT/OBC GB models        */
    int          *atomnumber;   /* Atomic number, used for QM/MM                */
} t_atomtypes;

!!block.h
#include "../../utility/gmxmpi.h"
#inlcude "idef.h"

typedef struct {
    int      nr;           /* The number of blocks			*/
    atom_id *index;        /* Array of indices (dim: nr+1)  */
    int      nalloc_index; /* The allocation size for index        */
} t_block;

typedef struct {
    int      nr;    /* The number of blocks			*/
    atom_id *index; /* Array of indices in a (dim: nr+1)	*/
    int      nra;   /* The number of atoms          */
    atom_id *a;     /* Array of atom numbers in each group  */
                    /* (dim: nra)				*/
                    /* Block i (0<=i<nr) runs from		*/
                    /* index[i] to index[i+1]-1. There will */
                    /* allways be an extra entry in index	*/
                    /* to terminate the table		*/
    int nalloc_index;           /* The allocation size for index        */
    int nalloc_a;               /* The allocation size for a            */
} t_blocka;

!!commrec.h

#define DD_MAXZONE  8
#define DD_MAXIZONE 4

typedef struct gmx_domdec_master *gmx_domdec_master_p_t;

typedef struct {
    int  j0;     /* j-zone start               */
    int  j1;     /* j-zone end                 */
    int  cg1;    /* i-charge-group end         */
    int  jcg0;   /* j-charge-group start       */
    int  jcg1;   /* j-charge-group end         */
    ivec shift0; /* Minimum shifts to consider */
    ivec shift1; /* Maximum shifts to consider */
} gmx_domdec_ns_ranges_t;

typedef struct {
    rvec x0;     /* Zone lower corner in triclinic coordinates         */
    rvec x1;     /* Zone upper corner in triclinic coordinates         */
    rvec bb_x0;  /* Zone bounding box lower corner in Cartesian coords */
    rvec bb_x1;  /* Zone bounding box upper corner in Cartesian coords */
} gmx_domdec_zone_size_t;

typedef struct {
    /* The number of zones including the home zone */
    int                    n;
    /* The shift of the zones with respect to the home zone */
    ivec                   shift[DD_MAXZONE];
    /* The charge group boundaries for the zones */
    int                    cg_range[DD_MAXZONE+1];
    /* The number of neighbor search zones with i-particles */
    int                    nizone;
    /* The neighbor search charge group ranges for each i-zone */
    gmx_domdec_ns_ranges_t izone[DD_MAXIZONE];
    /* Boundaries of the zones */
    gmx_domdec_zone_size_t size[DD_MAXZONE];
    /* The cg density of the home zone */
    real                   dens_zone0;
} gmx_domdec_zones_t;

typedef struct gmx_ga2la *gmx_ga2la_t;

typedef struct gmx_hash *gmx_hash_t;

typedef struct gmx_reverse_top *gmx_reverse_top_p_t;

typedef struct gmx_domdec_constraints *gmx_domdec_constraints_p_t;

typedef struct gmx_domdec_specat_comm *gmx_domdec_specat_comm_p_t;

typedef struct gmx_domdec_comm *gmx_domdec_comm_p_t;

typedef struct gmx_pme_comm_n_box *gmx_pme_comm_n_box_p_t;

typedef struct {
    int  npbcdim;
    int  nboundeddim;
    rvec box0;
    rvec box_size;
    /* Tells if the box is skewed for each of the three cartesian directions */
    ivec tric_dir;
    rvec skew_fac;
    /* Orthogonal vectors for triclinic cells, Cartesian index */
    rvec v[DIM][DIM];
    /* Normal vectors for the cells walls */
    rvec normal[DIM];
} gmx_ddbox_t;

typedef struct {
    /* these buffers are used as destination buffers if MPI_IN_PLACE isn't
       supported.*/
    int             *ibuf; /* for ints */
    int              ibuf_alloc;

    gmx_large_int_t *libuf;
    int              libuf_alloc;

    float           *fbuf; /* for floats */
    int              fbuf_alloc;

    double          *dbuf; /* for doubles */
    int              dbuf_alloc;
} mpi_in_place_buf_t;

typedef struct {
    /* The DD particle-particle nodes only */
    /* The communication setup within the communicator all
     * defined in dd->comm in domdec.c
     */
    int                    nnodes;
    MPI_Comm               mpi_comm_all;
    /* Use MPI_Sendrecv communication instead of non-blocking calls */
    gmx_bool               bSendRecv2;
    /* The local DD cell index and rank */
    ivec                   ci;
    int                    rank;
    ivec                   master_ci;
    int                    masterrank;
    /* Communication with the PME only nodes */
    int                    pme_nodeid;
    gmx_bool               pme_receive_vir_ener;
    gmx_pme_comm_n_box_p_t cnb;
    int                    nreq_pme;
    MPI_Request            req_pme[4];


    /* The communication setup, identical for each cell, cartesian index */
    ivec     nc;
    int      ndim;
    ivec     dim; /* indexed by 0 to ndim */
    gmx_bool bGridJump;

    /* PBC from dim 0 to npbcdim */
    int npbcdim;

    /* Screw PBC? */
    gmx_bool bScrewPBC;

    /* Forward and backward neighboring cells, indexed by 0 to ndim */
    int  neighbor[DIM][2];

    /* Only available on the master node */
    gmx_domdec_master_p_t ma;

    /* Are there inter charge group constraints */
    gmx_bool bInterCGcons;
    gmx_bool bInterCGsettles;

    /* Global atom number to interaction list */
    gmx_reverse_top_p_t reverse_top;
    int                 nbonded_global;
    int                 nbonded_local;

    /* The number of inter charge-group exclusions */
    int  n_intercg_excl;

    /* Vsite stuff */
    gmx_hash_t                 ga2la_vsite;
    gmx_domdec_specat_comm_p_t vsite_comm;

    /* Constraint stuff */
    gmx_domdec_constraints_p_t constraints;
    gmx_domdec_specat_comm_p_t constraint_comm;

    /* The local to gobal charge group index and local cg to local atom index */
    int   ncg_home;
    int   ncg_tot;
    int  *index_gl;
    int  *cgindex;
    int   cg_nalloc;
    /* Local atom to local cg index, only for special cases */
    int  *la2lc;
    int   la2lc_nalloc;

    /* The number of home atoms */
    int   nat_home;
    /* The total number of atoms: home and received zones */
    int   nat_tot;
    /* Index from the local atoms to the global atoms */
    int  *gatindex;
    int   gatindex_nalloc;

    /* Global atom number to local atom number list */
    gmx_ga2la_t ga2la;

    /* Communication stuff */
    gmx_domdec_comm_p_t comm;

    /* The partioning count, to keep track of the state */
    gmx_large_int_t ddp_count;


    /* gmx_pme_recv_f buffer */
    int   pme_recv_f_alloc;
    rvec *pme_recv_f_buf;

} gmx_domdec_t;

typedef struct gmx_partdec *gmx_partdec_p_t;

typedef struct {
    int       nsim;
    int       sim;
    MPI_Group mpi_group_masters;
    MPI_Comm  mpi_comm_masters;
    /* these buffers are used as destination buffers if MPI_IN_PLACE isn't
       supported.*/
    mpi_in_place_buf_t *mpb;
} gmx_multisim_t;

#define DUTY_PP  (1<<0)
#define DUTY_PME (1<<1)

typedef struct {
    int      bUse;
    MPI_Comm comm_intra;
    int      rank_intra;
    MPI_Comm comm_inter;

} gmx_nodecomm_t;

typedef struct {
    /* The nodeids in one sim are numbered sequentially from 0.
     * All communication within some simulation should happen
     * in mpi_comm_mysim, or its subset mpi_comm_mygroup.
     */
    int sim_nodeid, nnodes, npmenodes;

    /* thread numbers: */
    /* Not used yet: int threadid, nthreads; */
    /* The nodeid in the PP/PME, PP or PME group */
    int      nodeid;
    MPI_Comm mpi_comm_mysim;
    MPI_Comm mpi_comm_mygroup;

    /* MPI ranks within a physical node for hardware access */
    int            nrank_intranode;    /* nr of ranks on this physical node */
    int            rank_intranode;     /* our rank on this physical node */
    int            nrank_pp_intranode; /* as nrank_intranode, for particle-particle only */
    int            rank_pp_intranode;  /* as rank_intranode, for particle-particle only */

    gmx_nodecomm_t nc;

    /* For domain decomposition */
    gmx_domdec_t *dd;

    /* For particle decomposition */
    gmx_partdec_p_t pd;

    /* The duties of this node, see the defines above */
    int             duty;

    gmx_multisim_t *ms;

    /* these buffers are used as destination buffers if MPI_IN_PLACE isn't
       supported.*/
    mpi_in_place_buf_t *mpb;
} t_commrec;

#define MASTERNODE(cr)     (((cr)->nodeid == 0) || !PAR(cr))
/* #define MASTERTHREAD(cr)   ((cr)->threadid == 0) */
/* #define MASTER(cr)         (MASTERNODE(cr) && MASTERTHREAD(cr)) */
#define MASTER(cr)         MASTERNODE(cr)
#define SIMMASTER(cr)      ((MASTER(cr) && ((cr)->duty & DUTY_PP)) || !PAR(cr))
#define NODEPAR(cr)        ((cr)->nnodes > 1)
/* #define THREADPAR(cr)      ((cr)->nthreads > 1) */
/* #define PAR(cr)            (NODEPAR(cr) || THREADPAR(cr)) */
#define PAR(cr)            NODEPAR(cr)
#define RANK(cr, nodeid)    (nodeid)
#define MASTERRANK(cr)     (0)

#define DOMAINDECOMP(cr)   (((cr)->dd != NULL) && PAR(cr))
#define DDMASTER(dd)       ((dd)->rank == (dd)->masterrank)

#define PARTDECOMP(cr)     ((cr)->pd != NULL)

#define MULTISIM(cr)       ((cr)->ms)
#define MSRANK(ms, nodeid)  (nodeid)
#define MASTERSIM(ms)      ((ms)->sim == 0)

/* The master of all (the node that prints the remaining run time etc.) */
#define MULTIMASTER(cr)    (SIMMASTER(cr) && (!MULTISIM(cr) || MASTERSIM((cr)->ms)))

!!constr.h    //??relay on what file??

/* Abstract type for LINCS that is defined only in the file that uses it */
typedef struct gmx_lincsdata *gmx_lincsdata_t;

/* Abstract type for SHAKE that is defined only in the file that uses it */
typedef struct gmx_shakedata *gmx_shakedata_t;

/* Abstract type for SETTLE that is defined only in the file that uses it */
typedef struct gmx_settledata *gmx_settledata_t;

/* Abstract type for constraints */
typedef struct gmx_constr *gmx_constr_t;

/* Abstract type for essential dynamics that is defined only in edsam.c */
typedef struct gmx_edsam *gmx_edsam_t;

!!energy.h
#include "simple.h"

typedef struct {
    real   e;    /* The current energy.					        */
    double eav;  /* The running average                 */
    double esum; /* The sum of energies until now.			*/
} t_energy;

!!enums.h

~key types: eintmod; gmx_table_interation; gmx_table_format; gmx_nblist_kernel_geometry; gmx_nbkernel_elec; gmx_nbkernel_vdw;gmx_nblist_interation_type;

/* note: these enums should correspond to the names in gmxlib/names.c */

enum {
    epbcXYZ, epbcNONE, epbcXY, epbcSCREW, epbcNR
};

enum {
    etcNO, etcBERENDSEN, etcNOSEHOOVER, etcYES, etcANDERSEN, etcANDERSENMASSIVE, etcVRESCALE, etcNR
}; /* yes is an alias for berendsen */

#define ETC_ANDERSEN(e) (((e) == etcANDERSENMASSIVE) || ((e) == etcANDERSEN))

enum {
    epcNO, epcBERENDSEN, epcPARRINELLORAHMAN, epcISOTROPIC, epcMTTK, epcNR
}; /* isotropic is an alias for berendsen */

/* trotter decomposition extended variable parts */
enum {
    etrtNONE, etrtNHC, etrtBAROV, etrtBARONHC, etrtNHC2, etrtBAROV2, etrtBARONHC2,
    etrtVELOCITY1, etrtVELOCITY2, etrtPOSITION, etrtSKIPALL, etrtNR
};

/* sequenced parts of the trotter decomposition */
enum {
    ettTSEQ0,  ettTSEQ1,  ettTSEQ2,  ettTSEQ3,  ettTSEQ4, ettTSEQMAX
};

enum {
    epctISOTROPIC, epctSEMIISOTROPIC, epctANISOTROPIC,
    epctSURFACETENSION, epctNR
};

enum {
    erscNO, erscALL, erscCOM, erscNR
};

enum {
    ecutsGROUP, ecutsVERLET, ecutsNR
};

/* Coulomb / VdW interaction modifiers.
 * grompp replaces eintmodPOTSHIFT_VERLET by eintmodPOTSHIFT or eintmodNONE.
 * Exactcutoff is only used by Reaction-field-zero, and is not user-selectable.
 */
enum eintmod {
    eintmodPOTSHIFT_VERLET, eintmodPOTSHIFT, eintmodNONE, eintmodPOTSWITCH, eintmodEXACTCUTOFF, eintmodNR
};

/*
 * eelNOTUSED1 used to be GB, but to enable generalized born with different
 * forms of electrostatics (RF, switch, etc.) in the future it is now selected
 * separately (through the implicit_solvent option).
 */
enum {
    eelCUT,     eelRF,     eelGRF,   eelPME,  eelEWALD,  eelP3M_AD,
    eelPOISSON, eelSWITCH, eelSHIFT, eelUSER, eelGB_NOTUSED, eelRF_NEC, eelENCADSHIFT,
    eelPMEUSER, eelPMESWITCH, eelPMEUSERSWITCH, eelRF_ZERO, eelNR
};

/* Ewald geometry */
enum {
    eewg3D, eewg3DC, eewgNR
};

#define EEL_RF(e) ((e) == eelRF || (e) == eelGRF || (e) == eelRF_NEC || (e) == eelRF_ZERO )

#define EEL_PME(e)  ((e) == eelPME || (e) == eelPMESWITCH || (e) == eelPMEUSER || (e) == eelPMEUSERSWITCH || (e) == eelP3M_AD)
#define EEL_FULL(e) (EEL_PME(e) || (e) == eelPOISSON || (e) == eelEWALD)

#define EEL_SWITCHED(e) ((e) == eelSWITCH || (e) == eelSHIFT || (e) == eelENCADSHIFT || (e) == eelPMESWITCH || (e) == eelPMEUSERSWITCH)

#define EEL_USER(e) ((e) == eelUSER || (e) == eelPMEUSER || (e) == (eelPMEUSERSWITCH))

#define EEL_IS_ZERO_AT_CUTOFF(e) (EEL_SWITCHED(e) || (e) == eelRF_ZERO)

#define EEL_MIGHT_BE_ZERO_AT_CUTOFF(e) (EEL_IS_ZERO_AT_CUTOFF(e) || (e) == eelUSER || (e) == eelPMEUSER)

enum {
    evdwCUT, evdwSWITCH, evdwSHIFT, evdwUSER, evdwENCADSHIFT, evdwNR
};

#define EVDW_SWITCHED(e) ((e) == evdwSWITCH || (e) == evdwSHIFT || (e) == evdwENCADSHIFT)

#define EVDW_IS_ZERO_AT_CUTOFF(e) EVDW_SWITCHED(e)

#define EVDW_MIGHT_BE_ZERO_AT_CUTOFF(e) (EVDW_IS_ZERO_AT_CUTOFF(e) || (e) == evdwUSER)

enum {
    ensGRID, ensSIMPLE, ensNR
};

/* eiVV is normal velocity verlet -- eiVVAK uses 1/2*(KE(t-dt/2)+KE(t+dt/2)) as the kinetic energy, and the half step kinetic
   energy for temperature control */

enum {
    eiMD, eiSteep, eiCG, eiBD, eiSD2, eiNM, eiLBFGS, eiTPI, eiTPIC, eiSD1, eiVV, eiVVAK, eiNR
};
#define EI_VV(e) ((e) == eiVV || (e) == eiVVAK)
#define EI_MD(e) ((e) == eiMD || EI_VV(e))
#define EI_SD(e) ((e) == eiSD1 || (e) == eiSD2)
#define EI_RANDOM(e) (EI_SD(e) || (e) == eiBD)
/*above integrators may not conserve momenta*/
#define EI_DYNAMICS(e) (EI_MD(e) || EI_SD(e) || (e) == eiBD)
#define EI_ENERGY_MINIMIZATION(e) ((e) == eiSteep || (e) == eiCG || (e) == eiLBFGS)
#define EI_TPI(e) ((e) == eiTPI || (e) == eiTPIC)

#define EI_STATE_VELOCITY(e) (EI_MD(e) || EI_SD(e))

enum {
    econtLINCS, econtSHAKE, econtNR
};

enum {
    edrNone, edrSimple, edrEnsemble, edrNR
};

enum {
    edrwConservative, edrwEqual, edrwNR
};

/* Combination rule things */
enum {
    eCOMB_NONE, eCOMB_GEOMETRIC, eCOMB_ARITHMETIC, eCOMB_GEOM_SIG_EPS, eCOMB_NR
};

/* NBF selection */
enum {
    eNBF_NONE, eNBF_LJ, eNBF_BHAM, eNBF_NR
};

/* simulated tempering methods */
enum {
    esimtempGEOMETRIC, esimtempEXPONENTIAL, esimtempLINEAR, esimtempNR
};
/* FEP selection */
enum {
    efepNO, efepYES, efepSTATIC, efepSLOWGROWTH, efepEXPANDED, efepNR
};
/* if efepNO, there are no evaluations at other states.
   if efepYES, treated equivalently to efepSTATIC.
   if efepSTATIC, then lambdas do not change during the simulation.
   if efepSLOWGROWTH, then the states change monotonically throughout the simulation.
   if efepEXPANDED, then expanded ensemble simulations are occuring.
 */

/* FEP coupling types */
enum {
    efptFEP, efptMASS, efptCOUL, efptVDW, efptBONDED, efptRESTRAINT, efptTEMPERATURE, efptNR
};

/* How the lambda weights are calculated:
   elamstatsMETROPOLIS = using the metropolis criteria
   elamstatsBARKER = using the Barker critera for transition weights - also called unoptimized Bennett
   elamstatsMINVAR = using Barker + minimum variance for weights
   elamstatsWL = Wang-Landu (using visitation counts)
   elamstatsWWL = Weighted Wang-Landau (using optimized gibbs weighted visitation counts)
 */
enum {
    elamstatsNO, elamstatsMETROPOLIS, elamstatsBARKER, elamstatsMINVAR, elamstatsWL, elamstatsWWL, elamstatsNR
};

#define ELAMSTATS_EXPANDED(e) ((e) > elamstatsNO)

#define EWL(e) ((e) == elamstatsWL || (e) == elamstatsWWL)

/* How moves in lambda are calculated:
   elmovemcMETROPOLIS - using the Metropolis criteria, and 50% up and down
   elmovemcBARKER - using the Barker criteria, and 50% up and down
   elmovemcGIBBS - computing the transition using the marginalized probabilities of the lambdas
   elmovemcMETGIBBS - computing the transition using the metropolized version of Gibbs (Monte Carlo Strategies in Scientific computing, Liu, p. 134)
 */
enum {
    elmcmoveNO, elmcmoveMETROPOLIS, elmcmoveBARKER, elmcmoveGIBBS, elmcmoveMETGIBBS, elmcmoveNR
};

/* how we decide whether weights have reached equilibrium
   elmceqNO - never stop, weights keep going
   elmceqYES - fix the weights from the beginning; no movement
   elmceqWLDELTA - stop when the WL-delta falls below a certain level
   elmceqNUMATLAM - stop when we have a certain number of samples at every step
   elmceqSTEPS - stop when we've run a certain total number of steps
   elmceqSAMPLES - stop when we've run a certain total number of samples
   elmceqRATIO - stop when the ratio of samples (lowest to highest) is sufficiently large
 */
enum {
    elmceqNO, elmceqYES, elmceqWLDELTA, elmceqNUMATLAM, elmceqSTEPS, elmceqSAMPLES, elmceqRATIO, elmceqNR
};

/* separate_dhdl_file selection */
enum
{
    /* NOTE: YES is the first one. Do NOT interpret this one as a gmx_bool */
    esepdhdlfileYES, esepdhdlfileNO, esepdhdlfileNR
};

/* dhdl_derivatives selection */
enum
{
    /* NOTE: YES is the first one. Do NOT interpret this one as a gmx_bool */
    edhdlderivativesYES, edhdlderivativesNO, edhdlderivativesNR
};

/* Solvent model */
enum {
    esolNO, esolSPC, esolTIP4P, esolNR
};

/* Dispersion correction */
enum {
    edispcNO, edispcEnerPres, edispcEner, edispcAllEnerPres, edispcAllEner, edispcNR
};

/* Shell types, for completion stuff */
enum {
    eshellCSH, eshellBASH, eshellZSH, eshellNR
};

/* Center of mass motion selection */
enum {
    ecmLINEAR, ecmANGULAR, ecmNO, ecmNR
};

/* New version of simulated annealing */
enum {
    eannNO, eannSINGLE, eannPERIODIC, eannNR
};

/* Implicit solvent algorithms */
enum {
    eisNO, eisGBSA, eisNR
};

/* Algorithms for calculating GB radii */
enum {
    egbSTILL, egbHCT, egbOBC, egbNR
};

enum {
    esaAPPROX, esaNO, esaSTILL, esaNR
};

/* Wall types */
enum {
    ewt93, ewt104, ewtTABLE, ewt126, ewtNR
};

/* Pull stuff */
enum {
    epullNO, epullUMBRELLA, epullCONSTRAINT, epullCONST_F, epullNR
};

enum {
    epullgDIST, epullgDIR, epullgCYL, epullgPOS, epullgDIRPBC, epullgNR
};

#define PULL_CYL(pull) ((pull)->eGeom == epullgCYL)

/* Enforced rotation groups */
enum {
    erotgISO, erotgISOPF,
    erotgPM, erotgPMPF,
    erotgRM, erotgRMPF,
    erotgRM2, erotgRM2PF,
    erotgFLEX, erotgFLEXT,
    erotgFLEX2, erotgFLEX2T,
    erotgNR
};

enum {
    erotgFitRMSD, erotgFitNORM, erotgFitPOT, erotgFitNR
};

/* QMMM */
enum {
    eQMmethodAM1, eQMmethodPM3, eQMmethodRHF,
    eQMmethodUHF, eQMmethodDFT, eQMmethodB3LYP, eQMmethodMP2, eQMmethodCASSCF, eQMmethodB3LYPLAN,
    eQMmethodDIRECT, eQMmethodNR
};

enum {
    eQMbasisSTO3G, eQMbasisSTO3G2, eQMbasis321G,
    eQMbasis321Gp, eQMbasis321dGp, eQMbasis621G,
    eQMbasis631G, eQMbasis631Gp, eQMbasis631dGp,
    eQMbasis6311G, eQMbasisNR
};

enum {
    eQMMMschemenormal, eQMMMschemeoniom, eQMMMschemeNR
};

enum {
    eMultentOptName, eMultentOptNo, eMultentOptLast, eMultentOptNR
};

/* flat-bottom posres geometries */
enum {
    efbposresZERO, efbposresSPHERE, efbposresCYLINDER, efbposresX, efbposresY, efbposresZ,
    efbposresNR
};

enum {
    eAdressOff, eAdressConst, eAdressXSplit, eAdressSphere, eAdressNR
};

enum {
    eAdressICOff, eAdressICThermoForce, eAdressICNR
};

enum {
    eAdressSITEcom, eAdressSITEcog, eAdressSITEatom, eAdressSITEatomatom, eAdressSITENR
};


/* The interactions contained in a (possibly merged) table
 * for computing electrostatic, VDW repulsion and/or VDW dispersion
 * contributions.
 */
enum gmx_table_interaction
{
    GMX_TABLE_INTERACTION_ELEC,
    GMX_TABLE_INTERACTION_VDWREP_VDWDISP,
    GMX_TABLE_INTERACTION_VDWEXPREP_VDWDISP,
    GMX_TABLE_INTERACTION_VDWDISP,
    GMX_TABLE_INTERACTION_ELEC_VDWREP_VDWDISP,
    GMX_TABLE_INTERACTION_ELEC_VDWEXPREP_VDWDISP,
    GMX_TABLE_INTERACTION_ELEC_VDWDISP,
    GMX_TABLE_INTERACTION_NR
};

/* Different formats for table data. Cubic spline tables are typically stored
 * with the four Y,F,G,H intermediate values (check tables.c for format), which
 * makes it easy to load with a single 4-way SIMD instruction too.
 * Linear tables only need one value per table point, or two if both V and F
 * are calculated. However, with SIMD instructions this makes the loads unaligned,
 * and in that case we store the data as F, D=F(i+1)-F(i), V, and then a blank value,
 * which again makes it possible to load as a single instruction.
 */
enum gmx_table_format
{
    GMX_TABLE_FORMAT_CUBICSPLINE_YFGH,
    GMX_TABLE_FORMAT_LINEAR_VF,
    GMX_TABLE_FORMAT_LINEAR_V,
    GMX_TABLE_FORMAT_LINEAR_F,
    GMX_TABLE_FORMAT_LINEAR_FDV0,
    GMX_TABLE_FORMAT_NR
};

/* Neighborlist geometry type.
 * Kernels will compute interactions between two particles,
 * 3-center water, 4-center water or coarse-grained beads.
 */
enum gmx_nblist_kernel_geometry
{
    GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE,
    GMX_NBLIST_GEOMETRY_WATER3_PARTICLE,
    GMX_NBLIST_GEOMETRY_WATER3_WATER3,
    GMX_NBLIST_GEOMETRY_WATER4_PARTICLE,
    GMX_NBLIST_GEOMETRY_WATER4_WATER4,
    GMX_NBLIST_GEOMETRY_CG_CG,
    GMX_NBLIST_GEOMETRY_NR
};

/* Types of electrostatics calculations available inside nonbonded kernels.
 * Note that these do NOT necessarily correspond to the user selections in the MDP file;
 * many interactions for instance map to tabulated kernels.
 */
enum gmx_nbkernel_elec
{
    GMX_NBKERNEL_ELEC_NONE,
    GMX_NBKERNEL_ELEC_COULOMB,
    GMX_NBKERNEL_ELEC_REACTIONFIELD,
    GMX_NBKERNEL_ELEC_CUBICSPLINETABLE,
    GMX_NBKERNEL_ELEC_GENERALIZEDBORN,
    GMX_NBKERNEL_ELEC_EWALD,
    GMX_NBKERNEL_ELEC_NR
};

/* Types of vdw calculations available inside nonbonded kernels.
 * Note that these do NOT necessarily correspond to the user selections in the MDP file;
 * many interactions for instance map to tabulated kernels.
 */
enum gmx_nbkernel_vdw
{
    GMX_NBKERNEL_VDW_NONE,
    GMX_NBKERNEL_VDW_LENNARDJONES,
    GMX_NBKERNEL_VDW_BUCKINGHAM,
    GMX_NBKERNEL_VDW_CUBICSPLINETABLE,
    GMX_NBKERNEL_VDW_NR
};
/* Types of interactions inside the neighborlist
 */
enum gmx_nblist_interaction_type
{
    GMX_NBLIST_INTERACTION_STANDARD,
    GMX_NBLIST_INTERACTION_FREE_ENERGY,
    GMX_NBLIST_INTERACTION_ADRESS,
    GMX_NBLIST_INTERACTION_NR
};

!!fcdata.h

typedef real rvec5[5];

/* Distance restraining stuff */
typedef struct {
    int      dr_weighting; /* Weighting of pairs in one restraint              */
    gmx_bool dr_bMixed;    /* Use sqrt of the instantaneous times              *
                            * the time averaged violation                      */
    real     dr_fc;        /* Force constant for disres,                       *
                            * which is multiplied by a (possibly)              *
                            * different factor for each restraint              */
    real  dr_tau;          /* Time constant for disres		          */
    real  ETerm;           /* multiplication factor for time averaging         */
    real  ETerm1;          /* 1 - ETerm1                                       */
    real  exp_min_t_tau;   /* Factor for slowly switching on the force         */
    int   nres;            /* The number of distance restraints                */
    int   npair;           /* The number of distance restraint pairs           */
    real  sumviol;         /* The sum of violations                            */
    real *rt;              /* The calculated instantaneous distance (npr)      */
    real *rm3tav;          /* The calculated time averaged distance (npr)      */
    real *Rtl_6;           /* The calculated instantaneous r^-6 (nr)           */
    real *Rt_6;            /* The calculated inst. ens. averaged r^-6 (nr)     */
    real *Rtav_6;          /* The calculated time and ens. averaged r^-6 (nr)  */
    int   nsystems;        /* The number of systems for ensemble averaging     */
} t_disresdata;


/* Orientation restraining stuff */
typedef struct {
    real      fc;            /* Force constant for the restraints                  */
    real      edt;           /* Multiplication factor for time averaging           */
    real      edt_1;         /* 1 - edt                                            */
    real      exp_min_t_tau; /* Factor for slowly switching on the force         */
    int       nr;            /* The number of orientation restraints               */
    int       nex;           /* The number of experiments                          */
    int       nref;          /* The number of atoms for the fit                    */
    real     *mref;          /* The masses of the reference atoms                  */
    rvec     *xref;          /* The reference coordinates for the fit (nref)       */
    rvec     *xtmp;          /* Temporary array for fitting (nref)                 */
    matrix    R;             /* Rotation matrix to rotate to the reference coor.   */
    tensor   *S;             /* Array of order tensors for each experiment (nexp)  */
    rvec5    *Dinsl;         /* The order matrix D for all restraints (nr x 5)     */
    rvec5    *Dins;          /* The ensemble averaged D (nr x 5)                   */
    rvec5    *Dtav;          /* The time and ensemble averaged D (nr x 5)          */
    real     *oinsl;         /* The calculated instantaneous orientations          */
    real     *oins;          /* The calculated emsemble averaged orientations      */
    real     *otav;          /* The calculated time and ensemble averaged orient.  */
    real      rmsdev;        /* The weighted (using kfac) RMS deviation            */
    rvec5    *tmp;           /* An array of temporary 5-vectors (nex);             */
    real   ***TMP;           /* An array of temporary 5x5 matrices (nex);          */
    real     *eig;           /* Eigenvalues/vectors, for output only (nex x 12)    */

    /* variables for diagonalization with diagonalize_orires_tensors()*/
    double **M;
    double  *eig_diag;
    double **v;
} t_oriresdata;

/*
 * Data struct used in the force calculation routines
 * for storing the tables for bonded interactions and
 * for storing information which is needed in following steps
 * (for instance for time averaging in distance retraints)
 * or for storing output, since force routines only return the potential.
 */
typedef struct {
    bondedtable_t *bondtab;
    bondedtable_t *angletab;
    bondedtable_t *dihtab;

    t_disresdata   disres;
    t_oriresdata   orires;
} t_fcdata;

!!filenm.h

enum {
    efMDP, efGCT,
    efTRX, efTRO, efTRN, efTRR, efTRJ, efXTC, efG87,
    efEDR,
    efSTX, efSTO, efGRO, efG96, efPDB, efBRK, efENT, efESP, efPQR, efXYZ,
    efCPT,
    efLOG, efXVG, efOUT,
    efNDX,
    efTOP, efITP,
    efTPX, efTPS, efTPR, efTPA, efTPB,
    efTEX, efRTP, efATP, efHDB,
    efDAT, efDLG,
    efMAP, efEPS, efMAT, efM2P,
    efMTX,
    efEDI,
    efHAT,
    efCUB,
    efXPM,
    efRND,
    efNR
};

typedef struct {
    int           ftp;    /* File type (see enum above)		*/
    const char   *opt;    /* Command line option			*/
    const char   *fn;     /* File name (as set in source code)	*/
    unsigned long flag;   /* Flag for all kinds of info (see defs)*/
    int           nfiles; /* number of files			*/
    char        **fns;    /* File names				*/
} t_filenm;

#define ffSET   1<<0
#define ffREAD  1<<1
#define ffWRITE 1<<2
#define ffOPT   1<<3
#define ffLIB   1<<4
#define ffMULT  1<<5
#define ffRW    (ffREAD | ffWRITE)
#define ffOPTRD (ffREAD | ffOPT)
#define ffOPTWR (ffWRITE| ffOPT)
#define ffOPTRW (ffRW   | ffOPT)
#define ffLIBRD (ffREAD | ffLIB)
#define ffLIBOPTRD (ffOPTRD | ffLIB)
#define ffRDMULT   (ffREAD  | ffMULT)
#define ffOPTRDMULT   (ffRDMULT | ffOPT)
#define ffWRMULT   (ffWRITE  | ffMULT)
#define ffOPTWRMULT   (ffWRMULT | ffOPT)

!!force_flags.h

/* The state has changed */
#define GMX_FORCE_STATECHANGED (1<<0)
/* The box might have changed */
#define GMX_FORCE_DYNAMICBOX   (1<<1)
/* Do neighbor searching */
#define GMX_FORCE_NS           (1<<2)
/* Update long-range neighborlists */
#define GMX_FORCE_LRNS         (1<<3)
/* Calculate bonded energies/forces */
#define GMX_FORCE_BONDED       (1<<4)
/* Store long-range forces in a separate array */
#define GMX_FORCE_SEPLRF       (1<<5)
/* Calculate non-bonded energies/forces */
#define GMX_FORCE_NONBONDED    (1<<6)
/* Calculate forces (not only energies) */
#define GMX_FORCE_FORCES       (1<<7)
/* Calculate the virial */
#define GMX_FORCE_VIRIAL       (1<<8)
/* Calculate energies */
#define GMX_FORCE_ENERGY       (1<<9)
/* Calculate dHdl */
#define GMX_FORCE_DHDL         (1<<10)
/* Calculate long-range energies/forces */
#define GMX_FORCE_DO_LR        (1<<11)

/* Normally one want all energy terms and forces */
#define GMX_FORCE_ALLFORCES    (GMX_FORCE_BONDED | GMX_FORCE_NONBONDED | GMX_FORCE_FORCES)

!!forcerec.h

#include "ns.h"
#include "genborn.h"
#include "qmmmrec.h"
#include "idef.h"
#include "nb_verlet.h"
#include "interaction_const.h"
#include "hw_info.h"

typedef struct gmx_pme *gmx_pme_t;



/* Structure describing the data in a single table */
typedef struct
{
    enum gmx_table_interaction  interaction; /* Types of interactions stored in this table */
    enum gmx_table_format       format;      /* Interpolation type and data format */

    real                        r;           /* range of the table */
    int                         n;           /* n+1 is the number of table points */
    real                        scale;       /* distance (nm) between two table points */
    real                        scale_exp;   /* distance for exponential part of VdW table, not always used */
    real *                      data;        /* the actual table data */

    /* Some information about the table layout. This can also be derived from the interpolation
     * type and the table interactions, but it is convenient to have here for sanity checks, and it makes it
     * much easier to access the tables in the nonbonded kernels when we can set the data from variables.
     * It is always true that stride = formatsize*ninteractions
     */
    int                         formatsize;    /* Number of fp variables for each table point (1 for F, 2 for VF, 4 for YFGH, etc.) */
    int                         ninteractions; /* Number of interactions in table, 1 for coul-only, 3 for coul+rep+disp. */
    int                         stride;        /* Distance to next table point (number of fp variables per table point in total) */
} t_forcetable;

typedef struct
{
    t_forcetable   table_elec;
    t_forcetable   table_vdw;
    t_forcetable   table_elec_vdw;

    /* The actual neighbor lists, short and long range, see enum above
     * for definition of neighborlist indices.
     */
    t_nblist nlist_sr[eNL_NR];
    t_nblist nlist_lr[eNL_NR];
} t_nblists;

/* macros for the cginfo data in forcerec */
/* The maximum cg size in cginfo is 63
 * because we only have space for 6 bits in cginfo,
 * this cg size entry is actually only read with domain decomposition.
 * But there is a smaller limit due to the t_excl data structure
 * which is defined in nblist.h.
 */
#define SET_CGINFO_GID(cgi, gid)      (cgi) = (((cgi)  &  ~65535)  |  (gid)   )
#define GET_CGINFO_GID(cgi)        ( (cgi)            &   65535)
#define SET_CGINFO_EXCL_INTRA(cgi)   (cgi) =  ((cgi)  |  (1<<16))
#define GET_CGINFO_EXCL_INTRA(cgi) ( (cgi)            &  (1<<16))
#define SET_CGINFO_EXCL_INTER(cgi)   (cgi) =  ((cgi)  |  (1<<17))
#define GET_CGINFO_EXCL_INTER(cgi) ( (cgi)            &  (1<<17))
#define SET_CGINFO_SOLOPT(cgi, opt)   (cgi) = (((cgi)  & ~(3<<18)) | ((opt)<<18))
#define GET_CGINFO_SOLOPT(cgi)     (((cgi)>>18)       &   3)
#define SET_CGINFO_CONSTR(cgi)       (cgi) =  ((cgi)  |  (1<<20))
#define GET_CGINFO_CONSTR(cgi)     ( (cgi)            &  (1<<20))
#define SET_CGINFO_SETTLE(cgi)       (cgi) =  ((cgi)  |  (1<<21))
#define GET_CGINFO_SETTLE(cgi)     ( (cgi)            &  (1<<21))
/* This bit is only used with bBondComm in the domain decomposition */
#define SET_CGINFO_BOND_INTER(cgi)   (cgi) =  ((cgi)  |  (1<<22))
#define GET_CGINFO_BOND_INTER(cgi) ( (cgi)            &  (1<<22))
#define SET_CGINFO_HAS_VDW(cgi)      (cgi) =  ((cgi)  |  (1<<23))
#define GET_CGINFO_HAS_VDW(cgi)    ( (cgi)            &  (1<<23))
#define SET_CGINFO_HAS_Q(cgi)        (cgi) =  ((cgi)  |  (1<<24))
#define GET_CGINFO_HAS_Q(cgi)      ( (cgi)            &  (1<<24))
#define SET_CGINFO_NATOMS(cgi, opt)   (cgi) = (((cgi)  & ~(63<<25)) | ((opt)<<25))
#define GET_CGINFO_NATOMS(cgi)     (((cgi)>>25)       &   63)


/* Value to be used in mdrun for an infinite cut-off.
 * Since we need to compare with the cut-off squared,
 * this value should be slighlty smaller than sqrt(GMX_FLOAT_MAX).
 */
#define GMX_CUTOFF_INF 1E+18

/* enums for the neighborlist type */
enum {
    enbvdwNONE, enbvdwLJ, enbvdwBHAM, enbvdwTAB, enbvdwNR
};
/* OOR is "one over r" -- standard coul */
enum {
    enbcoulNONE, enbcoulOOR, enbcoulRF, enbcoulTAB, enbcoulGB, enbcoulFEWALD, enbcoulNR
};

enum {
    egCOULSR, egLJSR, egBHAMSR, egCOULLR, egLJLR, egBHAMLR,
    egCOUL14, egLJ14, egGB, egNR
};

typedef struct {
    int   nener;      /* The number of energy group pairs     */
    real *ener[egNR]; /* Energy terms for each pair of groups */
} gmx_grppairener_t;

typedef struct {
    real              term[F_NRE];         /* The energies for all different interaction types */
    gmx_grppairener_t grpp;
    double            dvdl_lin[efptNR];    /* Contributions to dvdl with linear lam-dependence */
    double            dvdl_nonlin[efptNR]; /* Idem, but non-linear dependence                  */
    int               n_lambda;
    int               fep_state;           /*current fep state -- just for printing */
    double           *enerpart_lambda;     /* Partial energy for lambda and flambda[] */
    real              foreign_term[F_NRE]; /* alternate array for storing foreign lambda energies */
    gmx_grppairener_t foreign_grpp;        /* alternate array for storing foreign lambda energies */
} gmx_enerdata_t;
/* The idea is that dvdl terms with linear lambda dependence will be added
 * automatically to enerpart_lambda. Terms with non-linear lambda dependence
 * should explicitly determine the energies at foreign lambda points
 * when n_lambda > 0.
 */

typedef struct {
    int  cg_start;
    int  cg_end;
    int  cg_mod;
    int *cginfo;
} cginfo_mb_t;


/* ewald table type */
typedef struct ewald_tab *ewald_tab_t;

typedef struct {
    rvec             *f;
    int               f_nalloc;
    unsigned          red_mask; /* Mask for marking which parts of f are filled */
    rvec             *fshift;
    real              ener[F_NRE];
    gmx_grppairener_t grpp;
    real              Vcorr;
    real              dvdl[efptNR];
    tensor            vir;
} f_thread_t;

typedef struct {
    interaction_const_t *ic;

    /* Domain Decomposition */
    gmx_bool bDomDec;

    /* PBC stuff */
    int                  ePBC;
    gmx_bool             bMolPBC;
    int                  rc_scaling;
    rvec                 posres_com;
    rvec                 posres_comB;

    const gmx_hw_info_t *hwinfo;
    gmx_bool             use_cpu_acceleration;

    /* Interaction for calculated in kernels. In many cases this is similar to
     * the electrostatics settings in the inputrecord, but the difference is that
     * these variables always specify the actual interaction in the kernel - if
     * we are tabulating reaction-field the inputrec will say reaction-field, but
     * the kernel interaction will say cubic-spline-table. To be safe we also
     * have a kernel-specific setting for the modifiers - if the interaction is
     * tabulated we already included the inputrec modification there, so the kernel
     * modification setting will say 'none' in that case.
     */
    int nbkernel_elec_interaction;
    int nbkernel_vdw_interaction;
    int nbkernel_elec_modifier;
    int nbkernel_vdw_modifier;

    /* Use special N*N kernels? */
    gmx_bool bAllvsAll;
    /* Private work data */
    void    *AllvsAll_work;
    void    *AllvsAll_workgb;

    /* Cut-Off stuff.
     * Infinite cut-off's will be GMX_CUTOFF_INF (unlike in t_inputrec: 0).
     */
    real rlist, rlistlong;

    /* Dielectric constant resp. multiplication factor for charges */
    real zsquare, temp;
    real epsilon_r, epsilon_rf, epsfac;

    /* Constants for reaction fields */
    real kappa, k_rf, c_rf;

    /* Charge sum and dipole for topology A/B ([0]/[1]) for Ewald corrections */
    double qsum[2];
    double q2sum[2];
    rvec   mu_tot[2];

    /* Dispersion correction stuff */
    int  eDispCorr;

    /* The shift of the shift or user potentials */
    real enershiftsix;
    real enershifttwelve;
    /* Integrated differces for energy and virial with cut-off functions */
    real enerdiffsix;
    real enerdifftwelve;
    real virdiffsix;
    real virdifftwelve;
    /* Constant for long range dispersion correction (average dispersion)
     * for topology A/B ([0]/[1]) */
    real avcsix[2];
    /* Constant for long range repulsion term. Relative difference of about
     * 0.1 percent with 0.8 nm cutoffs. But hey, it's cheap anyway...
     */
    real avctwelve[2];

    /* Fudge factors */
    real fudgeQQ;

    /* Table stuff */
    gmx_bool     bcoultab;
    gmx_bool     bvdwtab;
    /* The normal tables are in the nblists struct(s) below */
    t_forcetable tab14; /* for 1-4 interactions only */

    /* PPPM & Shifting stuff */
    int   coulomb_modifier;
    real  rcoulomb_switch, rcoulomb;
    real *phi;

    /* VdW stuff */
    int    vdw_modifier;
    double reppow;
    real   rvdw_switch, rvdw;
    real   bham_b_max;

    /* Free energy */
    int      efep;
    real     sc_alphavdw;
    real     sc_alphacoul;
    int      sc_power;
    real     sc_r_power;
    real     sc_sigma6_def;
    real     sc_sigma6_min;
    gmx_bool bSepDVDL;

    /* NS Stuff */
    int  eeltype;
    int  vdwtype;
    int  cg0, hcg;
    /* solvent_opt contains the enum for the most common solvent
     * in the system, which will be optimized.
     * It can be set to esolNO to disable all water optimization */
    int          solvent_opt;
    int          nWatMol;
    gmx_bool     bGrid;
    gmx_bool     bExcl_IntraCGAll_InterCGNone;
    cginfo_mb_t *cginfo_mb;
    int         *cginfo;
    rvec        *cg_cm;
    int          cg_nalloc;
    rvec        *shift_vec;

    /* The neighborlists including tables */
    int                 nnblists;
    int                *gid2nblists;
    t_nblists          *nblists;

    int                 cutoff_scheme; /* group- or Verlet-style cutoff */
    gmx_bool            bNonbonded;    /* true if nonbonded calculations are *not* turned off */
    nonbonded_verlet_t *nbv;

    /* The wall tables (if used) */
    int            nwall;
    t_forcetable **wall_tab;

    /* The number of charge groups participating in do_force_lowlevel */
    int ncg_force;
    /* The number of atoms participating in do_force_lowlevel */
    int natoms_force;
    /* The number of atoms participating in force and constraints */
    int natoms_force_constr;
    /* The allocation size of vectors of size natoms_force */
    int nalloc_force;

    /* Twin Range stuff, f_twin has size natoms_force */
    gmx_bool bTwinRange;
    int      nlr;
    rvec    *f_twin;

    /* Forces that should not enter into the virial summation:
     * PPPM/PME/Ewald/posres
     */
    gmx_bool bF_NoVirSum;
    int      f_novirsum_n;
    int      f_novirsum_nalloc;
    rvec    *f_novirsum_alloc;
    /* Pointer that points to f_novirsum_alloc when pressure is calcaluted,
     * points to the normal force vectors wen pressure is not requested.
     */
    rvec *f_novirsum;

    /* Long-range forces and virial for PPPM/PME/Ewald */
    gmx_pme_t pmedata;
    tensor    vir_el_recip;

    /* PME/Ewald stuff */
    gmx_bool    bEwald;
    real        ewaldcoeff;
    ewald_tab_t ewald_table;

    /* Virial Stuff */
    rvec *fshift;
    rvec  vir_diag_posres;
    dvec  vir_wall_z;

    /* Non bonded Parameter lists */
    int      ntype; /* Number of atom types */
    gmx_bool bBHAM;
    real    *nbfp;

    /* Energy group pair flags */
    int *egp_flags;

    /* xmdrun flexible constraints */
    real fc_stepsize;

    /* Generalized born implicit solvent */
    gmx_bool       bGB;
    /* Generalized born stuff */
    real           gb_epsilon_solvent;
    /* Table data for GB */
    t_forcetable   gbtab;
    /* VdW radius for each atomtype (dim is thus ntype) */
    real          *atype_radius;
    /* Effective radius (derived from effective volume) for each type */
    real          *atype_vol;
    /* Implicit solvent - surface tension for each atomtype */
    real          *atype_surftens;
    /* Implicit solvent - radius for GB calculation */
    real          *atype_gb_radius;
    /* Implicit solvent - overlap for HCT model */
    real          *atype_S_hct;
    /* Generalized born interaction data */
    gmx_genborn_t *born;

    /* Table scale for GB */
    real gbtabscale;
    /* Table range for GB */
    real gbtabr;
    /* GB neighborlists (the sr list will contain for each atom all other atoms
     * (for use in the SA calculation) and the lr list will contain
     * for each atom all atoms 1-4 or greater (for use in the GB calculation)
     */
    t_nblist gblist_sr;
    t_nblist gblist_lr;
    t_nblist gblist;

    /* Inverse square root of the Born radii for implicit solvent */
    real *invsqrta;
    /* Derivatives of the potential with respect to the Born radii */
    real *dvda;
    /* Derivatives of the Born radii with respect to coordinates */
    real *dadx;
    real *dadx_rawptr;
    int   nalloc_dadx; /* Allocated size of dadx */

    /* If > 0 signals Test Particle Insertion,
     * the value is the number of atoms of the molecule to insert
     * Only the energy difference due to the addition of the last molecule
     * should be calculated.
     */
    gmx_bool n_tpi;

    /* Neighbor searching stuff */
    gmx_ns_t ns;

    /* QMMM stuff */
    gmx_bool         bQMMM;
    t_QMMMrec       *qr;

    /* QM-MM neighborlists */
    t_nblist QMMMlist;

    /* Limit for printing large forces, negative is don't print */
    real print_force;

    /* coarse load balancing time measurement */
    double t_fnbf;
    double t_wait;
    int    timesteps;

    /* parameter needed for AdResS simulation */
    int             adress_type;
    gmx_bool        badress_tf_full_box;
    real            adress_const_wf;
    real            adress_ex_width;
    real            adress_hy_width;
    int             adress_icor;
    int             adress_site;
    rvec            adress_refs;
    int             n_adress_tf_grps;
    int           * adress_tf_table_index;
    int            *adress_group_explicit;
    t_forcetable *  atf_tabs;
    real            adress_ex_forcecap;
    gmx_bool        adress_do_hybridpairs;

    /* User determined parameters, copied from the inputrec */
    int  userint1;
    int  userint2;
    int  userint3;
    int  userint4;
    real userreal1;
    real userreal2;
    real userreal3;
    real userreal4;

    /* Thread local force and energy data */
    /* FIXME move to bonded_thread_data_t */
    int         nthreads;
    int         red_ashift;
    int         red_nblock;
    f_thread_t *f_t;

    /* Exclusion load distribution over the threads */
    int  *excl_load;
} t_forcerec;

/* Important: Starting with Gromacs-4.6, the values of c6 and c12 in the nbfp array have
 * been scaled by 6.0 or 12.0 to save flops in the kernels. We have corrected this everywhere
 * in the code, but beware if you are using these macros externally.
 */
#define C6(nbfp, ntp, ai, aj)     (nbfp)[2*((ntp)*(ai)+(aj))]
#define C12(nbfp, ntp, ai, aj)    (nbfp)[2*((ntp)*(ai)+(aj))+1]
#define BHAMC(nbfp, ntp, ai, aj)  (nbfp)[3*((ntp)*(ai)+(aj))]
#define BHAMA(nbfp, ntp, ai, aj)  (nbfp)[3*((ntp)*(ai)+(aj))+1]
#define BHAMB(nbfp, ntp, ai, aj)  (nbfp)[3*((ntp)*(ai)+(aj))+2]

!!genborn.h

#include"simple.h"
typedef struct
{
    int  nbonds;
    int  bond[10];
    real length[10];
} genborn_bonds_t;

typedef struct gbtmpnbls *gbtmpnbls_t;

/* Struct to hold all the information for GB */
typedef struct
{
    int nr;                   /* number of atoms, length of arrays below */
    int n12;                  /* number of 1-2 (bond) interactions       */
    int n13;                  /* number of 1-3 (angle) terms             */
    int n14;                  /* number of 1-4 (torsion) terms           */
    int nalloc;               /* Allocation of local arrays (with DD)    */


    /* Arrays below that end with _globalindex are used for setting up initial values of
     * all gb parameters and values. They all have length natoms, which for DD is the
     * global atom number.
     * Values are then taken from these arrays to local copies, that have names without
     * _globalindex, in the routine make_local_gb(), which is called once for single
     * node runs, and for DD at every call to dd_partition_system
     */

    real       *gpol;              /* Atomic polarisation energies */
    real       *gpol_globalindex;  /*  */
    real       *gpol_still_work;   /* Work array for Still model */
    real       *gpol_hct_work;     /* Work array for HCT/OBC models */
    real       *bRad;              /* Atomic Born radii */
    real       *vsolv;             /* Atomic solvation volumes */
    real       *vsolv_globalindex; /*  */
    real       *gb_radius;         /* Radius info, copied from atomtypes */
    real       *gb_radius_globalindex;

    int        *use;                /* Array that till if this atom does GB */
    int        *use_globalindex;    /* Global array for parallelization */

    real        es;                 /* Solvation energy and derivatives */
    real       *asurf;              /* Atomic surface area */
    rvec       *dasurf;             /* Surface area derivatives */
    real        as;                 /* Total surface area */

    real       *drobc;              /* Parameters for OBC chain rule calculation */
    real       *param;              /* Precomputed factor rai*atype->S_hct for HCT/OBC */
    real       *param_globalindex;  /*  */

    real       *log_table;          /* Table for logarithm lookup */

    real        obc_alpha;          /* OBC parameters */
    real        obc_beta;           /* OBC parameters */
    real        obc_gamma;          /* OBC parameters */
    real        gb_doffset;         /* Dielectric offset for Still/HCT/OBC */
    real        gb_epsilon_solvent; /*   */
    real        epsilon_r;          /* Used for inner dielectric */

    real        sa_surface_tension; /* Surface tension for non-polar solvation */

    real       *work;               /* Used for parallel summation and in the chain rule, length natoms         */
    real       *buf;                /* Used for parallel summation and in the chain rule, length natoms         */
    int        *count;              /* Used for setting up the special gb nblist, length natoms                 */
    gbtmpnbls_t nblist_work;        /* Used for setting up the special gb nblist, dim natoms*nblist_work_nalloc */
    int         nblist_work_nalloc; /* Length of second dimension of nblist_work                                */
} gmx_genborn_t;

!!globsig.h

enum {
    eglsNABNSB, eglsCHKPT, eglsSTOPCOND, eglsRESETCOUNTERS, eglsNR
};

typedef struct {
    int nstms;       /* The frequency for intersimulation communication */
    int sig[eglsNR]; /* The signal set by one process in do_md */
    int set[eglsNR]; /* The communicated signal, equal for all processes */
} globsig_t;

!!graph.h
#include "idef.h"

typedef enum {
    egcolWhite, egcolGrey, egcolBlack, egcolNR
} egCol;

typedef struct {
    int          at0;       /* The first atom the graph was constructed for */
    int          at1;       /* The last atom the graph was constructed for */
    int          nnodes;    /* The number of nodes, nnodes=at_end-at_start	*/
    int          nbound;    /* The number of nodes with edges		*/
    int          at_start;  /* The first connected atom in this graph	*/
    int          at_end;    /* The last+1 connected atom in this graph	*/
    int         *nedge;     /* For each node the number of edges		*/
    atom_id    **edge;      /* For each node, the actual edges (bidirect.)	*/
    gmx_bool     bScrewPBC; /* Screw boundary conditions                    */
    ivec        *ishift;    /* Shift for each particle                  */
    int          negc;
    egCol       *egc;       /* color of each node */
} t_graph;


#define SHIFT_IVEC(g, i) ((g)->ishift[i])

!!group.h
#include "simple.h"

typedef struct {
    real    Th;             /* Temperature at half step        */
    real    T;              /* Temperature at full step        */
    tensor  ekinh;          /* Kinetic energy at half step     */
    tensor  ekinh_old;      /* Kinetic energy at old half step */
    tensor  ekinf;          /* Kinetic energy at full step     */
    real    lambda;         /* Berendsen coupling lambda       */
    double  ekinscalef_nhc; /* Scaling factor for NHC- full step */
    double  ekinscaleh_nhc; /* Scaling factor for NHC- half step */
    double  vscale_nhc;     /* Scaling factor for NHC- velocity */
} t_grp_tcstat;

typedef struct {
    int     nat;    /* Number of atoms in this group		*/
    rvec    u;      /* Mean velocities of home particles        */
    rvec    uold;   /* Previous mean velocities of home particles   */
    double  mA;     /* Mass for topology A		                */
    double  mB;     /* Mass for topology B		                */
} t_grp_acc;

typedef struct {
    real    cos_accel;  /* The acceleration for the cosine profile      */
    real    mvcos;      /* The cos momenta of home particles            */
    real    vcos;       /* The velocity of the cosine profile           */
} t_cos_acc;

typedef struct {
    gmx_bool         bNEMD;
    int              ngtc;            /* The number of T-coupling groups      */
    t_grp_tcstat    *tcstat;          /* T-coupling data            */
    tensor         **ekin_work_alloc; /* Allocated locations for *_work members */
    tensor         **ekin_work;       /* Work arrays for tcstat per thread    */
    real           **dekindl_work;    /* Work location for dekindl per thread */
    int              ngacc;           /* The number of acceleration groups    */
    t_grp_acc       *grpstat;         /* Acceleration data			*/
    tensor           ekin;            /* overall kinetic energy               */
    tensor           ekinh;           /* overall 1/2 step kinetic energy      */
    real             dekindl;         /* dEkin/dlambda at half step           */
    real             dekindl_old;     /* dEkin/dlambda at old half step       */
    t_cos_acc        cosacc;          /* Cosine acceleration data             */
} gmx_ekindata_t;

#define GID(igid, jgid, gnr) ((igid < jgid) ? (igid*gnr+jgid) : (jgid*gnr+igid))

!!hw_info.h
#include "simple.h"
#include "nbnxn_cuda_types_ext.h"
#include "../gmx_cpuid.h"

typedef enum
{
    egpuCompatible = 0,  egpuNonexistent,  egpuIncompatible, egpuInsane
} e_gpu_detect_res_t;

/* Textual names of the GPU detection/check results (see e_gpu_detect_res_t). */
static const char * const gpu_detect_res_str[] =
{
    "compatible", "inexistent", "incompatible", "insane"
};

/* GPU device information -- for now with only CUDA devices.
 * The gmx_hardware_detect module initializes it. */
typedef struct
{
    gmx_bool             bUserSet;      /* true if the GPUs in cuda_dev_use are manually provided by the user */

    int                  ncuda_dev_use; /* number of devices selected to be used */
    int                 *cuda_dev_use;  /* index of the devices selected to be used */
    int                  ncuda_dev;     /* total number of devices detected */
    cuda_dev_info_ptr_t  cuda_dev;      /* devices detected in the system (per node) */
} gmx_gpu_info_t;

/* Hardware information structure with CPU and GPU information.
 * It is initialized by gmx_detect_hardware().
 * NOTE: this structure may only contain structures that are globally valid
 *       (i.e. must be able to be shared among all threads) */
typedef struct
{
    gmx_bool        bCanUseGPU;          /* True if compatible GPUs are detected during hardware detection */
    gmx_gpu_info_t  gpu_info;            /* Information about GPUs detected in the system */

    gmx_cpuid_t     cpuid_info;          /* CPUID information about CPU detected;
                                            NOTE: this will only detect the CPU thread 0 of the
                                            current process runs on. */
    int             nthreads_hw_avail;   /* Number of hardware threads available; this number
                                            is based on the number of CPUs reported as available
                                            by the OS at the time of detection. */
    gmx_bool        bConsistencyChecked; /* whether
                                            gmx_check_hw_runconf_consistency()
                                            has been run with this hw_info */
} gmx_hw_info_t;

!!idef.h

#include "simple.h"

/* check kernel/toppush.c when you change these numbers */
#define MAXATOMLIST 6
#define MAXFORCEPARAM   12
#define NR_RBDIHS   6
#define NR_FOURDIHS     4

typedef atom_id t_iatom;

/* this MUST correspond to the
   t_interaction_function[F_NRE] in gmxlib/ifunc.c */
enum {
    F_BONDS,
    F_G96BONDS,
    F_MORSE,
    F_CUBICBONDS,
    F_CONNBONDS,
    F_HARMONIC,
    F_FENEBONDS,
    F_TABBONDS,
    F_TABBONDSNC,
    F_RESTRBONDS,
    F_ANGLES,
    F_G96ANGLES,
    F_LINEAR_ANGLES,
    F_CROSS_BOND_BONDS,
    F_CROSS_BOND_ANGLES,
    F_UREY_BRADLEY,
    F_QUARTIC_ANGLES,
    F_TABANGLES,
    F_PDIHS,
    F_RBDIHS,
    F_FOURDIHS,
    F_IDIHS,
    F_PIDIHS,
    F_TABDIHS,
    F_CMAP,
    F_GB12,
    F_GB13,
    F_GB14,
    F_GBPOL,
    F_NPSOLVATION,
    F_LJ14,
    F_COUL14,
    F_LJC14_Q,
    F_LJC_PAIRS_NB,
    F_LJ,
    F_BHAM,
    F_LJ_LR,
    F_BHAM_LR,
    F_DISPCORR,
    F_COUL_SR,
    F_COUL_LR,
    F_RF_EXCL,
    F_COUL_RECIP,
    F_DPD,
    F_POLARIZATION,
    F_WATER_POL,
    F_THOLE_POL,
    F_ANHARM_POL,
    F_POSRES,
    F_FBPOSRES,
    F_DISRES,
    F_DISRESVIOL,
    F_ORIRES,
    F_ORIRESDEV,
    F_ANGRES,
    F_ANGRESZ,
    F_DIHRES,
    F_DIHRESVIOL,
    F_CONSTR,
    F_CONSTRNC,
    F_SETTLE,
    F_VSITE2,
    F_VSITE3,
    F_VSITE3FD,
    F_VSITE3FAD,
    F_VSITE3OUT,
    F_VSITE4FD,
    F_VSITE4FDN,
    F_VSITEN,
    F_COM_PULL,
    F_EQM,
    F_EPOT,
    F_EKIN,
    F_ETOT,
    F_ECONSERVED,
    F_TEMP,
    F_VTEMP_NOLONGERUSED,
    F_PDISPCORR,
    F_PRES,
    F_DVDL_CONSTR,
    F_DVDL,
    F_DKDL,
    F_DVDL_COUL,
    F_DVDL_VDW,
    F_DVDL_BONDED,
    F_DVDL_RESTRAINT,
    F_DVDL_TEMPERATURE, /* not calculated for now, but should just be the energy (NVT) or enthalpy (NPT), or 0 (NVE) */
    F_NRE               /* This number is for the total number of energies	*/
};

#define IS_RESTRAINT_TYPE(ifunc) (((ifunc == F_POSRES) || (ifunc == F_DISRES) || (ifunc == F_RESTRBONDS) || (ifunc == F_DISRESVIOL) || (ifunc == F_ORIRES) || (ifunc == F_ORIRESDEV) || (ifunc == F_ANGRES) || (ifunc == F_ANGRESZ) || (ifunc == F_DIHRES)))

/* A macro for checking if ftype is an explicit pair-listed LJ or COULOMB
 * interaction type:
 * bonded LJ (usually 1-4), or special listed non-bonded for FEP.
 */
#define IS_LISTED_LJ_C(ftype) ((ftype) >= F_LJ14 && (ftype) <= F_LJC_PAIRS_NB)

typedef union
{
    /* Some parameters have A and B values for free energy calculations.
     * The B values are not used for regular simulations of course.
     * Free Energy for nonbondeds can be computed by changing the atom type.
     * The harmonic type is used for all harmonic potentials:
     * bonds, angles and improper dihedrals
     */
    struct {
        real a, b, c;
    } bham;
    struct {
        real rA, krA, rB, krB;
    } harmonic;
    struct {
        real klinA, aA, klinB, aB;
    } linangle;
    struct {
        real lowA, up1A, up2A, kA, lowB, up1B, up2B, kB;
    } restraint;
    /* No free energy supported for cubic bonds, FENE, WPOL or cross terms */
    struct {
        real b0, kb, kcub;
    } cubic;
    struct {
        real bm, kb;
    } fene;
    struct {
        real r1e, r2e, krr;
    } cross_bb;
    struct {
        real r1e, r2e, r3e, krt;
    } cross_ba;
    struct {
        real thetaA, kthetaA, r13A, kUBA, thetaB, kthetaB, r13B, kUBB;
    } u_b;
    struct {
        real theta, c[5];
    } qangle;
    struct {
        real alpha;
    } polarize;
    struct {
        real alpha, drcut, khyp;
    } anharm_polarize;
    struct {
        real al_x, al_y, al_z, rOH, rHH, rOD;
    } wpol;
    struct {
        real a, alpha1, alpha2, rfac;
    } thole;
    struct {
        real c6, c12;
    } lj;
    struct {
        real c6A, c12A, c6B, c12B;
    } lj14;
    struct {
        real fqq, qi, qj, c6, c12;
    } ljc14;
    struct {
        real qi, qj, c6, c12;
    } ljcnb;
    /* Proper dihedrals can not have different multiplicity when
     * doing free energy calculations, because the potential would not
     * be periodic anymore.
     */
    struct {
        real phiA, cpA; int mult; real phiB, cpB;
    } pdihs;
    struct {
        real dA, dB;
    } constr;
    /* Settle can not be used for Free energy calculations of water bond geometry.
     * Use shake (or lincs) instead if you have to change the water bonds.
     */
    struct {
        real doh, dhh;
    } settle;
    struct {
        real b0A, cbA, betaA, b0B, cbB, betaB;
    } morse;
    struct {
        real pos0A[DIM], fcA[DIM], pos0B[DIM], fcB[DIM];
    } posres;
    struct {
        real pos0[DIM], r, k; int geom;
    } fbposres;
    struct {
        real rbcA[NR_RBDIHS], rbcB[NR_RBDIHS];
    } rbdihs;
    struct {
        real a, b, c, d, e, f;
    } vsite;
    struct {
        int  n; real a;
    } vsiten;
    /* NOTE: npair is only set after reading the tpx file */
    struct {
        real low, up1, up2, kfac; int type, label, npair;
    } disres;
    struct {
        real phiA, dphiA, kfacA, phiB, dphiB, kfacB;
    } dihres;
    struct {
        int  ex, power, label; real c, obs, kfac;
    } orires;
    struct {
        int  table; real kA; real kB;
    } tab;
    struct {
        real sar, st, pi, gbr, bmlt;
    } gb;
    struct {
        int cmapA, cmapB;
    } cmap;
    struct {
        real buf[MAXFORCEPARAM];
    } generic;                                               /* Conversion */
} t_iparams;

typedef int t_functype;

/*
 * The nonperturbed/perturbed interactions are now separated (sorted) in the
 * ilist, such that the first 0..(nr_nonperturbed-1) ones are exactly that, and
 * the remaining ones from nr_nonperturbed..(nr-1) are perturbed bonded
 * interactions.
 */
typedef struct
{
    int      nr;
    int      nr_nonperturbed;
    t_iatom *iatoms;
    int      nalloc;
} t_ilist;

/*
 * The struct t_ilist defines a list of atoms with their interactions.
 * General field description:
 *   int nr
 *	the size (nr elements) of the interactions array (iatoms[]).
 *   t_iatom *iatoms
 *  specifies which atoms are involved in an interaction of a certain
 *       type. The layout of this array is as follows:
 *
 *	  +-----+---+---+---+-----+---+---+-----+---+---+---+-----+---+---+...
 *	  |type1|at1|at2|at3|type2|at1|at2|type1|at1|at2|at3|type3|at1|at2|
 *	  +-----+---+---+---+-----+---+---+-----+---+---+---+-----+---+---+...
 *
 *  So for interaction type type1 3 atoms are needed, and for type2 and
 *      type3 only 2. The type identifier is used to select the function to
 *	calculate the interaction and its actual parameters. This type
 *	identifier is an index in a params[] and functype[] array.
 */

typedef struct
{
    real *cmap; /* Has length 4*grid_spacing*grid_spacing, */
    /* there are 4 entries for each cmap type (V,dVdx,dVdy,d2dVdxdy) */
} cmapdata_t;

typedef struct
{
    int         ngrid;        /* Number of allocated cmap (cmapdata_t ) grids */
    int         grid_spacing; /* Grid spacing */
    cmapdata_t *cmapdata;     /* Pointer to grid with actual, pre-interpolated data */
} gmx_cmap_t;


typedef struct
{
    int         ntypes;
    int         atnr;
    t_functype *functype;
    t_iparams  *iparams;
    double      reppow;    /* The repulsion power for VdW: C12*r^-reppow   */
    real        fudgeQQ;   /* The scaling factor for Coulomb 1-4: f*q1*q2  */
    gmx_cmap_t  cmap_grid; /* The dihedral correction maps                 */
} gmx_ffparams_t;

enum {
    ilsortUNKNOWN, ilsortNO_FE, ilsortFE_UNSORTED, ilsortFE_SORTED
};

typedef struct
{
    int         ntypes;
    int         atnr;
    t_functype *functype;
    t_iparams  *iparams;
    real        fudgeQQ;
    gmx_cmap_t  cmap_grid;
    t_iparams  *iparams_posres, *iparams_fbposres;
    int         iparams_posres_nalloc, iparams_fbposres_nalloc;

    t_ilist     il[F_NRE];
    int         ilsort;
} t_idef;

/*
 * The struct t_idef defines all the interactions for the complete
 * simulation. The structure is setup in such a way that the multinode
 * version of the program  can use it as easy as the single node version.
 * General field description:
 *   int ntypes
 *	defines the number of elements in functype[] and param[].
 *   int nodeid
 *      the node id (if parallel machines)
 *   int atnr
 *      the number of atomtypes
 *   t_functype *functype
 *	array of length ntypes, defines for every force type what type of
 *      function to use. Every "bond" with the same function but different
 *	force parameters is a different force type. The type identifier in the
 *	forceatoms[] array is an index in this array.
 *   t_iparams *iparams
 *	array of length ntypes, defines the parameters for every interaction
 *      type. The type identifier in the actual interaction list
 *      (ilist[ftype].iatoms[]) is an index in this array.
 *   gmx_cmap_t cmap_grid
 *      the grid for the dihedral pair correction maps.
 *   t_iparams *iparams_posres, *iparams_fbposres
 *	defines the parameters for position restraints only.
 *      Position restraints are the only interactions that have different
 *      parameters (reference positions) for different molecules
 *      of the same type. ilist[F_POSRES].iatoms[] is an index in this array.
 *   t_ilist il[F_NRE]
 *      The list of interactions for each type. Note that some,
 *      such as LJ and COUL will have 0 entries.
 */

typedef struct {
    int   n;      /* n+1 is the number of points */
    real  scale;  /* distance between two points */
    real *data;   /* the actual table data, per point there are 4 numbers */
} bondedtable_t;

!!ifunc.h

#include "idef.h"
#include "mdatom.h"
#include "fcdata.h"
#include "graph.h"
#include "pbc.h"

typedef real t_ifunc (int nbonds, const t_iatom iatoms[],
                      const t_iparams iparams[],
                      const rvec x[], rvec f[], rvec fshift[],
                      const t_pbc *pbc, const t_graph *g,
                      real lambda, real *dvdlambda,
                      const t_mdatoms *md, t_fcdata *fcd,
                      int *ddgatindex);

/*
 * The function type t_ifunc() calculates one interaction, using iatoms[]
 * and iparams. Within the function the number of atoms to be used is
 * known. Within the function only the atomid part of the iatoms[] array
 * is supplied, not the type field (see also t_ilist). The function
 * returns the potential energy. If pbc==NULL the coordinates in x are
 * assumed to be such that no calculation of PBC is necessary,
 * If pbc!=NULL a full PBC calculation is performed.
 * If g!=NULL it is used for determining the shift forces.
 * With domain decomposition ddgatindex can be used for getting global
 * atom numbers for warnings and error messages.
 * ddgatindex is NULL when domain decomposition is not used.
 */

#define IF_NULL       0
#define IF_BOND       1
#define IF_VSITE      1<<1
#define IF_CONSTRAINT 1<<2
#define IF_CHEMBOND   1<<3
#define IF_BTYPE      1<<4
#define IF_ATYPE      1<<5
#define IF_TABULATED  1<<6
#define IF_LIMZERO    1<<7
/* These flags tell to some of the routines what can be done with this
 * item in the list.
 * With IF_BOND a bonded interaction will be calculated.
 * With IF_BTYPE grompp can convert the bond to a Morse potential.
 * With IF_BTYPE or IF_ATYPE the bond/angle can be converted to
 * a constraint or used for vsite parameter determination by grompp.
 * IF_LIMZERO indicates that for a bonded interaction the potential
 * does goes to zero for large distances, thus if such an interaction
 * it not assigned to any node by the domain decompostion, the simulation
 * still continue, if mdrun has been told so.
 */
typedef struct
{
    const char *name;         /* the name of this function			*/
    const char *longname;     /* The name for printing etc.                   */
    int         nratoms;      /* nr of atoms needed for this function		*/
    int         nrfpA, nrfpB; /* number of parameters for this function.      */
                              /* this corresponds to the number of params in  */
                              /* iparams struct! (see idef.h)                 */
    /* A and B are for normal and free energy components respectively.    */
    unsigned long   flags;    /* Flags (see above)                            */
    int             nrnb_ind; /* index for nrnb (-1 if unknown)               */
    t_ifunc        *ifunc;    /* the function it self				*/
} t_interaction_function;

#define NRFPA(ftype) (interaction_function[(ftype)].nrfpA)
#define NRFPB(ftype) (interaction_function[(ftype)].nrfpB)
#define NRFP(ftype)  (NRFPA(ftype)+NRFPB(ftype))
#define NRAL(ftype) (interaction_function[(ftype)].nratoms)

#define IS_CHEMBOND(ftype) (interaction_function[(ftype)].nratoms == 2 && (interaction_function[(ftype)].flags & IF_CHEMBOND))
/* IS_CHEMBOND tells if function type ftype represents a chemical bond */

/* IS_ANGLE tells if a function type ftype represents an angle
 * Per Larsson, 2007-11-06
 */
#define IS_ANGLE(ftype) (interaction_function[(ftype)].nratoms == 3 && (interaction_function[(ftype)].flags & IF_ATYPE))
#define IS_VSITE(ftype) (interaction_function[(ftype)].flags & IF_VSITE)

#define IS_TABULATED(ftype) (interaction_function[(ftype)].flags & IF_TABULATED)

extern const t_interaction_function interaction_function[F_NRE];
/* initialised interaction functions descriptor		

!!inputrec.h

#include "simple.h"
#include "../sysstuff.h"

typedef struct {
    int   n;    /* Number of terms				*/
    real *a;    /* Coeffients (V / nm )                     */
    real *phi;  /* Phase angles					*/
} t_cosines;

typedef struct {
    real E0;            /* Field strength (V/nm)                        */
    real omega;         /* Frequency (1/ps)                             */
    real t0;            /* Centre of the Gaussian pulse (ps)            */
    real sigma;         /* Width of the Gaussian pulse (FWHM) (ps)      */
} t_efield;

#define EGP_EXCL  (1<<0)
#define EGP_TABLE (1<<1)

typedef struct {
    int       ngtc;           /* # T-Coupl groups                        */
    int       nhchainlength;  /* # of nose-hoover chains per group       */
    int       ngacc;          /* # Accelerate groups                     */
    int       ngfrz;          /* # Freeze groups                         */
    int       ngener;         /* # Ener groups			    */
    real     *nrdf;           /* Nr of degrees of freedom in a group	    */
    real     *ref_t;          /* Coupling temperature	per group   */
    int      *annealing;      /* No/simple/periodic SA for each group    */
    int      *anneal_npoints; /* Number of annealing time points per grp */
    real    **anneal_time;    /* For ea. group: Time points              */
    real    **anneal_temp;    /* For ea. grp: Temperature at these times */
                              /* Final temp after all intervals is ref_t */
    real     *tau_t;          /* Tau coupling time              */
    rvec     *acc;            /* Acceleration per group		    */
    ivec     *nFreeze;        /* Freeze the group in each direction ?    */
    int      *egp_flags;      /* Exclusions/tables of energy group pairs */

    /* QMMM stuff */
    int          ngQM;         /* nr of QM groups                              */
    int         *QMmethod;     /* Level of theory in the QM calculation        */
    int         *QMbasis;      /* Basisset in the QM calculation               */
    int         *QMcharge;     /* Total charge in the QM region                */
    int         *QMmult;       /* Spin multiplicicty in the QM region          */
    gmx_bool    *bSH;          /* surface hopping (diabatic hop only)          */
    int         *CASorbitals;  /* number of orbiatls in the active space       */
    int         *CASelectrons; /* number of electrons in the active space      */
    real        *SAon;         /* at which gap (A.U.) the SA is switched on    */
    real        *SAoff;
    int         *SAsteps;      /* in how many steps SA goes from 1-1 to 0.5-0.5*/
    gmx_bool    *bOPT;
    gmx_bool    *bTS;
} t_grpopts;

enum {
    epgrppbcNONE, epgrppbcREFAT, epgrppbcCOS
};

typedef struct {
    int         nat;        /* Number of atoms in the pull group */
    atom_id    *ind;        /* The global atoms numbers */
    int         nat_loc;    /* Number of local pull atoms */
    int         nalloc_loc; /* Allocation size for ind_loc and weight_loc */
    atom_id    *ind_loc;    /* Local pull indices */
    int         nweight;    /* The number of weights (0 or nat) */
    real       *weight;     /* Weights (use all 1 when weight==NULL) */
    real       *weight_loc; /* Weights for the local indices */
    int         epgrppbc;   /* The type of pbc for this pull group, see enum above */
    atom_id     pbcatom;    /* The reference atom for pbc (global number) */
    rvec        vec;        /* The pull vector, direction or position */
    rvec        init;       /* Initial reference displacement */
    real        rate;       /* Rate of motion (nm/ps) */
    real        k;          /* force constant */
    real        kB;         /* force constant for state B */
    real        wscale;     /* scaling factor for the weights: sum w m/sum w w m */
    real        invtm;      /* inverse total mass of the group: 1/wscale sum w m */
    dvec        x;          /* center of mass before update */
    dvec        xp;         /* center of mass after update before constraining */
    dvec        dr;         /* The distance from the reference group */
    double      f_scal;     /* Scalar force for directional pulling */
    dvec        f;          /* force due to the pulling/constraining */
} t_pullgrp;

typedef struct {
    int   eSimTempScale; /* simulated temperature scaling; linear or exponential */
    real  simtemp_low;   /* the low temperature for simulated tempering  */
    real  simtemp_high;  /* the high temperature for simulated tempering */
    real *temperatures;  /* the range of temperatures used for simulated tempering */
} t_simtemp;

typedef struct {
    int    nstdhdl;                 /* The frequency for calculating dhdl           */
    double init_lambda;             /* fractional value of lambda (usually will use
                                       init_fep_state, this will only be for slow growth,
                                       and for legacy free energy code. Only has a
                                       valid value if positive)   */
    int      init_fep_state;        /* the initial number of the state                   */
    double   delta_lambda;          /* change of lambda per time step (fraction of (0.1) */
    gmx_bool bPrintEnergy;          /* Whether to print the energy in the dhdl           */
    int      n_lambda;              /* The number of foreign lambda points               */
    double **all_lambda;            /* The array of all lambda values                    */
    int      lambda_neighbors;      /* The number of neighboring lambda states to
                                       calculate the energy for in up and down directions
                                       (-1 for all) */
    int      lambda_start_n;        /* The first lambda to calculate energies for */
    int      lambda_stop_n;         /* The last lambda +1 to calculate energies for */
    real     sc_alpha;              /* free energy soft-core parameter                   */
    int      sc_power;              /* lambda power for soft-core interactions           */
    real     sc_r_power;            /* r power for soft-core interactions                */
    real     sc_sigma;              /* free energy soft-core sigma when c6 or c12=0      */
    real     sc_sigma_min;          /* free energy soft-core sigma for ?????             */
    gmx_bool bScCoul;               /* use softcore for the coulomb portion as well (default FALSE) */
    gmx_bool separate_dvdl[efptNR]; /* whether to print the dvdl term associated with
                                       this term; if it is not specified as separate,
                                       it is lumped with the FEP term */
    int    separate_dhdl_file;      /* whether to write a separate dhdl.xvg file
                                       note: NOT a gmx_bool, but an enum */
    int    dhdl_derivatives;        /* whether to calculate+write dhdl derivatives
                                       note: NOT a gmx_bool, but an enum */
    int    dh_hist_size;            /* The maximum table size for the dH histogram */
    double dh_hist_spacing;         /* The spacing for the dH histogram */
} t_lambda;

typedef struct {
    int      nstexpanded;         /* The frequency of expanded ensemble state changes */
    int      elamstats;           /* which type of move updating do we use for lambda monte carlo (or no for none) */
    int      elmcmove;            /* what move set will be we using for state space moves */
    int      elmceq;              /* the method we use to decide of we have equilibrated the weights */
    int      equil_n_at_lam;      /* the minumum number of samples at each lambda for deciding whether we have reached a minimum */
    real     equil_wl_delta;      /* WL delta at which we stop equilibrating weights */
    real     equil_ratio;         /* use the ratio of weights (ratio of minimum to maximum) to decide when to stop equilibrating */
    int      equil_steps;         /* after equil_steps steps we stop equilibrating the weights */
    int      equil_samples;       /* after equil_samples total samples (steps/nstfep), we stop equilibrating the weights */
    int      lmc_seed;            /* random number seed for lambda mc switches */
    gmx_bool minvar;              /* whether to use minumum variance weighting */
    int      minvarmin;           /* the number of samples needed before kicking into minvar routine */
    real     minvar_const;        /* the offset for the variance in MinVar */
    int      c_range;             /* range of cvalues used for BAR */
    gmx_bool bSymmetrizedTMatrix; /* whether to print symmetrized matrices */
    int      nstTij;              /* How frequently to print the transition matrices */
    int      lmc_repeats;         /* number of repetitions in the MC lambda jumps */  /*MRS -- VERIFY THIS */
    int      lmc_forced_nstart;   /* minimum number of samples for each state before free sampling */ /* MRS -- VERIFY THIS! */
    int      gibbsdeltalam;       /* distance in lambda space for the gibbs interval */
    real     wl_scale;            /* scaling factor for wang-landau */
    real     wl_ratio;            /* ratio between largest and smallest number for freezing the weights */
    real     init_wl_delta;       /* starting delta for wang-landau */
    gmx_bool bWLoneovert;         /* use one over t convergece for wang-landau when the delta get sufficiently small */
    gmx_bool bInit_weights;       /* did we initialize the weights? */
    real     mc_temp;             /* To override the main temperature, or define it if it's not defined */
    real    *init_lambda_weights; /* user-specified initial weights to start with  */
} t_expanded;

typedef struct {
    int            ngrp;       /* number of groups */
    int            eGeom;      /* pull geometry */
    ivec           dim;        /* used to select components for constraint */
    real           cyl_r1;     /* radius of cylinder for dynamic COM */
    real           cyl_r0;     /* radius of cylinder including switch length */
    real           constr_tol; /* absolute tolerance for constraints in (nm) */
    int            nstxout;    /* Output frequency for pull x */
    int            nstfout;    /* Output frequency for pull f */
    int            ePBC;       /* the boundary conditions */
    int            npbcdim;    /* do pbc in dims 0 <= dim < npbcdim */
    gmx_bool       bRefAt;     /* do we need reference atoms for a group COM ? */
    int            cosdim;     /* dimension for cosine weighting, -1 if none */
    gmx_bool       bVirial;    /* do we need to add the pull virial? */
    t_pullgrp     *grp;        /* groups to pull/restrain/etc/ */
    t_pullgrp     *dyna;       /* dynamic groups for use with local constraints */
    rvec          *rbuf;       /* COM calculation buffer */
    dvec          *dbuf;       /* COM calculation buffer */
    double        *dbuf_cyl;   /* cylinder ref. groups COM calculation buffer */

    FILE          *out_x;      /* output file for pull data */
    FILE          *out_f;      /* output file for pull data */
} t_pull;


/* Abstract types for enforced rotation only defined in pull_rotation.c       */
typedef struct gmx_enfrot *gmx_enfrot_t;
typedef struct gmx_enfrotgrp *gmx_enfrotgrp_t;

typedef struct {
    int         eType;             /* Rotation type for this group                  */
    int         bMassW;            /* Use mass-weighed positions?                   */
    int         nat;               /* Number of atoms in the group                  */
    atom_id    *ind;               /* The global atoms numbers                      */
    rvec       *x_ref;             /* The reference positions                       */
    rvec        vec;               /* The normalized rotation vector                */
    real        rate;              /* Rate of rotation (degree/ps)                  */
    real        k;                 /* Force constant (kJ/(mol nm^2)                 */
    rvec        pivot;             /* Pivot point of rotation axis (nm)             */
    int         eFittype;          /* Type of fit to determine actual group angle   */
    int         PotAngle_nstep;    /* Number of angles around the reference angle
                                      for which the rotation potential is also
                                      evaluated (for fit type 'potential' only)     */
    real            PotAngle_step; /* Distance between two angles in degrees (for
                                      fit type 'potential' only)                    */
    real            slab_dist;     /* Slab distance (nm)                            */
    real            min_gaussian;  /* Minimum value the gaussian must have so that
                                      the force is actually evaluated               */
    real            eps;           /* Additive constant for radial motion2 and
                                      flexible2 potentials (nm^2)                   */
    gmx_enfrotgrp_t enfrotgrp;     /* Stores non-inputrec rotation data per group   */
} t_rotgrp;

typedef struct {
    int          ngrp;       /* Number of rotation groups                     */
    int          nstrout;    /* Output frequency for main rotation outfile    */
    int          nstsout;    /* Output frequency for per-slab data            */
    t_rotgrp    *grp;        /* Groups to rotate                              */
    gmx_enfrot_t enfrot;     /* Stores non-inputrec enforced rotation data    */
} t_rot;


typedef struct {
    int      type;           /* type of AdResS simulation                    */
    gmx_bool bnew_wf;        /* enable new AdResS weighting function         */
    gmx_bool bchempot_dx;    /*true:interaction table format input is F=-dmu/dx   false: dmu_dwp  */
    gmx_bool btf_full_box;   /* true: appy therm force everywhere in the box according to table false: only in hybrid region */
    real     const_wf;       /* value of weighting function for eAdressConst */
    real     ex_width;       /* center of the explicit zone                  */
    real     hy_width;       /* width of the hybrid zone                     */
    int      icor;           /* type of interface correction                 */
    int      site;           /* AdResS CG site location                      */
    rvec     refs;           /* Coordinates for AdResS reference             */
    real     ex_forcecap;    /* in the hybrid zone, cap forces large then this to adress_ex_forcecap */
    gmx_bool do_hybridpairs; /* If true pair interaction forces are also scaled in an adress way*/

    int    * tf_table_index; /* contains mapping of energy group index -> i-th adress tf table*/
    int      n_tf_grps;
    int     *group_explicit;
    int      n_energy_grps;
} t_adress;

typedef struct {
    int             eI;                   /* Integration method                 */
    gmx_large_int_t nsteps;               /* number of steps to be taken			*/
    int             simulation_part;      /* Used in checkpointing to separate chunks */
    gmx_large_int_t init_step;            /* start at a stepcount >0 (used w. tpbconv)    */
    int             nstcalcenergy;        /* frequency of energy calc. and T/P coupl. upd.	*/
    int             cutoff_scheme;        /* group or verlet cutoffs     */
    int             ns_type;              /* which ns method should we use?               */
    int             nstlist;              /* number of steps before pairlist is generated	*/
    int             ndelta;               /* number of cells per rlong			*/
    int             nstcomm;              /* number of steps after which center of mass	*/
    /* motion is removed				*/
    int             comm_mode;            /* Center of mass motion removal algorithm      */
    int             nstcheckpoint;        /* checkpointing frequency                      */
    int             nstlog;               /* number of steps after which print to logfile	*/
    int             nstxout;              /* number of steps after which X is output	*/
    int             nstvout;              /* id. for V					*/
    int             nstfout;              /* id. for F					*/
    int             nstenergy;            /* number of steps after which energies printed */
    int             nstxtcout;            /* id. for compressed trj (.xtc)		*/
    double          init_t;               /* initial time (ps)              */
    double          delta_t;              /* time step (ps)				*/
    real            xtcprec;              /* precision of xtc file                        */
    real            fourier_spacing;      /* requested fourier_spacing, when nk? not set  */
    int             nkx, nky, nkz;        /* number of k vectors in each spatial dimension*/
                                          /* for fourier methods for long range electrost.*/
    int             pme_order;            /* interpolation order for PME                  */
    real            ewald_rtol;           /* Real space tolerance for Ewald, determines   */
                                          /* the real/reciprocal space relative weight    */
    int             ewald_geometry;       /* normal/3d ewald, or pseudo-2d LR corrections */
    real            epsilon_surface;      /* Epsilon for PME dipole correction            */
    gmx_bool        bOptFFT;              /* optimize the fft plan at start               */
    int             ePBC;                 /* Type of periodic boundary conditions		*/
    int             bPeriodicMols;        /* Periodic molecules                           */
    gmx_bool        bContinuation;        /* Continuation run: starting state is correct	*/
    int             etc;                  /* temperature coupling               */
    int             nsttcouple;           /* interval in steps for temperature coupling   */
    gmx_bool        bPrintNHChains;       /* whether to print nose-hoover chains        */
    int             epc;                  /* pressure coupling                            */
    int             epct;                 /* pressure coupling type			*/
    int             nstpcouple;           /* interval in steps for pressure coupling      */
    real            tau_p;                /* pressure coupling time (ps)			*/
    tensor          ref_p;                /* reference pressure (kJ/(mol nm^3))		*/
    tensor          compress;             /* compressability ((mol nm^3)/kJ)        */
    int             refcoord_scaling;     /* How to scale absolute reference coordinates  */
    rvec            posres_com;           /* The COM of the posres atoms                  */
    rvec            posres_comB;          /* The B-state COM of the posres atoms          */
    int             andersen_seed;        /* Random seed for Andersen thermostat (obsolete) */
    real            verletbuf_drift;      /* Max. drift (kJ/mol/ps/atom) for list buffer  */
    real            rlist;                /* short range pairlist cut-off (nm)		*/
    real            rlistlong;            /* long range pairlist cut-off (nm)		*/
    int             nstcalclr;            /* Frequency of evaluating direct space long-range interactions */
    real            rtpi;                 /* Radius for test particle insertion           */
    int             coulombtype;          /* Type of electrostatics treatment             */
    int             coulomb_modifier;     /* Modify the Coulomb interaction              */
    real            rcoulomb_switch;      /* Coulomb switch range start (nm)		*/
    real            rcoulomb;             /* Coulomb cutoff (nm)		                */
    real            epsilon_r;            /* relative dielectric constant                 */
    real            epsilon_rf;           /* relative dielectric constant of the RF       */
    int             implicit_solvent;     /* No (=explicit water), or GBSA solvent models */
    int             gb_algorithm;         /* Algorithm to use for calculation Born radii  */
    int             nstgbradii;           /* Frequency of updating Generalized Born radii */
    real            rgbradii;             /* Cutoff for GB radii calculation              */
    real            gb_saltconc;          /* Salt concentration (M) for GBSA models       */
    real            gb_epsilon_solvent;   /* dielectric coeff. of implicit solvent     */
    real            gb_obc_alpha;         /* 1st scaling factor for Bashford-Case GB      */
    real            gb_obc_beta;          /* 2nd scaling factor for Bashford-Case GB      */
    real            gb_obc_gamma;         /* 3rd scaling factor for Bashford-Case GB      */
    real            gb_dielectric_offset; /* Dielectric offset for Still/HCT/OBC     */
    int             sa_algorithm;         /* Algorithm for SA part of GBSA                */
    real            sa_surface_tension;   /* Energy factor for SA part of GBSA */
    int             vdwtype;              /* Type of Van der Waals treatment              */
    int             vdw_modifier;         /* Modify the VdW interaction                   */
    real            rvdw_switch;          /* Van der Waals switch range start (nm)        */
    real            rvdw;                 /* Van der Waals cutoff (nm)	        */
    int             eDispCorr;            /* Perform Long range dispersion corrections    */
    real            tabext;               /* Extension of the table beyond the cut-off,   *
                                           * as well as the table length for 1-4 interac. */
    real            shake_tol;            /* tolerance for shake				*/
    int             efep;                 /* free energy calculations                     */
    t_lambda       *fepvals;              /* Data for the FEP state                       */
    gmx_bool        bSimTemp;             /* Whether to do simulated tempering            */
    t_simtemp      *simtempvals;          /* Variables for simulated tempering            */
    gmx_bool        bExpanded;            /* Whether expanded ensembles are used          */
    t_expanded     *expandedvals;         /* Expanded ensemble parameters              */
    int             eDisre;               /* Type of distance restraining                 */
    real            dr_fc;                /* force constant for ta_disre			*/
    int             eDisreWeighting;      /* type of weighting of pairs in one restraints	*/
    gmx_bool        bDisreMixed;          /* Use comb of time averaged and instan. viol's	*/
    int             nstdisreout;          /* frequency of writing pair distances to enx   */
    real            dr_tau;               /* time constant for memory function in disres    */
    real            orires_fc;            /* force constant for orientational restraints  */
    real            orires_tau;           /* time constant for memory function in orires    */
    int             nstorireout;          /* frequency of writing tr(SD) to enx           */
    real            dihre_fc;             /* force constant for dihedral restraints (obsolete)	*/
    real            em_stepsize;          /* The stepsize for updating			*/
    real            em_tol;               /* The tolerance				*/
    int             niter;                /* Number of iterations for convergence of      */
                                          /* steepest descent in relax_shells             */
    real            fc_stepsize;          /* Stepsize for directional minimization        */
                                          /* in relax_shells                              */
    int             nstcgsteep;           /* number of steps after which a steepest       */
                                          /* descents step is done while doing cg         */
    int             nbfgscorr;            /* Number of corrections to the hessian to keep */
    int             eConstrAlg;           /* Type of constraint algorithm                 */
    int             nProjOrder;           /* Order of the LINCS Projection Algorithm      */
    real            LincsWarnAngle;       /* If bond rotates more than %g degrees, warn   */
    int             nLincsIter;           /* Number of iterations in the final Lincs step */
    gmx_bool        bShakeSOR;            /* Use successive overrelaxation for shake      */
    real            bd_fric;              /* Friction coefficient for BD (amu/ps)         */
    int             ld_seed;              /* Random seed for SD and BD                    */
    int             nwall;                /* The number of walls                          */
    int             wall_type;            /* The type of walls                            */
    real            wall_r_linpot;        /* The potentail is linear for r<=wall_r_linpot */
    int             wall_atomtype[2];     /* The atom type for walls                      */
    real            wall_density[2];      /* Number density for walls                     */
    real            wall_ewald_zfac;      /* Scaling factor for the box for Ewald         */
    int             ePull;                /* Type of pulling: no, umbrella or constraint  */
    t_pull         *pull;                 /* The data for center of mass pulling          */
    gmx_bool        bRot;                 /* Calculate enforced rotation potential(s)?    */
    t_rot          *rot;                  /* The data for enforced rotation potentials    */
    real            cos_accel;            /* Acceleration for viscosity calculation       */
    tensor          deform;               /* Triclinic deformation velocities (nm/ps)     */
    int             userint1;             /* User determined parameters                   */
    int             userint2;
    int             userint3;
    int             userint4;
    real            userreal1;
    real            userreal2;
    real            userreal3;
    real            userreal4;
    t_grpopts       opts;          /* Group options				*/
    t_cosines       ex[DIM];       /* Electric field stuff	(spatial part)		*/
    t_cosines       et[DIM];       /* Electric field stuff	(time part)		*/
    gmx_bool        bQMMM;         /* QM/MM calculation                            */
    int             QMconstraints; /* constraints on QM bonds                      */
    int             QMMMscheme;    /* Scheme: ONIOM or normal                      */
    real            scalefactor;   /* factor for scaling the MM charges in QM calc.*/
                                   /* parameter needed for AdResS simulation       */
    gmx_bool        bAdress;       /* Is AdResS enabled ? */
    t_adress       *adress;        /* The data for adress simulations */
} t_inputrec;

#define DEFORM(ir) ((ir).deform[XX][XX] != 0 || (ir).deform[YY][YY] != 0 || (ir).deform[ZZ][ZZ] != 0 || (ir).deform[YY][XX] != 0 || (ir).deform[ZZ][XX] != 0 || (ir).deform[ZZ][YY] != 0)

#define DYNAMIC_BOX(ir) ((ir).epc != epcNO || (ir).eI == eiTPI || DEFORM(ir))

#define PRESERVE_SHAPE(ir) ((ir).epc != epcNO && (ir).deform[XX][XX] == 0 && ((ir).epct == epctISOTROPIC || (ir).epct == epctSEMIISOTROPIC))

#define NEED_MUTOT(ir) (((ir).coulombtype == eelEWALD || EEL_PME((ir).coulombtype)) && ((ir).ewald_geometry == eewg3DC || (ir).epsilon_surface != 0))

#define IR_TWINRANGE(ir) ((ir).rlist > 0 && ((ir).rlistlong == 0 || (ir).rlistlong > (ir).rlist))

#define IR_ELEC_FIELD(ir) ((ir).ex[XX].n > 0 || (ir).ex[YY].n > 0 || (ir).ex[ZZ].n > 0)

#define IR_EXCL_FORCES(ir) (EEL_FULL((ir).coulombtype) || (EEL_RF((ir).coulombtype) && (ir).coulombtype != eelRF_NEC) || (ir).implicit_solvent != eisNO)
/* use pointer definitions of ir here, since that's what's usually used in the code */
#define IR_NPT_TROTTER(ir) ((((ir)->eI == eiVV) || ((ir)->eI == eiVVAK)) && (((ir)->epc == epcMTTK) && ((ir)->etc == etcNOSEHOOVER)))

#define IR_NVT_TROTTER(ir) ((((ir)->eI == eiVV) || ((ir)->eI == eiVVAK)) && ((!((ir)->epc == epcMTTK)) && ((ir)->etc == etcNOSEHOOVER)))

#define IR_NPH_TROTTER(ir) ((((ir)->eI == eiVV) || ((ir)->eI == eiVVAK)) && (((ir)->epc == epcMTTK) && (!(((ir)->etc == etcNOSEHOOVER)))))

!! interaction_const.h 
typedef struct {
    /* VdW */
    real rvdw;
    real sh_invrc6; /* For shifting the LJ potential */

    /* type of electrostatics (defined in enums.h) */
    int  eeltype;

    /* Coulomb */
    real rcoulomb;

    /* Cut-off */
    real rlist;
    real rlistlong;

    /* PME/Ewald */
    real ewaldcoeff;
    real sh_ewald;   /* For shifting the Ewald potential */

    /* Dielectric constant resp. multiplication factor for charges */
    real epsilon_r;
    real epsfac;

    /* Constants for reaction-field or plain cut-off */
    real epsilon_rf;
    real k_rf;
    real c_rf;

    /* Force/energy interpolation tables, linear in force, quadratic in V */
    real  tabq_scale;
    int   tabq_size;
    /* Coulomb force table, size of array is tabq_size (when used) */
    real *tabq_coul_F;
    /* Coulomb energy table, size of array is tabq_size (when used) */
    real *tabq_coul_V;
    /* Coulomb force+energy table, size of array is tabq_size*4,
       entry quadruplets are: F[i], F[i+1]-F[i], V[i], 0,
       this is used with single precision x86 SIMD for aligned loads */
    real *tabq_coul_FDV0;
} interaction_const_t;

!! ishift.h

/* not really neccesary, right now: */
#define D_BOX_Z 1
#define D_BOX_Y 1
#define D_BOX_X 2
#define N_BOX_Z (2*D_BOX_Z+1)
#define N_BOX_Y (2*D_BOX_Y+1)
#define N_BOX_X (2*D_BOX_X+1)
#define N_IVEC  (N_BOX_Z*N_BOX_Y*N_BOX_X)
#define CENTRAL (N_IVEC/2)
#define SHIFTS  N_IVEC

#define XYZ2IS(x, y, z) (N_BOX_X*(N_BOX_Y*((z)+D_BOX_Z)+(y)+D_BOX_Y)+(x)+D_BOX_X)
#define IVEC2IS(iv)   (XYZ2IS((iv)[XX], (iv)[YY], (iv)[ZZ]))
#define IS2X(iv)      (((iv) % N_BOX_X) - D_BOX_X)
#define IS2Y(iv)      ((((iv) / N_BOX_X) % N_BOX_Y) - D_BOX_Y)
#define IS2Z(iv)      ((iv) / (N_BOX_X*N_BOX_Y) - D_BOX_Z)

!! iteratedconstraints.h

/* iterate constraints up to 50 times  */
#define MAXITERCONST       50

/* data type */
typedef struct
{
    real     f, fprev, x, xprev;
    int      iter_i;
    gmx_bool bIterationActive;
    real     allrelerr[MAXITERCONST+2];
    int      num_close; /* number of "close" violations, caused by limited precision. */
} gmx_iterate_t;

void gmx_iterate_init(gmx_iterate_t *iterate, gmx_bool bIterate);

gmx_bool done_iterating(const t_commrec *cr, FILE *fplog, int nsteps, gmx_iterate_t *iterate, gmx_bool bFirstIterate, real fom, real *newf);

!!matrix.h

#include "simple.h"

typedef struct {
    real r, g, b;
} t_rgb;

typedef struct {
    char c1; /* should all be non-zero (and printable and not '"') */
    char c2; /*
              * should all be zero (single char color names: smaller xpm's)
              * or should all be non-zero (double char color names: more colors)
              */
} t_xpmelmt;

typedef short t_matelmt;

typedef struct {
    t_xpmelmt   code; /* see comment for t_xpmelmt */
    const char *desc;
    t_rgb       rgb;
} t_mapping;

#define MAT_SPATIAL_X (1<<0)
#define MAT_SPATIAL_Y (1<<1)
/* Defines if x and y are spatial dimensions,
 * when not, there are n axis ticks at the middle of the elements,
 * when set, there are n+1 axis ticks at the edges of the elements.
 */

typedef struct {
    unsigned int flags; /* The possible flags are defined above */
    int          nx, ny;
    int          y0;
    char         title[256];
    char         legend[256];
    char         label_x[256];
    char         label_y[256];
    gmx_bool     bDiscrete;
    real        *axis_x;
    real        *axis_y;
    t_matelmt  **matrix;
    int          nmap;
    t_mapping   *map;
} t_matrix;
/* title      matrix title
 * legend     label for the continuous legend
 * label_x    label for the x-axis
 * label_y    label for the y-axis
 * nx, ny     size of the matrix
 * axis_x[]   the x-ticklabels
 * axis_y[]   the y-ticklables
 * *matrix[]  element x,y is matrix[x][y]
 * nmap       number of color levels for the output(?)
 */
 
 !!mdatom.h
 #include "simple.h"
 #define  NO_TF_TABLE 255
 #define  DEFAULT_TF_TABLE 0

typedef struct {
    real                   tmassA, tmassB, tmass;
    int                    nr;
    int                    nalloc;
    int                    nenergrp;
    gmx_bool               bVCMgrps;
    int                    nPerturbed;
    int                    nMassPerturbed;
    int                    nChargePerturbed;
    gmx_bool               bOrires;
    real                  *massA, *massB, *massT, *invmass;
    real                  *chargeA, *chargeB;
    gmx_bool              *bPerturbed;
    int                   *typeA, *typeB;
    unsigned short        *ptype;
    unsigned short        *cTC, *cENER, *cACC, *cFREEZE, *cVCM;
    unsigned short        *cU1, *cU2, *cORF;
    /* for QMMM, atomnumber contains atomic number of the atoms */
    gmx_bool              *bQM;
    /* The range of home atoms */
    int                    start;
    int                    homenr;
    /* The lambda value used to create the contents of the struct */
    real                   lambda;
    /* The AdResS weighting function */
    real                  *wf;
    unsigned short        *tf_table_index; /* The tf table that will be applied, if thermodyn, force enabled*/
} t_mdatoms;

!!membedt.h
/* abstract data type for membed variables needed in do_md */
typedef struct membed *gmx_membed_t;

!!nb_verlet.h

#include "nbnxn_pairlist.h"
#include "nbnxn_cuda_types_ext.h"

/* For testing the reference plain-C SIMD kernels, uncomment the next lines,
 * as well as the GMX_SIMD_REFERENCE_PLAIN_C define in gmx_simd_macros.h
 * The actual SIMD width is set in gmx_simd_macros.h
 * The 4xN reference kernels support 2-, 4- and 8-way SIMD.
 * The 2x(N+N) reference kernels support 8- and 16-way SIMD.
 */
/* #define GMX_NBNXN_SIMD */
/* #define GMX_NBNXN_SIMD_4XN */
/* #define GMX_NBNXN_SIMD_2XNN */


#ifdef GMX_X86_SSE2
/* Use SIMD accelerated nbnxn search and kernels */
#define GMX_NBNXN_SIMD

/* Uncomment the next line to use, slower, 128-bit SIMD with AVX-256 */
/* #define GMX_NBNXN_HALF_WIDTH_SIMD */

/* The nbnxn SIMD 4xN and 2x(N+N) kernels can be added independently.
 * Currently the 2xNN SIMD kernels only make sense with:
 *  8-way SIMD: 4x4 setup, works with AVX-256 in single precision
 * 16-way SIMD: 4x8 setup, not used, but most of the kernel code is there
 */
#define GMX_NBNXN_SIMD_4XN
#if defined GMX_X86_AVX_256 && !(defined GMX_DOUBLE || defined GMX_NBNXN_HALF_WIDTH_SIMD)
#define GMX_NBNXN_SIMD_2XNN
#endif

#endif


/*! Nonbonded NxN kernel types: plain C, CPU SIMD, GPU CUDA, GPU emulation */
typedef enum
{
    nbnxnkNotSet = 0,
    nbnxnk4x4_PlainC,
    nbnxnk4xN_SIMD_4xN,
    nbnxnk4xN_SIMD_2xNN,
    nbnxnk8x8x8_CUDA,
    nbnxnk8x8x8_PlainC,
    nbnxnkNR
} nbnxn_kernel_type;

/*! Return a string indentifying the kernel type */
const char *lookup_nbnxn_kernel_name(int kernel_type);

enum {
    ewaldexclTable, ewaldexclAnalytical
};

/* Atom locality indicator: local, non-local, all, used for calls to:
   gridding, pair-search, force calculation, x/f buffer operations */
enum {
    eatLocal = 0, eatNonlocal = 1, eatAll
};

#define LOCAL_A(x)               ((x) == eatLocal)
#define NONLOCAL_A(x)            ((x) == eatNonlocal)
#define LOCAL_OR_NONLOCAL_A(x)   (LOCAL_A(x) || NONLOCAL_A(x))

/* Interaction locality indicator (used in pair-list search/calculations):
    - local interactions require local atom data and affect local output only;
    - non-local interactions require both local and non-local atom data and
      affect both local- and non-local output. */
enum {
    eintLocal = 0, eintNonlocal = 1
};

#define LOCAL_I(x)               ((x) == eintLocal)
#define NONLOCAL_I(x)            ((x) == eintNonlocal)

enum {
    enbvClearFNo, enbvClearFYes
};

typedef struct {
    nbnxn_pairlist_set_t  nbl_lists;   /* pair list(s)                       */
    nbnxn_atomdata_t     *nbat;        /* atom data                          */
    int                   kernel_type; /* non-bonded kernel - see enum above */
    int                   ewald_excl;  /* Ewald exclusion - see enum above   */
} nonbonded_verlet_group_t;

/* non-bonded data structure with Verlet-type cut-off */
typedef struct {
    nbnxn_search_t           nbs;             /* n vs n atom pair searching data       */
    int                      ngrp;            /* number of interaction groups          */
    nonbonded_verlet_group_t grp[2];          /* local and non-local interaction group */

    gmx_bool                 bUseGPU;         /* TRUE when GPU acceleration is used */
    nbnxn_cuda_ptr_t         cu_nbv;          /* pointer to CUDA nb verlet data     */
    int                      min_ci_balanced; /* pair list balancing parameter
                                                 used for the 8x8x8 CUDA kernels    */
} nonbonded_verlet_t;

!!nblist.h

typedef unsigned long t_excl;

/* The maximum charge group size because of minimum size of t_excl
 * could be 32 bits.
 */
#define MAX_CHARGEGROUP_SIZE 32

/* The maximum charge group size for CG-CG nblists.
 * The excl entry in t_nblist uses blocks of this size.
 */
#define MAX_CGCGSIZE 32

typedef struct
{
    int             igeometry;    /* The type of list (atom, water, etc.)  */
    int             ielec;        /* Coulomb loop type index for kernels   */
    int             ielecmod;     /* Coulomb modifier (e.g. switch/shift)  */
    int             ivdw;         /* VdW loop type index for kernels       */
    int             ivdwmod;      /* VdW modifier (e.g. switch/shift)      */
    int             type;         /* Type of interaction, listed in
                                     gmx_nblist_interaction_type           */

    int             nri, maxnri;  /* Current/max number of i particles	   */
    int             nrj, maxnrj;  /* Current/max number of j particles	   */
    int             maxlen;       /* maxnr of j atoms for a single i atom  */
    int *           iinr;         /* The i-elements                        */
    int *           iinr_end;     /* The end atom, only with enlistCG      */
    int *           gid;          /* Index in energy arrays                */
    int *           shift;        /* Shift vector index                    */
    int *           jindex;       /* Index in jjnr                         */
    int *           jjnr;         /* The j-atom list                       */
    int *           jjnr_end;     /* The end atom, only with enltypeCG     */
    t_excl *        excl;         /* Exclusions, only with enltypeCG       */

    /* We use separate pointers for kernels that compute both potential
     * and force (vf suffix), only potential (v) or only force (f)
     */
    void *          kernelptr_vf;
    void *          kernelptr_v;
    void *          kernelptr_f;

    /* Pad the list of neighbors for each i atom with "-1" entries up to the
     * simd_padding_width, if it is larger than 0. This is necessary for many
     * accelerated kernels using single-instruction multiple-data operations
     * internally.
     */
    int             simd_padding_width;

} t_nblist;


/* For atom I =  nblist->iinr[N] (0 <= N < nblist->nri) there can be
 * several neighborlists (N's), for different energy groups (gid) and
 * different shifts (shift).
 * For corresponding J atoms for each list start at:
 * nblist->jjnr[JI]
 * with nblist->jindex[N] <= JI < nblist->jindex[N+1]
 *
 * enlist is of the form enlistUNIT1_UNIT2:
 * UNIT ATOM:  there is one atom: iinr[N] or jjnr[JI]
 * UNIT SPC:   there are 3 atoms: iinr[N],iinr[N]+1,iinr[N]+2, jjnr analog.
 * UNIT TIP4P: there are 4 atoms: iinr[N],...,iinr[N]+3, jjnr analog.
 * UNIT CG:    there are N atoms: iinr[N],...,iinr_end[N]-1, jjnr analog.
 *
 * Clear?
 */
 
!!nbnxn_cuda_types_ext.h
 /* Abstract types */
/* CUDA nonbonded structure */
typedef struct nbnxn_cuda *nbnxn_cuda_ptr_t;
/* CUDA GPU device info */
typedef struct cuda_dev_info *cuda_dev_info_ptr_t;

/* Types defined for the structs below. */
typedef struct wallclock_gpu wallclock_gpu_t;
typedef struct nbnxn_cuda_ktime nbnxn_cuda_ktime_t;

/* Nonbonded kernel time and call count. */
struct nbnxn_cuda_ktime
{
    double  t;
    int     c;
};

/* GPU timings for kernels and H2d/D2H transfers. */
struct wallclock_gpu
{
    nbnxn_cuda_ktime_t ktime[2][2]; /* table containing the timings of the four
                                       version of the nonbonded kernels: force-only,
                                       force+energy, force+pruning, and force+energy+pruning */
    double  nb_h2d_t;               /* host to device transfer time in nb calculation  */
    double  nb_d2h_t;               /* device to host transfer time in nb calculation */
    int     nb_c;                   /* total call count of the nonbonded gpu operations */
    double  pl_h2d_t;               /* pair search step host to device transfer time */
    int     pl_h2d_c;               /* pair search step  host to device transfer call count */
};

!!nbnxn_pairlist.h

/* A buffer data structure of 64 bytes
 * to be placed at the beginning and end of structs
 * to avoid cache invalidation of the real contents
 * of the struct by writes to neighboring memory.
 */
typedef struct {
    int dummy[16];
} gmx_cache_protect_t;

/* Abstract type for pair searching data */
typedef struct nbnxn_search * nbnxn_search_t;

/* Function that should return a pointer *ptr to memory
 * of size nbytes.
 * Error handling should be done within this function.
 */
typedef void nbnxn_alloc_t (void **ptr, size_t nbytes);

/* Function that should free the memory pointed to by *ptr.
 * NULL should not be passed to this function.
 */
typedef void nbnxn_free_t (void *ptr);

/* This is the actual cluster-pair list j-entry.
 * cj is the j-cluster.
 * The interaction bits in excl are indexed i-major, j-minor.
 * The cj entries are sorted such that ones with exclusions come first.
 * This means that once a full mask (=NBNXN_INTERACTION_MASK_ALL)
 * is found, all subsequent j-entries in the i-entry also have full masks.
 */
typedef struct {
    int      cj;    /* The j-cluster                             */
    unsigned excl;  /* The topology exclusion (interaction) bits */
} nbnxn_cj_t;

/* In nbnxn_ci_t the integer shift contains the shift in the lower 7 bits.
 * The upper bits contain information for non-bonded kernel optimization.
 * Simply calculating LJ and Coulomb for all pairs in a cluster pair is fine.
 * But three flags can be used to skip interactions, currently only for subc=0
 * !(shift & NBNXN_CI_DO_LJ(subc))   => we can skip LJ for all pairs
 * shift & NBNXN_CI_HALF_LJ(subc)    => we can skip LJ for the second half of i
 * !(shift & NBNXN_CI_DO_COUL(subc)) => we can skip Coulomb for all pairs
 */
#define NBNXN_CI_SHIFT          127
#define NBNXN_CI_DO_LJ(subc)    (1<<(7+3*(subc)))
#define NBNXN_CI_HALF_LJ(subc)  (1<<(8+3*(subc)))
#define NBNXN_CI_DO_COUL(subc)  (1<<(9+3*(subc)))

/* Simple pair-list i-unit */
typedef struct {
    int ci;             /* i-cluster             */
    int shift;          /* Shift vector index plus possible flags, see above */
    int cj_ind_start;   /* Start index into cj   */
    int cj_ind_end;     /* End index into cj     */
} nbnxn_ci_t;

/* Grouped pair-list i-unit */
typedef struct {
    int sci;            /* i-super-cluster       */
    int shift;          /* Shift vector index plus possible flags */
    int cj4_ind_start;  /* Start index into cj4  */
    int cj4_ind_end;    /* End index into cj4    */
} nbnxn_sci_t;

typedef struct {
    unsigned imask;        /* The i-cluster interactions mask for 1 warp  */
    int      excl_ind;     /* Index into the exclusion array for 1 warp   */
} nbnxn_im_ei_t;

typedef struct {
    int           cj[4];   /* The 4 j-clusters                            */
    nbnxn_im_ei_t imei[2]; /* The i-cluster mask data       for 2 warps   */
} nbnxn_cj4_t;

typedef struct {
    unsigned pair[32];     /* Topology exclusion interaction bits for one warp,
                            * each unsigned has bitS for 4*8 i clusters
                            */
} nbnxn_excl_t;

typedef struct {
    gmx_cache_protect_t cp0;

    nbnxn_alloc_t      *alloc;
    nbnxn_free_t       *free;

    gmx_bool            bSimple;         /* Simple list has na_sc=na_s and uses cj   *
                                          * Complex list uses cj4                    */

    int                     na_ci;       /* The number of atoms per i-cluster        */
    int                     na_cj;       /* The number of atoms per j-cluster        */
    int                     na_sc;       /* The number of atoms per super cluster    */
    real                    rlist;       /* The radius for constructing the list     */
    int                     nci;         /* The number of i-clusters in the list     */
    nbnxn_ci_t             *ci;          /* The i-cluster list, size nci             */
    int                     ci_nalloc;   /* The allocation size of ci                */
    int                     nsci;        /* The number of i-super-clusters in the list */
    nbnxn_sci_t            *sci;         /* The i-super-cluster list                 */
    int                     sci_nalloc;  /* The allocation size of sci               */

    int                     ncj;         /* The number of j-clusters in the list     */
    nbnxn_cj_t             *cj;          /* The j-cluster list, size ncj             */
    int                     cj_nalloc;   /* The allocation size of cj                */

    int                     ncj4;        /* The total number of 4*j clusters         */
    nbnxn_cj4_t            *cj4;         /* The 4*j cluster list, size ncj4          */
    int                     cj4_nalloc;  /* The allocation size of cj4               */
    int                     nexcl;       /* The count for excl                       */
    nbnxn_excl_t           *excl;        /* Atom interaction bits (non-exclusions)   */
    int                     excl_nalloc; /* The allocation size for excl             */
    int                     nci_tot;     /* The total number of i clusters           */

    struct nbnxn_list_work *work;

    gmx_cache_protect_t     cp1;
} nbnxn_pairlist_t;

typedef struct {
    int                nnbl;        /* number of lists */
    nbnxn_pairlist_t **nbl;         /* lists */
    gmx_bool           bCombined;   /* TRUE if lists get combined into one (the 1st) */
    gmx_bool           bSimple;     /* TRUE if the list of of type "simple"
                                       (na_sc=na_s, no super-clusters used) */
    int                natpair_ljq; /* Total number of atom pairs for LJ+Q kernel */
    int                natpair_lj;  /* Total number of atom pairs for LJ kernel   */
    int                natpair_q;   /* Total number of atom pairs for Q kernel    */
} nbnxn_pairlist_set_t;

enum {
    nbatXYZ, nbatXYZQ, nbatX4, nbatX8
};

typedef struct {
    real *f;      /* f, size natoms*fstride                             */
    real *fshift; /* Shift force array, size SHIFTS*DIM                 */
    int   nV;     /* The size of *Vvdw and *Vc                          */
    real *Vvdw;   /* Temporary Van der Waals group energy storage       */
    real *Vc;     /* Temporary Coulomb group energy storage             */
    int   nVS;    /* The size of *VSvdw and *VSc                        */
    real *VSvdw;  /* Temporary SIMD Van der Waals group energy storage  */
    real *VSc;    /* Temporary SIMD Coulomb group energy storage        */
} nbnxn_atomdata_output_t;

/* Block size in atoms for the non-bonded thread force-buffer reduction,
 * should be a multiple of all cell and x86 SIMD sizes (i.e. 2, 4 and 8).
 * Should be small to reduce the reduction and zeroing cost,
 * but too small will result in overhead.
 * Currently the block size is NBNXN_BUFFERFLAG_SIZE*3*sizeof(real)=192 bytes.
 */
#ifdef GMX_DOUBLE
#define NBNXN_BUFFERFLAG_SIZE   8
#else
#define NBNXN_BUFFERFLAG_SIZE  16
#endif

/* We currently store the reduction flags as bits in an unsigned int.
 * In most cases this limits the number of flags to 32.
 * The reduction will automatically disable the flagging and do a full
 * reduction when the flags won't fit, but this will lead to very slow
 * reduction. As we anyhow don't expect reasonable performance with
 * more than 32 threads, we put in this hard limit.
 * You can increase this number, but the reduction will be very slow.
 */
#define NBNXN_BUFFERFLAG_MAX_THREADS  32

/* Flags for telling if threads write to force output buffers */
typedef struct {
    int       nflag;       /* The number of flag blocks                         */
    unsigned *flag;        /* Bit i is set when thread i writes to a cell-block */
    int       flag_nalloc; /* Allocation size of cxy_flag                       */
} nbnxn_buffer_flags_t;

/* LJ combination rules: geometric, Lorentz-Berthelot, none */
enum {
    ljcrGEOM, ljcrLB, ljcrNONE, ljcrNR
};

typedef struct {
    nbnxn_alloc_t           *alloc;
    nbnxn_free_t            *free;
    int                      ntype;           /* The number of different atom types                 */
    real                    *nbfp;            /* Lennard-Jones 6*C6 and 12*C12 params, size ntype^2*2 */
    int                      comb_rule;       /* Combination rule, see enum above                   */
    real                    *nbfp_comb;       /* LJ parameter per atom type, size ntype*2           */
    real                    *nbfp_s4;         /* As nbfp, but with stride 4, size ntype^2*4. This
                                               * might suit 4-wide SIMD loads of two values (e.g.
                                               * two floats in single precision on x86).            */
    int                      natoms;          /* Number of atoms                                    */
    int                      natoms_local;    /* Number of local atoms                           */
    int                     *type;            /* Atom types                                         */
    real                    *lj_comb;         /* LJ parameters per atom for combining for pairs     */
    int                      XFormat;         /* The format of x (and q), enum                      */
    int                      FFormat;         /* The format of f, enum                              */
    real                    *q;               /* Charges, can be NULL if incorporated in x          */
    int                      na_c;            /* The number of atoms per cluster                    */
    int                      nenergrp;        /* The number of energy groups                        */
    int                      neg_2log;        /* Log2 of nenergrp                                   */
    int                     *energrp;         /* The energy groups per cluster, can be NULL         */
    gmx_bool                 bDynamicBox;     /* Do we need to update shift_vec every step?    */
    rvec                    *shift_vec;       /* Shift vectors, copied from t_forcerec              */
    int                      xstride;         /* stride for a coordinate in x (usually 3 or 4)      */
    int                      fstride;         /* stride for a coordinate in f (usually 3 or 4)      */
    real                    *x;               /* x and possibly q, size natoms*xstride              */

    /* j-atom minus i-atom index for generating self and Newton exclusions
     * cluster-cluster pairs of the diagonal, for 4xn and 2xnn kernels.
     */
    real                    *simd_4xn_diagonal_j_minus_i;
    real                    *simd_2xnn_diagonal_j_minus_i;
    /* Filters for topology exclusion masks for the SIMD kernels.
     * filter2 is the same as filter1, but with each element duplicated.
     */
    unsigned                *simd_exclusion_filter1;
    unsigned                *simd_exclusion_filter2;

    int                      nout;            /* The number of force arrays                         */
    nbnxn_atomdata_output_t *out;             /* Output data structures               */
    int                      nalloc;          /* Allocation size of all arrays (for x/f *x/fstride) */
    gmx_bool                 bUseBufferFlags; /* Use the flags or operate on all atoms     */
    nbnxn_buffer_flags_t     buffer_flags;    /* Flags for buffer zeroing+reduc.  */
} nbnxn_atomdata_t;

!!nlistheuristics.h
typedef struct {
    gmx_bool        bGStatEveryStep;
    gmx_large_int_t step_ns;
    gmx_large_int_t step_nscheck;
    gmx_large_int_t nns;
    matrix          scale_tot;
    int             nabnsb;
    double          s1;
    double          s2;
    double          ab;
    double          lt_runav;
    double          lt_runav2;
} gmx_nlheur_t;

void reset_nlistheuristics(gmx_nlheur_t *nlh, gmx_large_int_t step);

void init_nlistheuristics(gmx_nlheur_t *nlh,
                          gmx_bool bGStatEveryStep, gmx_large_int_t step);

void update_nliststatistics(gmx_nlheur_t *nlh, gmx_large_int_t step);

void set_nlistheuristics(gmx_nlheur_t *nlh, gmx_bool bReset, gmx_large_int_t step);

!!nrnb.h
#define eNR_NBKERNEL_NONE -1

enum
{
    eNR_NBKERNEL_VDW_VF,
    eNR_NBKERNEL_VDW_F,
    eNR_NBKERNEL_ELEC_VF,
    eNR_NBKERNEL_ELEC_F,
    eNR_NBKERNEL_ELEC_W3_VF,
    eNR_NBKERNEL_ELEC_W3_F,
    eNR_NBKERNEL_ELEC_W3W3_VF,
    eNR_NBKERNEL_ELEC_W3W3_F,
    eNR_NBKERNEL_ELEC_W4_VF,
    eNR_NBKERNEL_ELEC_W4_F,
    eNR_NBKERNEL_ELEC_W4W4_VF,
    eNR_NBKERNEL_ELEC_W4W4_F,
    eNR_NBKERNEL_ELEC_VDW_VF,
    eNR_NBKERNEL_ELEC_VDW_F,
    eNR_NBKERNEL_ELEC_VDW_W3_VF,
    eNR_NBKERNEL_ELEC_VDW_W3_F,
    eNR_NBKERNEL_ELEC_VDW_W3W3_VF,
    eNR_NBKERNEL_ELEC_VDW_W3W3_F,
    eNR_NBKERNEL_ELEC_VDW_W4_VF,
    eNR_NBKERNEL_ELEC_VDW_W4_F,
    eNR_NBKERNEL_ELEC_VDW_W4W4_VF,
    eNR_NBKERNEL_ELEC_VDW_W4W4_F,

    eNR_NBKERNEL_NR,                        /* Total number of interaction-specific kernel entries */

    eNR_NBKERNEL_GENERIC = eNR_NBKERNEL_NR, /* Reuse number; KERNEL_NR is not an entry itself */
    eNR_NBKERNEL_GENERIC_CG,
    eNR_NBKERNEL_GENERIC_ADRESS,
    eNR_NBKERNEL_FREE_ENERGY,               /* Add other generic kernels _before_ the free energy one */

    eNR_NBKERNEL_ALLVSALL,
    eNR_NBKERNEL_ALLVSALLGB,

    eNR_NBNXN_DIST2,
    eNR_NBNXN_LJ_RF,    eNR_NBNXN_LJ_RF_E,
    eNR_NBNXN_LJ_TAB,   eNR_NBNXN_LJ_TAB_E,
    eNR_NBNXN_LJ_EWALD, eNR_NBNXN_LJ_EWALD_E,
    eNR_NBNXN_LJ,       eNR_NBNXN_LJ_E,
    eNR_NBNXN_RF,       eNR_NBNXN_RF_E,
    eNR_NBNXN_TAB,      eNR_NBNXN_TAB_E,
    eNR_NBNXN_EWALD,    eNR_NBNXN_EWALD_E,
    eNR_NB14,
    eNR_BORN_RADII_STILL,     eNR_BORN_RADII_HCT_OBC,
    eNR_BORN_CHAINRULE,
    eNR_BORN_AVA_RADII_STILL, eNR_BORN_AVA_RADII_HCT_OBC,
    eNR_BORN_AVA_CHAINRULE,
    eNR_WEIGHTS,              eNR_SPREADQ,              eNR_SPREADQBSP,
    eNR_GATHERF,              eNR_GATHERFBSP,           eNR_FFT,
    eNR_CONV,                 eNR_SOLVEPME, eNR_NS,      eNR_RESETX,
    eNR_SHIFTX,               eNR_CGCM,                 eNR_FSUM,
    eNR_BONDS,                eNR_G96BONDS,             eNR_FENEBONDS,
    eNR_TABBONDS,             eNR_RESTRBONDS,           eNR_LINEAR_ANGLES,
    eNR_ANGLES,               eNR_G96ANGLES,            eNR_QANGLES,
    eNR_TABANGLES,            eNR_PROPER,               eNR_IMPROPER,
    eNR_RB,                   eNR_FOURDIH,              eNR_TABDIHS,
    eNR_DISRES,               eNR_ORIRES,               eNR_DIHRES,
    eNR_POSRES,               eNR_FBPOSRES,
    eNR_ANGRES,               eNR_ANGRESZ,
    eNR_MORSE,                eNR_CUBICBONDS,           eNR_WALLS,
    eNR_POLARIZE,             eNR_ANHARM_POL,
    eNR_WPOL,                 eNR_THOLE,                eNR_VIRIAL,
    eNR_UPDATE,               eNR_EXTUPDATE,            eNR_STOPCM,
    eNR_PCOUPL,               eNR_EKIN,                 eNR_LINCS,
    eNR_LINCSMAT,             eNR_SHAKE,                eNR_CONSTR_V,
    eNR_SHAKE_RIJ,            eNR_CONSTR_VIR,           eNR_SETTLE,
    eNR_VSITE2,               eNR_VSITE3,               eNR_VSITE3FD,
    eNR_VSITE3FAD,            eNR_VSITE3OUT,            eNR_VSITE4FD,
    eNR_VSITE4FDN,            eNR_VSITEN,               eNR_GB,
    eNR_CMAP,
    eNRNB
};


typedef struct
{
    double n[eNRNB];
}
t_nrnb;


typedef struct gmx_wallcycle *gmx_wallcycle_t;

!!ns.h

#include "nsgrid.h"
#include "nblist.h"

enum {
    eNL_VDWQQ, eNL_VDW, eNL_QQ,
    eNL_VDWQQ_FREE, eNL_VDW_FREE, eNL_QQ_FREE,
    eNL_VDWQQ_WATER, eNL_QQ_WATER,
    eNL_VDWQQ_WATERWATER, eNL_QQ_WATERWATER,
    eNL_NR
};

#define MAX_CG 1024

typedef struct {
    int     ncg;
    int     nj;
    atom_id jcg[MAX_CG];
} t_ns_buf;

typedef struct {
    gmx_bool      bCGlist;
    atom_id      *simple_aaj;
    t_grid       *grid;
    t_excl       *bexcl;
    gmx_bool     *bHaveVdW;
    t_ns_buf    **ns_buf;
    gmx_bool     *bExcludeAlleg;
    int           nra_alloc;
    int           cg_alloc;
    atom_id     **nl_sr;
    int          *nsr;
    atom_id     **nl_lr_ljc;
    atom_id     **nl_lr_one;
    int          *nlr_ljc;
    int          *nlr_one;
    /* the nblists should probably go in here */
    gmx_bool      nblist_initialized; /* has the nblist been initialized?  */
    int           dump_nl;            /* neighbour list dump level (from env. var. GMX_DUMP_NL)*/
} gmx_ns_t;

!!nsgrid.h

#include "simple.h"

typedef struct {
    int     nr;           /* Total number of charge groups	*/
    int     nboundeddim;  /* The number of bounded dimensions     */
    int     npbcdim;      /* The number of dimensions with pbc    */
    int     ncg_ideal;    /* The ideal number of cg's per cell    */
    ivec    n;            /* The dimension of the grid		*/
    int     ncells;       /* Total number of cells		*/
    int     cells_nalloc; /* Allocation size of index and nra       */
    ivec    ncpddc;       /* The number of cells per DD cell      */
    rvec    cell_size;    /* The size of the cells                */
    rvec    cell_offset;  /* The offset of the cell (0,0,0)       */
    int    *cell_index;   /* The cell number of each cg		*/
    int    *index;        /* The index into a for each cell	*/
    /* The location of the cell in the index*/
    /* array can be found by calling xyz2ci	*/
    int    *nra;    /* The number of entries in a cell	*/
    int     icg0;   /* The start of the i-cg range          */
    int     icg1;   /* The end of the i-cg range            */
    rvec   *os0;
    rvec   *os1;
    int    *a;         /* The grid of cgs			*/
    int     nr_alloc;  /* Allocation size of cell_index and a  */
    real   *dcx2;      /* Squared distance from atom to j-cell */
    real   *dcy2;      /* Squared distance from atom to j-cell */
    real   *dcz2;      /* Squared distance from atom to j-cell */
    int     dc_nalloc; /* Allocation size of dcx2, dyc2, dcz2  */
} t_grid;

!!oenv.h

/* output options opaque type, for functions in statutil.h and oenv.h */
typedef struct output_env *output_env_t;

!!pbc.h

#include "simple.h"
#define MAX_NTRICVEC 12

typedef struct {
    int        ePBC;
    int        ndim_ePBC;
    int        ePBCDX;
    int        dim;
    matrix     box;
    rvec       fbox_diag;
    rvec       hbox_diag;
    rvec       mhbox_diag;
    real       max_cutoff2;
    gmx_bool   bLimitDistance;
    real       limit_distance2;
    int        ntric_vec;
    ivec       tric_shift[MAX_NTRICVEC];
    rvec       tric_vec[MAX_NTRICVEC];
} t_pbc;

!!qmmmrec.h

#include "simple.h"

typedef struct {
    int                nrQMatoms;      /* total nr of QM atoms              */
    rvec              *xQM;            /* shifted to center of box          */
    int               *indexQM;        /* atom i = atom indexQM[i] in mdrun */
    int               *atomicnumberQM; /* atomic numbers of QM atoms        */
    real              *QMcharges;      /* atomic charges of QM atoms(ONIOM) */
    int               *shiftQM;
    int                QMcharge;       /* charge of the QM system           */
    int                multiplicity;   /* multipicity (no of unpaired eln)  */
    int                QMmethod;       /* see enums.h for all methods       */
    int                QMbasis;        /* see enums.h for all bases         */
    int                nelectrons;     /* total number of elecs in QM region*/
    gmx_bool           bTS;            /* Optimize a TS, only steep, no md  */
    gmx_bool           bOPT;           /* Optimize QM subsys, only steep, no md  */
    gmx_bool          *frontatoms;     /* qm atoms on the QM side of a QM-MM bond */
    /* Gaussian specific stuff */
    int                nQMcpus;        /* no. of CPUs used for the QM calc. */
    int                QMmem;          /* memory for the gaussian calc.     */
    int                accuracy;       /* convergence criterium (E(-x))     */
    gmx_bool           cpmcscf;        /* using cpmcscf(l1003)*/
    char              *gauss_dir;
    char              *gauss_exe;
    char              *devel_dir;
    char              *orca_basename; /* basename for I/O with orca        */
    char              *orca_dir;      /* directory for ORCA                */
    real              *c6;
    real              *c12;
    /* Surface hopping stuff */
    gmx_bool           bSH;     /* surface hopping (diabatic only)   */
    real               SAon;    /* at which energy gap the SA starts */
    real               SAoff;   /* at which energy gap the SA stops  */
    int                SAsteps; /* stepwise switchinng on the SA     */
    int                SAstep;  /* current state of SA               */
    int                CIdim;
    real              *CIvec1;
    real              *CIvec2;
    real              *CIvec1old;
    real              *CIvec2old;
    ivec               SHbasis;
    int                CASelectrons;
    int                CASorbitals;
} t_QMrec;

typedef struct {
    int            nrMMatoms;   /* nr of MM atoms, updated every step*/
    rvec          *xMM;         /* shifted to center of box          */
    int           *indexMM;     /* atom i = atom indexMM[I] in mdrun */
    real          *MMcharges;   /* MM point charges in std QMMM calc.*/
    int           *shiftMM;
    int           *MMatomtype;  /* only important for semi-emp.      */
    real           scalefactor;
    /* gaussian specific stuff */
    real          *c6;
    real          *c12;
} t_MMrec;


typedef struct {
    int             QMMMscheme; /* ONIOM (multi-layer) or normal          */
    int             nrQMlayers; /* number of QM layers (total layers +1 (MM)) */
    t_QMrec       **qm;         /* atoms and run params for each QM group */
    t_MMrec        *mm;         /* there can only be one MM subsystem !   */
} t_QMMMrec;






























