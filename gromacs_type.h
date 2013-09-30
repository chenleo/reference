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






































