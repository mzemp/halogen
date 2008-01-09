/* 
** definitions.h
*/

#define VERSION "1.1"
#define NINTMIN 5
#define NINTMAX 28
#define NINTMINDF 8
#define NINTMAXDF 28
#define TOL 1e-10
#define TOLDF 1e-3
#define TOLLININT 1e-10
#define DFFAILUREMAX 1e20
#define SBI 1e100
#define CutoffFac 0.3
#define MU 2.2229621e5
#define G 1
#define VelConvertFac 1.0227122
#define STRINGSIZE 50
#define INT int
#define FLOAT float
#define DOUBLE double
#define CHAR char
#define OFI1 "%d"
#define OFI2 "%-3d"
#define OFI3 "%-10d"
#define OFD1 "%g"
#define OFD2 "%8.7e"
#define OFD3 "%8.7e"
#define OFD4 "%+8.7e"

typedef struct gridr {

    DOUBLE *r;
    DOUBLE *logr;
    DOUBLE *rho;
    DOUBLE *logrho;
    DOUBLE *rhoHalo;
    DOUBLE *logrhoHalo;
    DOUBLE *rhoenc;
    DOUBLE *logrhoenc;
    DOUBLE *rhoencHalo;
    DOUBLE *logrhoencHalo;
    DOUBLE *Menc;
    DOUBLE *logMenc;
    DOUBLE *MencHalo;
    DOUBLE *logMencHalo;
    DOUBLE *Pot;
    DOUBLE *logPot;
    DOUBLE *Potoutr;
    DOUBLE *eqrvcmax;
    } GRIDR;

typedef struct griddf {

    DOUBLE *r;
    DOUBLE *logr;
    DOUBLE *E;
    DOUBLE *logE;
    DOUBLE *fE;
    DOUBLE *logfE;
    } GRIDDF;

typedef struct systemparameters {

    DOUBLE alpha;
    DOUBLE beta;
    DOUBLE gamma;
    DOUBLE delta;
    DOUBLE rho0;
    DOUBLE rs;
    DOUBLE rcutoff;
    DOUBLE rdecay;
    DOUBLE M;
    DOUBLE cvir;
    DOUBLE rvir;
    DOUBLE vvir;
    DOUBLE rhalf;
    DOUBLE rvcmax;
    DOUBLE vcmax;
    } SP;
    
typedef struct particle {

    DOUBLE r[4];
    DOUBLE v[4];
    DOUBLE L[4];
    DOUBLE mass;
    DOUBLE soft;
    DOUBLE Ekin;
    DOUBLE Epot;
    DOUBLE Etot;
    DOUBLE rperi;
    DOUBLE rapo;
    DOUBLE ecc;
    } PARTICLE;

typedef struct shell {

    INT N;
    INT Ninitial;
    INT Nnosplit;
    INT Nnew;
    INT massfac;
    DOUBLE rinner;
    DOUBLE router;
    DOUBLE Mp;
    DOUBLE Menc;
    DOUBLE mass;
    DOUBLE soft;
    PARTICLE *p;
    } SHELL;

typedef struct stuff {

    INT Ntot;
    INT Ninitialtot;
    INT Nnosplittot;
    INT Nnewtot;
    DOUBLE Mp;
    DOUBLE Ekin;
    DOUBLE Epot;
    DOUBLE Etot;
    DOUBLE Cr[4];
    DOUBLE Cv[4];
    DOUBLE Ltot[4];
    DOUBLE Nfemm;
    DOUBLE Nfesm;
    } STUFF;

typedef struct systeminfo {

    INT rsi_in_rs_units;
    INT rsi_in_rvir_units;
    INT rsi_in_rcutoff_units;
    INT rso_in_rs_units;
    INT rso_in_rvir_units;
    INT rso_in_rcutoff_units;
    INT rmor_in_rs_units;
    INT rmor_in_rvir_units;
    INT rmor_in_rcutoff_units;
    INT Nshell;
    INT Ismor;
    INT N0;
    INT DRMmax;
    DOUBLE soft0;
    DOUBLE rsi;
    DOUBLE rso;
    DOUBLE rmor;
    DOUBLE rimp;
    DOUBLE r1;
    DOUBLE r100;
    DOUBLE dfsf;
    SP *sp;
    SHELL *shell;
    GRIDDF *griddf;
    CHAR systemname[STRINGSIZE];
    } SI;

typedef struct generalinfo {

    INT output_gridr;
    INT output_griddf;
    INT output_tipsy_standard;
    INT output_tipsy_standard_dpp;
    INT output_gadget_binary;
    INT positionsonly;
    INT Ngridr;
    INT Ngriddf;
    DOUBLE OmegaM0;
    DOUBLE OmegaK0;
    DOUBLE OmegaL0;
    DOUBLE h0;
    DOUBLE z;
    DOUBLE rhocritz;
    DOUBLE Deltavirz;
    DOUBLE OmegaMz;
    DOUBLE rinner;
    DOUBLE router;
    DOUBLE factor_rinner;
    DOUBLE factor_router;
    DOUBLE randomseed;
    DOUBLE t[10];
    STUFF *stuff;
    GRIDR *gridr;
    CHAR inputname[STRINGSIZE];
    } GI;
