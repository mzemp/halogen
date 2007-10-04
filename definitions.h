/* 
** definitions.h
**
** Some definitions for halogen
*/

#define VERSION "1.0"
#define NGRIDR 2001
#define NGRIDDF 101
#define NINTMIN 5
#define NINTMAX 28
#define NINTMINDF 8
#define NINTMAXDF 28
#define TOL 1e-10
#define TOLDF 1e-3
#define TOLLININT 1e-10
#define DFFAILUREMAX 1e20
#define FACTORRINNER 1e-6
#define FACTORROUTER 1e20
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

    DOUBLE r[NGRIDR];
    DOUBLE logr[NGRIDR];
    DOUBLE rho[NGRIDR];
    DOUBLE logrho[NGRIDR];
    DOUBLE rhoHalo[NGRIDR];
    DOUBLE logrhoHalo[NGRIDR];
    DOUBLE rhoenc[NGRIDR];
    DOUBLE logrhoenc[NGRIDR];
    DOUBLE rhoencHalo[NGRIDR];
    DOUBLE logrhoencHalo[NGRIDR];
    DOUBLE Menc[NGRIDR];
    DOUBLE logMenc[NGRIDR];
    DOUBLE MencHalo[NGRIDR];
    DOUBLE logMencHalo[NGRIDR];
    DOUBLE Pot[NGRIDR];
    DOUBLE logPot[NGRIDR];
    DOUBLE Potoutr[NGRIDR];
    DOUBLE eqrvcmax[NGRIDR];
    } GRIDR;

typedef struct griddf {

    DOUBLE r[NGRIDDF];
    DOUBLE logr[NGRIDDF];
    DOUBLE E[NGRIDDF];
    DOUBLE logE[NGRIDDF];
    DOUBLE fE[NGRIDDF];
    DOUBLE logfE[NGRIDDF];
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
    } STUFF;

typedef struct systeminfo {

    INT set_rsi_to_rs;
    INT set_rsi_to_rvir;
    INT set_rsi_to_rcutoff;
    INT set_rso_to_rs;
    INT set_rso_to_rvir;
    INT set_rso_to_rcutoff;
    INT set_rmor_to_rs;
    INT set_rmor_to_rvir;
    INT set_rmor_to_rcutoff;
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
    DOUBLE logr[NGRIDR];
    DOUBLE logMenc[NGRIDR];
    DOUBLE logrhoenc[NGRIDR];
    SP *sp;
    SHELL *shell;
    GRIDDF *griddf;
    CHAR systemname[STRINGSIZE];
    } SI;

typedef struct generalinfo {

    DOUBLE OmegaM0;
    DOUBLE OmegaK0;
    DOUBLE OmegaL0;
    DOUBLE h0;
    DOUBLE z;
    DOUBLE rhocritz;
    DOUBLE Deltavirz;
    DOUBLE rinner;
    DOUBLE router;
    STUFF *stuff;
    GRIDR *gridr;
    } GI;
