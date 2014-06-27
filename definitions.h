/* 
** definitions.h
*/

#define NINTMIN 5
#define NINTMAX 28
#define NINTMINSIGMA 5
#define NINTMAXSIGMA 28
#define NINTMINDF 8
#define NINTMAXDF 28
#define TOL 1e-10
#define TOLSIGMA 1e-10
#define TOLDF 1e-3
#define TOLLININT 1e-10
#define DFFAILUREMAX 1e20
#define SBI 1e100
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
	DOUBLE *rhoBulge;
	DOUBLE *logrhoBulge;
	DOUBLE *rhoHalo;
	DOUBLE *logrhoHalo;
	DOUBLE *rhoenc;
	DOUBLE *logrhoenc;
	DOUBLE *rhoencBulge;
	DOUBLE *logrhoencBulge;
	DOUBLE *rhoencHalo;
	DOUBLE *logrhoencHalo;
	DOUBLE *Menc;
	DOUBLE *logMenc;
	DOUBLE *MencBulge;
	DOUBLE *logMencBulge;
	DOUBLE *MencHalo;
	DOUBLE *logMencHalo;
	DOUBLE *Pot;
	DOUBLE *logPot;
	DOUBLE *Potoutr;
	DOUBLE *sigma;
	DOUBLE *logsigma;
	DOUBLE *sigmaintoutr;
	DOUBLE *eqrvcmax;
	DOUBLE *eqrvcmaxBulge;
	DOUBLE *eqrvcmaxHalo;
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
	DOUBLE rba_r0;
	DOUBLE rba_at_r0;
	DOUBLE rba_slope;
	DOUBLE rba_min;
	DOUBLE rba_max;
	DOUBLE rca_r0;
	DOUBLE rca_at_r0;
	DOUBLE rca_slope;
	DOUBLE rca_min;
	DOUBLE rca_max;
	DOUBLE alpha_r0;
	DOUBLE alpha_at_r0;
	DOUBLE alpha_slope;
	DOUBLE alpha_min;
	DOUBLE alpha_max;
	DOUBLE beta_r0;
	DOUBLE beta_at_r0;
	DOUBLE beta_slope;
	DOUBLE beta_min;
	DOUBLE beta_max;
	DOUBLE gamma_r0;
	DOUBLE gamma_at_r0;
	DOUBLE gamma_slope;
	DOUBLE gamma_min;
	DOUBLE gamma_max;
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

typedef struct samplinginfo {

	INT Ntot;
	INT Ninitialtot;
	INT Nnosplittot;
	INT Nnewtot;
	DOUBLE Neff;
	DOUBLE Mp;
	DOUBLE Ekin;
	DOUBLE Epot;
	DOUBLE Etot;
	DOUBLE Nfemm;
	DOUBLE Nfesm;
	DOUBLE Cr[4];
	DOUBLE Cv[4];
	DOUBLE Ltot[4];
	} SAMP;

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
	INT soft0_in_rs_units;
	INT soft0_in_rvir_units;
	INT soft0_in_rcutoff_units;
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
	DOUBLE *logrhoenc;
	DOUBLE *logMenc;
	DOUBLE *eqrvcmax;
	SP *sp;
	SHELL *shell;
	GRIDDF *griddf;
	SAMP *samp;
	CHAR systemname[STRINGSIZE];
	} SI;

typedef struct generalinfo {

	INT do_bh;
	INT do_bulge;
	INT do_halo;
	INT output_gridr;
	INT output_griddf;
	INT output_tipsy_standard;
	INT output_tipsy_standard_dpp;
	INT coordinates;
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
	DOUBLE f_rinner;
	DOUBLE f_router;
	DOUBLE f_cutoff;
	DOUBLE randomseed;
	DOUBLE rvcmax;
	DOUBLE vcmax;
	DOUBLE t[10];
	GRIDR *gridr;
	SAMP *samp;
	CHAR outputname[STRINGSIZE];
	} GI;
