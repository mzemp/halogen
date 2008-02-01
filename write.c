/*
** write.c
**
** Write routines for HALOGEN
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <time.h>
#include <iof.h>
#include "definitions.h"
#include "functions.h"
#include "routines.h"
#include "write.h"

/*
** Routines for writing out grids
*/

void write_gridr_total(FILE *file, const GI *gi) {

    INT i;
    GRIDR *gridr;

    gridr = gi->gridr;
    for (i = 0; i < gi->Ngridr; i++) {
	fprintf(file,OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2"\n",
		gridr->r[i],gridr->logr[i],
		gridr->rho[i],gridr->logrho[i],
		gridr->rhoenc[i],gridr->logrhoenc[i],
		gridr->Menc[i],gridr->logMenc[i],
		gridr->sigma[i],gridr->logsigma[i],
		gridr->Pot[i],gridr->logPot[i]);
	}
    }

void write_gridr_system(FILE *file, const GI *gi, const SI *si) {

    INT i;
    GRIDR *gridr;

    gridr = gi->gridr;
    if (strcmp(si->systemname,"bulge") == 0) {
	for (i = 0; i < gi->Ngridr; i++) {
	    fprintf(file,OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2"\n",
		    gridr->r[i],gridr->logr[i],
		    gridr->rhoBulge[i],gridr->logrhoBulge[i],
		    gridr->rhoencBulge[i],gridr->logrhoencBulge[i],
		    gridr->MencBulge[i],gridr->logMencBulge[i]);
	    }
	}
    else if (strcmp(si->systemname,"halo") == 0) {
	for (i = 0; i < gi->Ngridr; i++) {
	    fprintf(file,OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2"\n",
		    gridr->r[i],gridr->logr[i],
		    gridr->rhoHalo[i],gridr->logrhoHalo[i],
		    gridr->rhoencHalo[i],gridr->logrhoencHalo[i],
		    gridr->MencHalo[i],gridr->logMencHalo[i]);
	    }
	}
    }

void write_griddf_system(FILE *file, const GI *gi, const SI *si) {

    INT i;
    GRIDDF *griddf;

    griddf = si->griddf;
    for (i = 0; i < gi->Ngriddf; i++) {
	fprintf(file,OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2"\n",
		griddf->r[i],griddf->logr[i],
		griddf->E[i],griddf->logE[i],
		griddf->fE[i],griddf->logfE[i]);
	}
    }

/*
** Routine for writing tipsy standard format memory efficient
*/

void write_tipsy_standard_halogen(FILE *fp, const GI *gi, const PARTICLE *bh, const SI *bulge, const SI *halo) {

    INT i, j, k, Ntotal, Ndark, Nstar;
    TIPSY_HEADER th;
    DARK_PARTICLE dp;
    STAR_PARTICLE sp;
    PARTICLE *p;
    XDR xdrs;

    xdrstdio_create(&xdrs,fp,XDR_ENCODE);
    /*
    ** Write out header
    */
    Ntotal = 0;
    Ndark = 0;
    Nstar = 0;
    if (bh->mass != 0) {
	Ndark++;
	}
    if (gi->do_bulge == 1) {
	for (j = 0; j < (bulge->Nshell+2); j++) {
	    Nstar +=  bulge->shell[j].N;
	    }
	}
    if (gi->do_halo == 1) {
	for (j = 0; j < (halo->Nshell+2); j++) {
	    Ndark +=  halo->shell[j].N;
	    }
	}
    Ntotal = Ndark + Nstar;
    th.time = 0;
    th.ntotal = Ntotal;
    th.ndim = 3;
    th.ngas = 0;
    th.ndark = Ndark;
    th.nstar = Nstar;
    write_tipsy_standard_header(&xdrs,&th);
    /*
    ** Write out dark matter particles
    */
    if (bh->mass != 0) {
	for (k = 0; k < 3; k++) {
	    dp.pos[k] = bh->r[k+1];
	    dp.vel[k] = bh->v[k+1];
	    }
	dp.mass = bh->mass;
	dp.eps = bh->soft;
	dp.phi = bh->Epot;
	write_tipsy_standard_dark(&xdrs,&dp);
	}
    if (gi->do_halo == 1) {
	for (j = 0; j < (halo->Nshell+2); j++) {
	    p = halo->shell[j].p;
	    for (i = 0; i < halo->shell[j].N; i++) {
		for (k = 0; k < 3; k++) {
		    dp.pos[k] = p[i].r[k+1];
		    dp.vel[k] = p[i].v[k+1];
		    }
		dp.mass = p[i].mass;
		dp.eps = p[i].soft;
		dp.phi = p[i].Epot;
		write_tipsy_standard_dark(&xdrs,&dp);
		}
	    }
	}
    if (gi->do_bulge == 1) {
	for (j = 0; j < (bulge->Nshell+2); j++) {
	    p = bulge->shell[j].p;
	    for (i = 0; i < bulge->shell[j].N; i++) {
		for (k = 0; k < 3; k++) {
		    sp.pos[k] = p[i].r[k+1];
		    sp.vel[k] = p[i].v[k+1];
		    }
		sp.mass = p[i].mass;
		sp.eps = p[i].soft;
		sp.phi = p[i].Epot;
		sp.metals = 0;
		sp.tform = 0;
		write_tipsy_standard_star(&xdrs,&sp);
		}
	    }
	}
    xdr_destroy(&xdrs);
    }

/*
** Routine for writing tipsy standard format with double precision positions memory efficient
*/

void write_tipsy_standard_dpp_halogen(FILE *fp, const GI *gi, const PARTICLE *bh, const SI *bulge, const SI *halo) {

    INT i, j, k, Ntotal, Ndark, Nstar;
    TIPSY_HEADER th;
    DARK_PARTICLE_DPP dpdpp;
    STAR_PARTICLE_DPP spdpp;
    PARTICLE *p;
    XDR xdrs;

    xdrstdio_create(&xdrs,fp,XDR_ENCODE);
    /*
    ** Write out header
    */
    Ntotal = 0;
    Ndark = 0;
    Nstar = 0;
    if (bh->mass != 0) {
	Ndark++;
	}
    if (gi->do_bulge == 1) {
	for (j = 0; j < (bulge->Nshell+2); j++) {
	    Nstar +=  bulge->shell[j].N;
	    }
	}
    if (gi->do_halo == 1) {
	for (j = 0; j < (halo->Nshell+2); j++) {
	    Ndark +=  halo->shell[j].N;
	    }
	}
    Ntotal = Ndark + Nstar;
    th.time = 0;
    th.ntotal = Ntotal;
    th.ndim = 3;
    th.ngas = 0;
    th.ndark = Ndark;
    th.nstar = Nstar;
    write_tipsy_standard_header(&xdrs,&th);
    /*
    ** Write out dark matter particles
    */
    if (bh->mass != 0) {
	for (k = 0; k < 3; k++) {
	    dpdpp.pos[k] = bh->r[k+1];
	    dpdpp.vel[k] = bh->v[k+1];
	    }
	dpdpp.mass = bh->mass;
	dpdpp.eps = bh->soft;
	dpdpp.phi = bh->Epot;
	write_tipsy_standard_dark_dpp(&xdrs,&dpdpp);
	}
    if (gi->do_halo == 1) {
	for (j = 0; j < (halo->Nshell+2); j++) {
	    p = halo->shell[j].p;
	    for (i = 0; i < halo->shell[j].N; i++) {
		for (k = 0; k < 3; k++) {
		    dpdpp.pos[k] = p[i].r[k+1];
		    dpdpp.vel[k] = p[i].v[k+1];
		    }
		dpdpp.mass = p[i].mass;
		dpdpp.eps = p[i].soft;
		dpdpp.phi = p[i].Epot;
		write_tipsy_standard_dark_dpp(&xdrs,&dpdpp);
		}
	    }
	}
    if (gi->do_bulge == 1) {
	for (j = 0; j < (bulge->Nshell+2); j++) {
	    p = bulge->shell[j].p;
	    for (i = 0; i < bulge->shell[j].N; i++) {
		for (k = 0; k < 3; k++) {
		    spdpp.pos[k] = p[i].r[k+1];
		    spdpp.vel[k] = p[i].v[k+1];
		    }
		spdpp.mass = p[i].mass;
		spdpp.eps = p[i].soft;
		spdpp.phi = p[i].Epot;
		spdpp.metals = 0;
		spdpp.tform = 0;
		write_tipsy_standard_star_dpp(&xdrs,&spdpp);
		}
	    }
	}
    xdr_destroy(&xdrs);
    }

/*
** Routine for writing some general output
*/

void write_general_output(FILE *file, int argc, char **argv, GI *gi, const PARTICLE *bh, const SI *bulge, const SI *halo) {

    INT i;
    DOUBLE temp;

    fprintf(file,"HALOGEN "VERSION" by Marcel Zemp\n\n");
    fprintf(file,"This version works in units where [G] = 1 LU^3 TU^-2 MU^-1.\n");
    fprintf(file,"LU = kpc / TU = Gyr / VU = kpc Gyr^-1 / MU = %e Mo\n\n",MU);
    fprintf(file,"Command line\n\n");
    for (i = 0; i < argc; i++) fprintf(file,"%s ",argv[i]);
    fprintf(file,"\n\n");
    fprintf(file,"Cosmological parameters\n\n");
    fprintf(file,"OmegaM0   = "OFD1"\n",gi->OmegaM0);
    fprintf(file,"OmegaK0   = "OFD1"\n",gi->OmegaK0);
    fprintf(file,"OmegaL0   = "OFD1"\n",gi->OmegaL0);
    fprintf(file,"h0        = "OFD1"\n",gi->h0);
    fprintf(file,"z         = "OFD1"\n",gi->z);
    fprintf(file,"Deltavirz = "OFD3"\n",gi->Deltavirz);
    fprintf(file,"rhocritz  = "OFD3" MU LU^-3 = "OFD3" Mo LU^-3\n",gi->rhocritz,gi->rhocritz*MU);
    fprintf(file,"rhovirz   = "OFD3" MU LU^-3 = "OFD3" Mo LU^-3\n\n",gi->rhocritz*gi->Deltavirz,gi->rhocritz*gi->Deltavirz*MU);
    if (bh->mass > 0) {
	fprintf(file,"Black hole properties\n\n");
	fprintf(file,"MBH    = "OFD3" MU = "OFD3" Mo\n",bh->mass,bh->mass*MU);
	fprintf(file,"softBH = "OFD3" LU\n\n",bh->soft);
	}
    if (gi->do_bulge == 1) {
	write_output_system(file,gi,bulge);
	}
    if (gi->do_halo == 1) {
	write_output_system(file,gi,halo);
	}
    fprintf(file,"General properties\n\n");
    if (gi->gridr->eqrvcmax[0] < 0) {
	temp = 2*pow(1000*gi->vcmax/(100*gi->h0*Ecosmo(gi)*gi->rvcmax*VelConvertFac),2);
	fprintf(file,"cV     = "OFD3"\n",temp);
	fprintf(file,"vcmax  = "OFD3" LU TU^-1\n",gi->vcmax);
	fprintf(file,"rvcmax = "OFD3" LU\n\n",gi->rvcmax);
	temp = Menc_total(gi->rvcmax,gi);
	fprintf(file,"M(rvcmax)        = "OFD3" MU = "OFD3" Mo\n",temp,temp*MU);
	}
    if (gi->do_bulge == 1) {
	temp = Menc_total(bulge->r1,gi);
	fprintf(file,"M(r1_bulge)      = "OFD3" MU = "OFD3" Mo\n",temp,temp*MU);
	temp = Menc_total(bulge->sp->rs,gi);
	fprintf(file,"M(rs_bulge)      = "OFD3" MU = "OFD3" Mo\n",temp,temp*MU);
	if (bulge->eqrvcmax[0] < 0) {
	    temp = Menc_total(bulge->sp->rvcmax,gi);
	    fprintf(file,"M(rvcmax_bulge)  = "OFD3" MU = "OFD3" Mo\n",temp,temp*MU);
	    }
	temp = Menc_total(bulge->sp->rhalf,gi);
	fprintf(file,"M(rhalf_bulge)   = "OFD3" MU = "OFD3" Mo\n",temp,temp*MU);
	if (bulge->sp->rcutoff != SBI) {
	    temp = Menc_total(bulge->sp->rcutoff,gi);
	    fprintf(file,"M(rcutoff_bulge) = "OFD3" MU = "OFD3" Mo\n",temp,temp*MU);
	    }
	temp = Menc_total(bulge->sp->rvir,gi);
	fprintf(file,"M(rvir_bulge)    = "OFD3" MU = "OFD3" Mo\n",temp,temp*MU);
	}
    if (gi->do_halo == 1) {
	temp = Menc_total(halo->r1,gi);
	fprintf(file,"M(r1_halo)       = "OFD3" MU = "OFD3" Mo\n",temp,temp*MU);
	temp = Menc_total(halo->sp->rs,gi);
	fprintf(file,"M(rs_halo)       = "OFD3" MU = "OFD3" Mo\n",temp,temp*MU);
	if (halo->eqrvcmax[0] < 0) {
	    temp = Menc_total(halo->sp->rvcmax,gi);
	    fprintf(file,"M(rvcmax_halo)   = "OFD3" MU = "OFD3" Mo\n",temp,temp*MU);
	    }
	temp = Menc_total(halo->sp->rhalf,gi);
	fprintf(file,"M(rhalf_halo)    = "OFD3" MU = "OFD3" Mo\n",temp,temp*MU);
	if (halo->sp->rcutoff != SBI) {
	    temp = Menc_total(halo->sp->rcutoff,gi);
	    fprintf(file,"M(rcutoff_halo)  = "OFD3" MU = "OFD3" Mo\n",temp,temp*MU);
	    }
	temp = Menc_total(halo->sp->rvir,gi);
	fprintf(file,"M(rvir_halo)     = "OFD3" MU = "OFD3" Mo\n",temp,temp*MU);
	}
    temp = Menc_total(gi->rinner,gi);
    fprintf(file,"M(rinner)        = "OFD3" MU = "OFD3" Mo\n",temp,temp*MU);
    temp = Menc_total(gi->router,gi);
    fprintf(file,"M(router)        = "OFD3" MU = "OFD3" Mo\n",temp,temp*MU);
    fprintf(file,"\n");
    if (gi->gridr->eqrvcmax[0] < 0) {
	temp = Tdyn_total(gi->rvcmax,gi);
	fprintf(file,"Tdyn(rvcmax)        = "OFD3" TU\n",temp);
	}
    if (gi->do_bulge == 1) {
	temp = Tdyn_total(bulge->r1,gi);
	fprintf(file,"Tdyn(r1_bulge)      = "OFD3" TU\n",temp);
	temp = Tdyn_total(bulge->sp->rs,gi);
	fprintf(file,"Tdyn(rs_bulge)      = "OFD3" TU\n",temp);
	if (bulge->eqrvcmax[0] < 0) {
	    temp = Tdyn_total(bulge->sp->rvcmax,gi);
	    fprintf(file,"Tdyn(rvcmax_bulge)  = "OFD3" TU\n",temp);
	    }
	temp = Tdyn_total(bulge->sp->rhalf,gi);
	fprintf(file,"Tdyn(rhalf_bulge)   = "OFD3" TU\n",temp);
	if (bulge->sp->rcutoff != SBI) {
	    temp = Tdyn_total(bulge->sp->rcutoff,gi);
	    fprintf(file,"Tdyn(rcutoff_bulge) = "OFD3" TU\n",temp);
	    }
	temp = Tdyn_total(bulge->sp->rvir,gi);
	fprintf(file,"Tdyn(rvir_bulge)    = "OFD3" TU\n",temp);
	}
    if (gi->do_halo == 1) {
	temp = Tdyn_total(halo->r1,gi);
	fprintf(file,"Tdyn(r1_halo)       = "OFD3" TU\n",temp);
	temp = Tdyn_total(halo->sp->rs,gi);
	fprintf(file,"Tdyn(rs_halo)       = "OFD3" TU\n",temp);
	if (halo->eqrvcmax[0] < 0) {
	    temp = Tdyn_total(halo->sp->rvcmax,gi);
	    fprintf(file,"Tdyn(rvcmax_halo)   = "OFD3" TU\n",temp);
	    }
	temp = Tdyn_total(halo->sp->rhalf,gi);
	fprintf(file,"Tdyn(rhalf_halo)    = "OFD3" TU\n",temp);
	if (halo->sp->rcutoff != SBI) {
	    temp = Tdyn_total(halo->sp->rcutoff,gi);
	    fprintf(file,"Tdyn(rcutoff_halo)  = "OFD3" TU\n",temp);
	    }
	temp = Tdyn_total(halo->sp->rvir,gi);
	fprintf(file,"Tdyn(rvir_halo)     = "OFD3" TU\n",temp);
	}
    fprintf(file,"\n");
    fprintf(file,"|Cr|      = "OFD3" LU            Cr      = ("OFD4", "OFD4", "OFD4") LU\n",gi->samp->Cr[0],gi->samp->Cr[1],gi->samp->Cr[2],gi->samp->Cr[3]);
    fprintf(file,"|Cv|      = "OFD3" LU TU^-1      Cv      = ("OFD4", "OFD4", "OFD4") LU TU^-1\n",gi->samp->Cv[0],gi->samp->Cv[1],gi->samp->Cv[2],gi->samp->Cv[3]);
 fprintf(file,"|L|       = "OFD3" MU LU^2 TU^-1 L       = ("OFD4", "OFD4", "OFD4") MU LU^2 TU^-1\n",gi->samp->Ltot[0],gi->samp->Ltot[1],gi->samp->Ltot[2],gi->samp->Ltot[3]);
    fprintf(file,"|L|       = "OFD3" Mo LU^2 TU^-1 L       = ("OFD4", "OFD4", "OFD4") Mo LU^2 TU^-1\n",gi->samp->Ltot[0]*MU,gi->samp->Ltot[1]*MU,gi->samp->Ltot[2]*MU,gi->samp->Ltot[3]*MU);
    fprintf(file,"|L|/Msamp = "OFD3" LU^2 TU^-1    L/Msamp = ("OFD4", "OFD4", "OFD4") LU^2 TU^-1\n",
	    gi->samp->Ltot[0]/gi->samp->Mp,gi->samp->Ltot[1]/gi->samp->Mp,gi->samp->Ltot[2]/gi->samp->Mp,gi->samp->Ltot[3]/gi->samp->Mp);
    fprintf(file,"\n");
    fprintf(file,"Etot = "OFD4" MU LU^2 TU^-2 = "OFD4" Mo LU^2 TU^-2\n",gi->samp->Etot,gi->samp->Etot*MU);
    fprintf(file,"Ekin = "OFD4" MU LU^2 TU^-2 = "OFD4" Mo LU^2 TU^-2\n",gi->samp->Ekin,gi->samp->Ekin*MU);
    fprintf(file,"Epot = "OFD4" MU LU^2 TU^-2 = "OFD4" Mo LU^2 TU^-2\n",gi->samp->Epot,gi->samp->Epot*MU);
    fprintf(file,"Rvir = |2*Ekin/Epot| = %g\n",fabs(2*gi->samp->Ekin/gi->samp->Epot));
    fprintf(file,"\n");
    fprintf(file,"Ntot        = "OFD3" = "OFI1"\n",(DOUBLE)gi->samp->Ntot,gi->samp->Ntot);
    fprintf(file,"Neff        = "OFD3"\n",(DOUBLE)gi->samp->Neff);
    fprintf(file,"Ninitialtot = "OFD3" = "OFI1"\n",(DOUBLE)gi->samp->Ninitialtot,gi->samp->Ninitialtot);
    fprintf(file,"Nnosplittot = "OFD3" = "OFI1"\n",(DOUBLE)gi->samp->Nnosplittot,gi->samp->Nnosplittot);
    fprintf(file,"Nnewtot     = "OFD3" = "OFI1"\n",(DOUBLE)gi->samp->Nnewtot,gi->samp->Nnewtot);
    fprintf(file,"Nfesm       = "OFD3" TU^-1\n",gi->samp->Nfesm);
    fprintf(file,"Nfemm       = "OFD3" TU^-1\n",gi->samp->Nfemm);
    fprintf(file,"sth         = "OFD3"\n",gi->samp->Nfesm/gi->samp->Nfemm);
    fprintf(file,"\n");
    fprintf(file,"Random seed = "OFD3"\n",gi->randomseed);
    fprintf(file,"Ngridr      = "OFI1"\n",gi->Ngridr);
    fprintf(file,"Ngriddf     = "OFI1"\n",gi->Ngriddf);
    fprintf(file,"rinner      = "OFD3" LU\n",gi->rinner);
    fprintf(file,"router      = "OFD3" LU\n\n",gi->router);
    fprintf(file,"Times for individual steps\n\n");
    fprintf(file,"Calculation of halo properties and initialisation of grid in r: "OFD1" seconds.\n",gi->t[1]-gi->t[0]);
    fprintf(file,"Initialisation of grid for distribution function: "OFD1" seconds.\n",gi->t[2]-gi->t[1]);
    fprintf(file,"Initialisation of shells: "OFD1" seconds.\n",gi->t[3]-gi->t[2]);
    fprintf(file,"Setting particle positions: "OFD1" seconds\n",gi->t[4]-gi->t[3]);
    fprintf(file,"Setting particle velocities: "OFD1" seconds\n",gi->t[5]-gi->t[4]);
    fprintf(file,"Setting remaining particle attributes: "OFD1" seconds\n",gi->t[6]-gi->t[5]);
    fprintf(file,"Orbit refinement: "OFD1" seconds\n",gi->t[7]-gi->t[6]);
    fprintf(file,"Calculating a few things and correct center of mass: "OFD1" seconds\n",gi->t[8]-gi->t[7]);
    gi->t[9] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(file,"Writing output: "OFD1" seconds\n",gi->t[9]-gi->t[8]);
    fprintf(file,"Total time: "OFD1" seconds\n",gi->t[9]-gi->t[0]);
    }

void write_output_system(FILE *file, const GI *gi, const SI *si) {

    INT i, k;
    DOUBLE temp;

    fprintf(file,"Model properties for %s\n\n",si->systemname);
    fprintf(file,"alpha = "OFD1"\n",si->sp->alpha);
    fprintf(file,"beta  = "OFD1"\n",si->sp->beta);
    fprintf(file,"gamma = "OFD1"\n",si->sp->gamma);
    fprintf(file,"rho0  = "OFD3" MU LU^-3 = "OFD3" Mo LU^-3\n",si->sp->rho0,si->sp->rho0*MU);
    fprintf(file,"cvir  = "OFD3"\n",si->sp->cvir);
    if (si->eqrvcmax[0] < 0) {
	temp = 2*pow(1000*si->sp->vcmax/(100*gi->h0*Ecosmo(gi)*si->sp->rvcmax*VelConvertFac),2);
	fprintf(file,"cV    = "OFD3"\n",temp);
	fprintf(file,"vcmax = "OFD3" LU TU^-1\n",si->sp->vcmax);
	}
    fprintf(file,"\n");
    fprintf(file,"rs      = "OFD3" LU\n",si->sp->rs);
    if (si->eqrvcmax[0] < 0) {
	fprintf(file,"rvcmax  = "OFD3" LU\n",si->sp->rvcmax);
	}
    fprintf(file,"rhalf   = "OFD3" LU\n",si->sp->rhalf);
    if (si->sp->rcutoff != SBI) {
	fprintf(file,"rcutoff = "OFD3" LU\n",si->sp->rcutoff);
	}
    fprintf(file,"rvir    = "OFD3" LU\n",si->sp->rvir);
    fprintf(file,"rsi     = "OFD3" LU\n",si->rsi);
    fprintf(file,"rso     = "OFD3" LU\n",si->rso);
    fprintf(file,"rmor    = "OFD3" LU\n\n",si->rmor);
    temp = Menc_system(si->sp->rs,gi,si);
    fprintf(file,"M_%s(rs)      = "OFD3" MU = "OFD3" Mo\n",si->systemname,temp,temp*MU);
    if (gi->gridr->eqrvcmax[0] < 0) {
	temp = Menc_system(si->sp->rvcmax,gi,si);
	fprintf(file,"M_%s(rvcmax)  = "OFD3" MU = "OFD3" Mo\n",si->systemname,temp,temp*MU);
	}
    temp = Menc_system(si->sp->rhalf,gi,si);
    fprintf(file,"M_%s(rhalf)   = "OFD3" MU = "OFD3" Mo\n",si->systemname,temp,temp*MU);
    if (si->sp->rcutoff != SBI) {
	temp = Menc_system(si->sp->rcutoff,gi,si);
	fprintf(file,"M_%s(rcutoff) = "OFD3" MU = "OFD3" Mo\n",si->systemname,temp,temp*MU);
	}
    temp = Menc_system(si->sp->rvir,gi,si);
    fprintf(file,"M_%s(rvir)    = "OFD3" MU = "OFD3" Mo\n",si->systemname,temp,temp*MU);
    temp = Menc_system(si->rsi,gi,si);
    fprintf(file,"M_%s(rsi)     = "OFD3" MU = "OFD3" Mo\n",si->systemname,temp,temp*MU);
    temp = Menc_system(si->rso,gi,si);
    fprintf(file,"M_%s(rso)     = "OFD3" MU = "OFD3" Mo\n",si->systemname,temp,temp*MU);
    if (si->rmor != 0) {
	temp = Menc_system(si->rmor,gi,si);
	fprintf(file,"M_%s(rmor)    = "OFD3" MU = "OFD3" Mo\n",si->systemname,temp,temp*MU);
	}
    temp = Menc_system(gi->rinner,gi,si);
    fprintf(file,"M_%s(rinner)  = "OFD3" MU = "OFD3" Mo\n",si->systemname,temp,temp*MU);
    temp = Menc_system(gi->router,gi,si);
    fprintf(file,"M_%s(router)  = "OFD3" MU = "OFD3" Mo\n",si->systemname,temp,temp*MU);
    fprintf(file,"\n");
    fprintf(file,"Sampling properties for %s\n\n",si->systemname);
    fprintf(file,"|Cr|      = "OFD3" LU            Cr      = ("OFD4", "OFD4", "OFD4") LU\n",si->samp->Cr[0],si->samp->Cr[1],si->samp->Cr[2],si->samp->Cr[3]);
    fprintf(file,"|Cv|      = "OFD3" LU TU^-1      Cv      = ("OFD4", "OFD4", "OFD4") LU TU^-1\n",si->samp->Cv[0],si->samp->Cv[1],si->samp->Cv[2],si->samp->Cv[3]);
    fprintf(file,"|L|       = "OFD3" MU LU^2 TU^-1 L       = ("OFD4", "OFD4", "OFD4") MU LU^2 TU^-1\n",si->samp->Ltot[0],si->samp->Ltot[1],si->samp->Ltot[2],si->samp->Ltot[3]);
    fprintf(file,"|L|       = "OFD3" Mo LU^2 TU^-1 L       = ("OFD4", "OFD4", "OFD4") Mo LU^2 TU^-1\n",si->samp->Ltot[0]*MU,si->samp->Ltot[1]*MU,si->samp->Ltot[2]*MU,si->samp->Ltot[3]*MU);
    fprintf(file,"|L|/Msamp = "OFD3" LU^2 TU^-1    L/Msamp = ("OFD4", "OFD4", "OFD4") LU^2 TU^-1\n",
	    si->samp->Ltot[0]/si->samp->Mp,si->samp->Ltot[1]/si->samp->Mp,si->samp->Ltot[2]/si->samp->Mp,si->samp->Ltot[3]/si->samp->Mp);
    fprintf(file,"\n");
    fprintf(file,"i   rshelli [LU]  rshello [LU]  N          Ninitial   Nnosplit   Nnew       Mtheo [MU]    Msamp [MU]    massmax [MU]  softmax [LU]\n");
    if ((si->sp->beta > 3) && (si->rsi == gi->router)) {
	k = si->Nshell+1;
	}
    else {
	k = si->Nshell+2;
	}
    for (i = 0; i < k; i++) {
	fprintf(file,OFI2" "OFD3" "OFD3" "OFI3" "OFI3" "OFI3" "OFI3" "OFD3" "OFD3" "OFD3" "OFD3"\n",
		i,si->shell[i].rinner,si->shell[i].router,si->shell[i].N,si->shell[i].Ninitial,si->shell[i].Nnosplit,si->shell[i].Nnew,(si->shell[i+1].Menc-si->shell[i].Menc),si->shell[i].Mp,si->shell[i].mass,si->shell[i].soft);
	}
    fprintf(file,"\n");
    fprintf(file,"Nshell              = "OFI1"\n",si->Nshell);
    fprintf(file,"Ismor               = "OFI1"\n",si->Ismor);
    fprintf(file,"DRMmax              = "OFI1"\n",si->DRMmax);
    fprintf(file,"Ntot                = "OFD3" = "OFI1"\n",(DOUBLE)si->samp->Ntot,si->samp->Ntot);
    fprintf(file,"Neff                = "OFD3"\n",si->samp->Neff);
    if (si->sp->beta <= 3) {
	fprintf(file,"Neff within rcutoff = "OFD3"\n",Menc_system(si->sp->rcutoff,gi,si)/si->shell[0].mass);
	}
    fprintf(file,"Ninitialtot         = "OFD3" = "OFI1"\n",(DOUBLE)si->samp->Ninitialtot,si->samp->Ninitialtot);
    fprintf(file,"Nnosplittot         = "OFD3" = "OFI1"\n",(DOUBLE)si->samp->Nnosplittot,si->samp->Nnosplittot);
    fprintf(file,"Nnewtot             = "OFD3" = "OFI1"\n",(DOUBLE)si->samp->Nnewtot,si->samp->Nnewtot);
    fprintf(file,"Nfesm               = "OFD3" TU^-1\n",si->samp->Nfesm);
    fprintf(file,"Nfemm               = "OFD3" TU^-1\n",si->samp->Nfemm);
    fprintf(file,"sth                 = "OFD3"\n",si->samp->Nfesm/si->samp->Nfemm);
    fprintf(file,"rimp                = "OFD3" LU\n",si->rimp);
    fprintf(file,"r1                  = "OFD3" LU\n",si->r1);
    fprintf(file,"r100                = "OFD3" LU\n",si->r100);
    temp = Menc_system(gi->router,gi,si);
    fprintf(file,"Mtheo               = "OFD3" MU = "OFD3" Mo\n",temp,temp*MU);
    fprintf(file,"Msamp               = "OFD3" MU = "OFD3" Mo\n",si->samp->Mp,si->samp->Mp*MU);
    fprintf(file,"(Msamp-Mtheo)/Mtheo = "OFD3"\n",si->samp->Mp/temp-1.0);
    fprintf(file,"DF split factor     = "OFD3"\n\n",si->dfsf);
    }
