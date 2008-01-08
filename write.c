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

/*
** Routine for writing tipsy standard format memory efficient
*/

void write_tipsy_standard_2(FILE *fp, const PARTICLE *bh, const SI *halo) {

    INT i, j, k, Ntotal;
    TIPSY_HEADER th;
    DARK_PARTICLE dp;
    PARTICLE *p;
    XDR xdrs;

    xdrstdio_create(&xdrs,fp,XDR_ENCODE);
    /*
    ** Write out header
    */
    Ntotal = 0;
    if (bh->mass != 0) {
	Ntotal++;
	}
    for (j = 0; j < (halo->Nshell+2); j++) {
	Ntotal = Ntotal + halo->shell[j].N;
	}
    th.time = 0;
    th.ntotal = Ntotal;
    th.ndim = 3;
    th.ngas = 0;
    th.ndark = Ntotal;
    th.nstar = 0;
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
    xdr_destroy(&xdrs);
    }

/*
** Routines for writing out grids
*/

void write_gridr(FILE *file, const GI *gi) {

    INT i;
    GRIDR *gridr;

    gridr = gi->gridr;
    for (i = 0; i < gi->Ngridr; i++) {
	fprintf(file,OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2"\n",
		gridr->r[i],gridr->logr[i],gridr->rho[i],gridr->logrho[i],gridr->rhoenc[i],gridr->logrhoenc[i],
		gridr->Menc[i],gridr->logMenc[i],gridr->Pot[i],gridr->logPot[i]);
	}
    }

void write_griddf(FILE *file, const GI *gi, const SI *si) {

    INT i;
    GRIDDF *griddf;

    griddf = si->griddf;
    for (i = 0; i < gi->Ngriddf; i++) {
	fprintf(file,OFD2" "OFD2" "OFD2" "OFD2" "OFD2" "OFD2"\n",
		griddf->r[i],griddf->logr[i],griddf->E[i],griddf->logE[i],griddf->fE[i],griddf->logfE[i]);
	}
    }

void write_general_output(FILE *file, int argc, char **argv, GI *gi, const PARTICLE *bh, const SI *halo) {

    INT i, k;

    fprintf(file,"HALOGEN "VERSION" by Marcel Zemp\n");
    fprintf(file,"This version works in units where [G] = 1 [L]^3 [T]^-2 [M]^-1.\n");
    fprintf(file,"[L] = kpc / [T] = Gyr / [V] = kpc Gyr^-1 / [M] = MU = %e Mo\n\n",MU);
    fprintf(file,"Command line\n\n");
    for (i = 0; i < argc; i++) fprintf(file,"%s ",argv[i]);
    fprintf(file,"\n\n");
    fprintf(file,"Cosmological parameters\n\n");
    fprintf(file,"OmegaM0   = "OFD1"\n",gi->OmegaM0);
    fprintf(file,"OmegaK0   = "OFD1"\n",gi->OmegaK0);
    fprintf(file,"OmegaL0   = "OFD1"\n",gi->OmegaL0);
    fprintf(file,"h0        = "OFD1"\n",gi->h0);
    fprintf(file,"z         = "OFD1"\n",gi->z);
    fprintf(file,"rhocritz  = "OFD3" Mo kpc^-3 = "OFD3" MU kpc^-3\n",gi->rhocritz*MU,gi->rhocritz);
    fprintf(file,"Deltavirz = "OFD3"\n",gi->Deltavirz);
    fprintf(file,"rhovirz   = "OFD3" Mo kpc^-3 = "OFD3" MU kpc^-3\n\n",gi->rhocritz*gi->Deltavirz*MU,gi->rhocritz*gi->Deltavirz);
    fprintf(file,"Model properties\n\n");
    fprintf(file,"alpha = "OFD1"\n",halo->sp->alpha);
    fprintf(file,"beta  = "OFD1"\n",halo->sp->beta);
    fprintf(file,"gamma = "OFD1"\n",halo->sp->gamma);
    fprintf(file,"rho0  = "OFD3" Mo kpc^-3 = "OFD3" MU kpc^-3\n",halo->sp->rho0*MU,halo->sp->rho0);
    fprintf(file,"cvir  = "OFD3"\n",halo->sp->cvir);
    if (halo->sp->gamma < 2) {
	fprintf(file,"cV    = "OFD3"\n",2*pow(1000*halo->sp->vcmax/(100*gi->h0*Ecosmo(gi)*halo->sp->rvcmax*VelConvertFac),2));
	fprintf(file,"vcmax = "OFD3" kpc Gyr^-1 = "OFD3" km s^-1\n",halo->sp->vcmax,halo->sp->vcmax/VelConvertFac);
	}
    fprintf(file,"\n");
    if (bh->mass > 0) {
	fprintf(file,"MBH    = "OFD3" Mo = "OFD3" MU\n",bh->mass*MU,bh->mass);
	fprintf(file,"softBH = "OFD3" kpc\n",bh->soft);
	fprintf(file,"\n");
	}
    fprintf(file,"rs      = "OFD3" kpc\n",halo->sp->rs);
    if (gi->gridr->eqrvcmax[0] < 0) {
	fprintf(file,"rvcmax  = "OFD3" kpc\n",halo->sp->rvcmax);
	}
    fprintf(file,"rhalf   = "OFD3" kpc\n",halo->sp->rhalf);
    if (halo->sp->rcutoff != SBI) {
	fprintf(file,"rcutoff = "OFD3" kpc\n",halo->sp->rcutoff);
	}
    fprintf(file,"rvir    = "OFD3" kpc\n",halo->sp->rvir);
    fprintf(file,"rinner  = "OFD3" kpc\n",gi->rinner);
    fprintf(file,"router  = "OFD3" kpc\n",gi->router);
    fprintf(file,"rsi     = "OFD3" kpc\n",halo->rsi);
    fprintf(file,"rso     = "OFD3" kpc\n",halo->rso);
    fprintf(file,"rmor    = "OFD3" kpc\n",halo->rmor);
    fprintf(file,"rdecay  = "OFD3" kpc\n",halo->sp->rdecay);
    fprintf(file,"\n");
    fprintf(file,"Tdyn(rs)      = "OFD3" Gyr\n",Tdyn(halo->sp->rs,gi));
    if (gi->gridr->eqrvcmax[0] < 0) {
	fprintf(file,"Tdyn(rvcmax)  = "OFD3" Gyr\n",Tdyn(halo->sp->rvcmax,gi));
	}
    fprintf(file,"Tdyn(rhalf)   = "OFD3" Gyr\n",Tdyn(halo->sp->rhalf,gi));
    if (halo->sp->rcutoff != SBI) {
	fprintf(file,"Tdyn(rcutoff) = "OFD3" Gyr\n",Tdyn(halo->sp->rcutoff,gi));
	}
    fprintf(file,"Tdyn(rvir)    = "OFD3" Gyr\n",Tdyn(halo->sp->rvir,gi));
    fprintf(file,"Tdyn(rinner)  = "OFD3" Gyr\n",Tdyn(gi->rinner,gi));
    fprintf(file,"Tdyn(router)  = "OFD3" Gyr\n",Tdyn(gi->router,gi));
    fprintf(file,"Tdyn(rsi)     = "OFD3" Gyr\n",Tdyn(halo->rsi,gi));
    fprintf(file,"Tdyn(rso)     = "OFD3" Gyr\n",Tdyn(halo->rso,gi));
    if (halo->rmor != 0) {
	fprintf(file,"Tdyn(rmor)    = "OFD3" Gyr\n",Tdyn(halo->rmor,gi));
	}
    fprintf(file,"\n");
    fprintf(file,"M(rs)      = "OFD3" Mo = "OFD3" MU\n",Menc(halo->sp->rs,gi)*MU,Menc(halo->sp->rs,gi));
    if (gi->gridr->eqrvcmax[0] < 0) {
	fprintf(file,"M(rvcmax)  = "OFD3" Mo = "OFD3" MU\n",Menc(halo->sp->rvcmax,gi)*MU,Menc(halo->sp->rvcmax,gi));
	}
    fprintf(file,"M(rhalf)   = "OFD3" Mo = "OFD3" MU\n",Menc(halo->sp->rhalf,gi)*MU,Menc(halo->sp->rhalf,gi));
    if (halo->sp->rcutoff != SBI) {
	fprintf(file,"M(rcutoff) = "OFD3" Mo = "OFD3" MU\n",Menc(halo->sp->rcutoff,gi)*MU,Menc(halo->sp->rcutoff,gi));
	}
    fprintf(file,"M(rvir)    = "OFD3" Mo = "OFD3" MU\n",Menc(halo->sp->rvir,gi)*MU,Menc(halo->sp->rvir,gi));
    fprintf(file,"M(rinner)  = "OFD3" Mo = "OFD3" MU\n",Menc(gi->rinner,gi)*MU,Menc(gi->rinner,gi));
    fprintf(file,"M(router)  = "OFD3" Mo = "OFD3" MU\n",Menc(gi->router,gi)*MU,Menc(gi->router,gi));
    fprintf(file,"M(rsi)     = "OFD3" Mo = "OFD3" MU\n",Menc(halo->rsi,gi)*MU,Menc(halo->rsi,gi));
    fprintf(file,"M(rso)     = "OFD3" Mo = "OFD3" MU\n",Menc(halo->rso,gi)*MU,Menc(halo->rso,gi));
    if (halo->rmor != 0) {
	fprintf(file,"M(rmor)    = "OFD3" Mo = "OFD3" MU\n",Menc(halo->rmor,gi)*MU,Menc(halo->rmor,gi));
	}
    fprintf(file,"\n");
    fprintf(file,"Sampling properties\n\n");
    fprintf(file,"|Cr|      = "OFD3" kpc             Cr      = ("OFD4", "OFD4", "OFD4") kpc\n",gi->stuff->Cr[0],gi->stuff->Cr[1],gi->stuff->Cr[2],gi->stuff->Cr[3]);
    fprintf(file,"|Cv|      = "OFD3" kpc Gyr^-1      Cv      = ("OFD4", "OFD4", "OFD4") kpc Gyr^-1\n",gi->stuff->Cv[0],gi->stuff->Cv[1],gi->stuff->Cv[2],gi->stuff->Cv[3]);
    fprintf(file,"|L|       = "OFD3" Mo kpc^2 Gyr^-1 L       = ("OFD4", "OFD4", "OFD4") Mo kpc^2 Gyr^-1\n",gi->stuff->Ltot[0]*MU,gi->stuff->Ltot[1]*MU,gi->stuff->Ltot[2]*MU,gi->stuff->Ltot[3]*MU);
    fprintf(file,"|L|       = "OFD3" MU kpc^2 Gyr^-1 L       = ("OFD4", "OFD4", "OFD4") MU kpc^2 Gyr^-1\n",gi->stuff->Ltot[0],gi->stuff->Ltot[1],gi->stuff->Ltot[2],gi->stuff->Ltot[3]);
    fprintf(file,"|L|/Msamp = "OFD3" kpc^2 Gyr^-1    L/Msamp = ("OFD4", "OFD4", "OFD4") kpc^2 Gyr^-1\n",
	    gi->stuff->Ltot[0]/gi->stuff->Mp,gi->stuff->Ltot[1]/gi->stuff->Mp,gi->stuff->Ltot[2]/gi->stuff->Mp,gi->stuff->Ltot[3]/gi->stuff->Mp);
    fprintf(file,"\n");
    fprintf(file,"Etot = "OFD4" Mo kpc^2 Gyr^-2 = "OFD4" MU kpc^2 Gyr^-2\n",gi->stuff->Etot*MU,gi->stuff->Etot);
    fprintf(file,"Ekin = "OFD4" Mo kpc^2 Gyr^-2 = "OFD4" MU kpc^2 Gyr^-2\n",gi->stuff->Ekin*MU,gi->stuff->Ekin);
    fprintf(file,"Epot = "OFD4" Mo kpc^2 Gyr^-2 = "OFD4" MU kpc^2 Gyr^-2\n",gi->stuff->Epot*MU,gi->stuff->Epot);
    fprintf(file,"Rvir = |2*Ekin/Epot| = %g\n",fabs(2*gi->stuff->Ekin/gi->stuff->Epot));
    fprintf(file,"\n");
    fprintf(file,"i   rshelli [kpc] rshello [kpc] N          Ninitial   Nnosplit   Nnew       Mtheo [MU]    Msamp [MU]    massmax [MU]  softmax [kpc]\n");
    if ((halo->sp->beta > 3) && (halo->rsi == gi->router)) {
	k = halo->Nshell+1;
	}
    else {
	k = halo->Nshell+2;
	}
    for (i = 0; i < k; i++) {
	fprintf(file,OFI2" "OFD3" "OFD3" "OFI3" "OFI3" "OFI3" "OFI3" "OFD3" "OFD3" "OFD3" "OFD3"\n",
		i,halo->shell[i].rinner,halo->shell[i].router,halo->shell[i].N,halo->shell[i].Ninitial,halo->shell[i].Nnosplit,halo->shell[i].Nnew,(halo->shell[i+1].Menc-halo->shell[i].Menc),halo->shell[i].Mp,halo->shell[i].mass,halo->shell[i].soft);
	}
    fprintf(file,"\n");
    fprintf(file,"Nshell              = "OFI1"\n",halo->Nshell);
    fprintf(file,"Ismor               = "OFI1"\n",halo->Ismor);
    fprintf(file,"DRMmax              = "OFI1"\n",halo->DRMmax);
    fprintf(file,"Ntot                = "OFD3" = "OFI1"\n",(DOUBLE)gi->stuff->Ntot,gi->stuff->Ntot);
    fprintf(file,"Neff                = "OFD3"\n",Menc(gi->router,gi)/halo->shell[0].mass);
    if (halo->sp->beta <= 3) {
	fprintf(file,"Neff within rcutoff = "OFD3"\n",Menc(halo->sp->rcutoff,gi)/halo->shell[0].mass);
	}
    fprintf(file,"Ninitialtot         = "OFD3" = "OFI1"\n",(DOUBLE)gi->stuff->Ninitialtot,gi->stuff->Ninitialtot);
    fprintf(file,"Nnosplittot         = "OFD3" = "OFI1"\n",(DOUBLE)gi->stuff->Nnosplittot,gi->stuff->Nnosplittot);
    fprintf(file,"Nnewtot             = "OFD3" = "OFI1"\n",(DOUBLE)gi->stuff->Nnewtot,gi->stuff->Nnewtot);
    fprintf(file,"Nfesm               = "OFD3" Gyr^-1\n",gi->stuff->Nfesm);
    fprintf(file,"Nfemm               = "OFD3" Gyr^-1\n",gi->stuff->Nfemm);
    fprintf(file,"sth                 = "OFD3"\n",gi->stuff->Nfesm/gi->stuff->Nfemm);
    fprintf(file,"rimp                = "OFD3" kpc\n",halo->rimp);
    fprintf(file,"r1                  = "OFD3" kpc\n",halo->r1);
    fprintf(file,"r100                = "OFD3" kpc\n",halo->r100);
    fprintf(file,"Tdyn(rimp)          = "OFD3" Gyr\n",Tdyn(halo->rimp,gi));
    fprintf(file,"Tdyn(r1)            = "OFD3" Gyr\n",Tdyn(halo->r1,gi));
    fprintf(file,"Tdyn(r100)          = "OFD3" Gyr\n",Tdyn(halo->r100,gi));
    fprintf(file,"Mtheo               = "OFD3" Mo = "OFD3" MU\n",Menc(gi->router,gi)*MU,Menc(gi->router,gi));
    fprintf(file,"Msamp               = "OFD3" Mo = "OFD3" MU\n",gi->stuff->Mp*MU,gi->stuff->Mp);
    fprintf(file,"(Msamp-Mtheo)/Mtheo = "OFD3"\n",gi->stuff->Mp/Menc(gi->router,gi)-1.0);
    fprintf(file,"DF split factor     = "OFD3"\n",halo->dfsf);
    fprintf(file,"Random seed         = "OFD3"\n",gi->randomseed);
    fprintf(file,"\n");
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
