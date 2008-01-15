/*
** initalise.c
**
** Initialising routines for HALOGEN
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <time.h>
#include "definitions.h"
#include "functions.h"
#include "usage.h"

/*
** Routine for initialising general info structure
*/

void initialise_general(GI *gi) {

    INT i;

    gi->do_bh = 0;
    gi->do_bulge = 0;
    gi->do_halo = 0;
    gi->output_gridr = 0;
    gi->output_griddf = 0;
    gi->output_tipsy_standard = 0;
    gi->output_tipsy_standard_dpp = 0;
    gi->output_gadget_binary = 0;
    gi->positionsonly = 0;
    gi->Ngridr = 2001;
    gi->Ngriddf = 101;
    gi->OmegaM0 = 0.3;
    gi->OmegaK0 = 0;
    gi->OmegaL0 = 0.7;
    gi->h0 = 0.7;
    gi->z = 0;
    gi->rhocritz = -1;
    gi->Deltavirz = -1;
    gi->OmegaMz = -1;
    gi->rinner = -1;
    gi->router = -1;
    gi->f_rinner = 1e-6;
    gi->f_router = 1e20;
    gi->f_cutoff = 0.3;
    gi->rvcmax = -1;
    gi->vcmax = -1;
    gi->randomseed = time(NULL);
    /*
    ** Sampling info stuff
    */
    gi->samp = malloc(sizeof(SAMP));
    assert(gi->samp != NULL);
    gi->samp->Ntot = 0;
    gi->samp->Ninitialtot = 0;
    gi->samp->Nnosplittot = 0;
    gi->samp->Nnewtot = 0;
    gi->samp->Neff = 0;
    gi->samp->Mp = 0;
    gi->samp->Ekin = 0;
    gi->samp->Epot = 0;
    gi->samp->Nfemm = 0;
    gi->samp->Nfesm = 0;
    for(i = 0; i < 4; i++) {
	gi->samp->Cr[i] = 0;
	gi->samp->Cv[i] = 0;
	gi->samp->Ltot[i] = 0;
	}
    }

/*
** Routine for initialising a system
*/

void initialise_system(SI *si) {

    INT i;

    /*
    ** General stuff
    */
    si->rsi_in_rs_units = 0;
    si->rsi_in_rvir_units = 0;
    si->rsi_in_rcutoff_units = 0;
    si->rso_in_rs_units = 0;
    si->rso_in_rvir_units = 0;
    si->rso_in_rcutoff_units = 0;
    si->rmor_in_rs_units = 0;
    si->rmor_in_rvir_units = 0;
    si->rmor_in_rcutoff_units = 0;
    si->Nshell = 0;
    si->Ismor = -1;
    si->N0 = -1;
    si->DRMmax = 1;
    si->soft0 = -1;
    si->rsi = -1;
    si->rso = -1;
    si->rmor = -1;
    si->rimp = SBI;
    si->r1 = -1;
    si->r100 = -1;
    si->dfsf = 0.01;
    si->logrhoenc = NULL;
    si->logMenc = NULL;
    /*
    ** Particle stuff
    */
    si->sp = malloc(sizeof(SP));
    assert(si->sp != NULL);
    si->sp->alpha = -1;
    si->sp->beta = -1;
    si->sp->gamma = -1;
    si->sp->delta = -1;
    si->sp->rho0 = -1;
    si->sp->rs = -1;
    si->sp->rcutoff = -1;
    si->sp->rdecay = -1;
    si->sp->M = -1;
    si->sp->cvir = -1;
    si->sp->rvir = -1;
    si->sp->vvir = -1;
    si->sp->rhalf = -1;
    si->sp->rvcmax = -1;
    si->sp->vcmax = -1;
    /*
    ** Sampling info stuff
    */
    si->samp = malloc(sizeof(SAMP));
    assert(si->samp != NULL);
    si->samp->Ntot = 0;
    si->samp->Ninitialtot = 0;
    si->samp->Nnosplittot = 0;
    si->samp->Nnewtot = 0;
    si->samp->Neff = 0;
    si->samp->Mp = 0;
    si->samp->Ekin = 0;
    si->samp->Epot = 0;
    si->samp->Nfemm = 0;
    si->samp->Nfesm = 0;
    for(i = 0; i < 4; i++) {
	si->samp->Cr[i] = 0;
	si->samp->Cv[i] = 0;
	si->samp->Ltot[i] = 0;
	}
    }

/*
** Routine for initialising a particle
*/

void initialise_particle(PARTICLE *p) {

    INT i;

    for (i = 0; i < 4; i++) {
	p->r[i] = 0;
	p->v[i] = 0;
	p->L[i] = 0;
	}
    p->mass = 0;
    p->soft = 0;
    p->Ekin = 0;
    p->Epot = 0;
    p->Etot = 0;
    p->rperi = 0;
    p->rapo = 0;
    p->ecc = 0;
    }

/* 
** Routine for initialising gridr 
*/

void initialise_gridr(GI *gi, PARTICLE *bh, SI *bulge, SI *halo) {

    INT i;
    DOUBLE dlogr, logr, r, r3;
    DOUBLE Mencr, MencHalor, DeltaMencHalor, MencBulger, DeltaMencBulger;
    DOUBLE rhor, rhoHalor, rhoBulger;
    DOUBLE rhoencr, rhoencHalor, rhoencBulger;
    DOUBLE Potr, Potoutr;
    GRIDR *gridr;
    SP *bp, *hp;

    gridr = gi->gridr;
    bp = bulge->sp;
    hp = halo->sp;
    /*
    ** Determine gridr size
    */
    if (gi->do_bulge == 1) {
	if (gi->do_halo == 1) {
	    if (bp->rs < hp->rs) {
		gi->rinner = gi->f_rinner*bp->rs;
		}
	    else {
		gi->rinner = gi->f_rinner*hp->rs;
		}
	    }
	else {
	    gi->rinner = gi->f_rinner*bp->rs;
	    }
	}
    else if (gi->do_halo == 1) {
	gi->rinner = gi->f_rinner*hp->rs;
	}
    if (gi->rinner == -1) {
	fprintf(stderr,"HALOGEN is confused!\n");
	fprintf(stderr,"Missing or bad parameter.\n");
	fprintf(stderr,"Could not determine rinner for grid in r.\n");
	usage();
	}
    r = sqrt(10); /* use r as factor here temporarily */
    if (gi->do_bulge == 1) {
	if (gi->do_halo == 1) {
	    /*
	    ** Now we have both: bulge & halo and we have to be 
	    ** careful with the choice of router. We in general assume
	    ** that the bulge is less extended than the halo.
	    */
	    if ((bp->beta <= 3) && (hp->beta > 3)) {
		gi->router = bp->rs;
		while (rho(bp->rs,bulge)/rho(gi->router,bulge) < gi->f_router) {
		    gi->router = gi->router*r;
		    }
		}
	    else {
		gi->router = hp->rs;
		while (rho(hp->rs,halo)/rho(gi->router,halo) < gi->f_router) {
		    gi->router = gi->router*r;
		    }
		}
	    }
	else {
	    gi->router = bp->rs;
	    while (rho(bp->rs,bulge)/rho(gi->router,bulge) < gi->f_router) {
		gi->router = gi->router*r;
		}
	    }
	}
    else if (gi->do_halo == 1) {
	gi->router = hp->rs;
	while (rho(hp->rs,halo)/rho(gi->router,halo) < gi->f_router) {
	    gi->router = gi->router*r;
	    }
	}
    if (gi->router == -1) {
	fprintf(stderr,"HALOGEN is confused!\n");
	fprintf(stderr,"Missing or bad parameter.\n");
	fprintf(stderr,"Could not determine router for grid in r.\n");
	usage();
	}
    /*
    ** Initialise gridr
    */
    dlogr = (log(gi->router)-log(gi->rinner))/(gi->Ngridr-1);
    i = 0;
    logr = log(gi->rinner);
    r = exp(logr);
    r3 = r*r*r;
    if (gi->do_bulge == 1) {
	rhoBulger = rho(r,bulge);
	DeltaMencBulger = 4*M_PI*rhoBulger*r3/(3-bp->gamma); /* approximate inner gridpoint by analytical calculation */
	}
    else {
	rhoBulger = 0;
	DeltaMencBulger = 0;
	}
    if (gi->do_halo == 1) {
	rhoHalor = rho(r,halo);
	DeltaMencHalor = 4*M_PI*rhoHalor*r3/(3-hp->gamma); /* approximate inner gridpoint by analytical calculation */
	}
    else {
	rhoHalor = 0;
	DeltaMencHalor = 0;
	}
    rhor = rhoHalor + rhoBulger;
    MencBulger = DeltaMencBulger;
    MencHalor = DeltaMencHalor;
    Mencr = bh->mass + DeltaMencBulger + DeltaMencHalor;
    rhoencBulger = MencBulger/(4*M_PI*r3/3.0);
    rhoencHalor = MencHalor/(4*M_PI*r3/3.0);
    rhoencr = Mencr/(4*M_PI*r3/3.0);
    gridr->r[i] = r;
    gridr->logr[i] = logr;
    gridr->rho[i] = rhor;
    gridr->logrho[i] = log(rhor);
    gridr->rhoBulge[i] = rhoBulger;
    gridr->logrhoBulge[i] = log(rhoBulger);
    gridr->rhoHalo[i] = rhoHalor;
    gridr->logrhoHalo[i] = log(rhoHalor);
    gridr->rhoenc[i] = rhoencr;
    gridr->logrhoenc[i] = log(rhoencr);
    gridr->rhoencBulge[i] = rhoencBulger;
    gridr->logrhoencBulge[i] = log(rhoencBulger);
    gridr->rhoencHalo[i] = rhoencHalor;
    gridr->logrhoencHalo[i] = log(rhoencHalor);
    gridr->Menc[i] = Mencr;
    gridr->logMenc[i] = log(Mencr);
    gridr->MencBulge[i] = MencBulger;
    gridr->logMencBulge[i] = log(MencBulger);
    gridr->MencHalo[i] = MencHalor;
    gridr->logMencHalo[i] = log(MencHalor);
    gridr->eqrvcmax[i] = Mencr - 4*M_PI*rhor*r3;
    gridr->eqrvcmaxBulge[i] = MencBulger - 4*M_PI*rhoBulger*r3;
    gridr->eqrvcmaxHalo[i] = MencHalor - 4*M_PI*rhoHalor*r3;
    for (i = 1; i < gi->Ngridr; i++) {
	logr = log(gi->rinner) + i*dlogr;
	r = exp(logr);
	r3 = r*r*r;
	if (gi->do_bulge == 1) {
	    rhoBulger = rho(r,bulge);
	    DeltaMencBulger = integral(integrandMenc,gridr->r[i-1],r,bulge);
	    }
	else {
	    rhoBulger = 0;
	    DeltaMencBulger = 0;
	    }
	if (gi->do_halo == 1) {
	    rhoHalor = rho(r,halo);
	    DeltaMencHalor = integral(integrandMenc,gridr->r[i-1],r,halo);
	    }
	else {
	    rhoHalor = 0;
	    DeltaMencHalor = 0;
	    }
	rhor = rhoHalor + rhoBulger;
	MencBulger = gridr->MencBulge[i-1] + DeltaMencBulger;
	MencHalor = gridr->MencHalo[i-1] + DeltaMencHalor;
	Mencr = gridr->Menc[i-1] + DeltaMencBulger + DeltaMencHalor;
	rhoencBulger = MencBulger/(4*M_PI*r3/3.0);
	rhoencHalor = MencHalor/(4*M_PI*r3/3.0);
	rhoencr = Mencr/(4*M_PI*r3/3.0);
	gridr->r[i] = r;
	gridr->logr[i] = logr;
	gridr->rho[i] = rhor;
	gridr->logrho[i] = log(rhor);
	gridr->rhoBulge[i] = rhoBulger;
	gridr->logrhoBulge[i] = log(rhoBulger);
	gridr->rhoHalo[i] = rhoHalor;
	gridr->logrhoHalo[i] = log(rhoHalor);
	gridr->rhoenc[i] = rhoencr;
	gridr->logrhoenc[i] = log(rhoencr);
	gridr->rhoencBulge[i] = rhoencBulger;
	gridr->logrhoencBulge[i] = log(rhoencBulger);
	gridr->rhoencHalo[i] = rhoencHalor;
	gridr->logrhoencHalo[i] = log(rhoencHalor);
	gridr->Menc[i] = Mencr;
	gridr->logMenc[i] = log(Mencr);
	gridr->MencBulge[i] = MencBulger;
	gridr->logMencBulge[i] = log(MencBulger);
	gridr->MencHalo[i] = MencHalor;
	gridr->logMencHalo[i] = log(MencHalor);
	gridr->eqrvcmax[i] = Mencr - 4*M_PI*rhor*r3;
	gridr->eqrvcmaxBulge[i] = MencBulger - 4*M_PI*rhoBulger*r3;
	gridr->eqrvcmaxHalo[i] = MencHalor - 4*M_PI*rhoHalor*r3;
	}
    i = gi->Ngridr-1;
    /* 
    ** Analytic approximation for models without exponential cutoff.
    ** Analytic approximation for models with exponential cutoff
    ** would neet WhittakerM function. Hence, set Potoutr = 0.
    */
    Potoutr = 0;
    if (gi->do_bulge == 1) {
	if (bp->beta > 3) {
	    Potoutr += 4*M_PI*G*rho(gi->router,bulge)*(gi->router*gi->router)/(2-bp->beta); /* minus cancel */
	    }
	}
    if (gi->do_halo == 1) {
	if (hp->beta > 3) {
	    Potoutr += 4*M_PI*G*rho(gi->router,halo)*(gi->router*gi->router)/(2-hp->beta); /* minus cancel */
	    }
	}
    Potr = (-1)*G*(gridr->Menc[i]/gridr->r[i]+Potoutr);
    gridr->Pot[i] = Potr;
    gridr->logPot[i] = log(-Potr);
    gridr->Potoutr[i] = Potoutr;
    for (i = (gi->Ngridr-2); i >= 0; i--) {
	Potoutr = gridr->Potoutr[i+1];
	if (gi->do_bulge == 1) {
	    Potoutr += integral(integrandPot,gridr->r[i],gridr->r[i+1],bulge);
	    }
	if (gi->do_halo == 1) {
	    Potoutr += integral(integrandPot,gridr->r[i],gridr->r[i+1],halo);
	    }
	Potr = (-1)*G*(gridr->Menc[i]/gridr->r[i]+Potoutr);
	gridr->Pot[i] = Potr;
	gridr->logPot[i] = log(-Potr);
	gridr->Potoutr[i] = Potoutr;
	}
    }

/*
** Routine for initialising griddf
*/

void initialise_griddf(const GI *gi, SI *si) {

    INT i, j, dj;
    DOUBLE fE;
    GRIDR *gridr;
    GRIDDF *griddf;

    gridr = gi->gridr;
    griddf = si->griddf;
    dj = (gi->Ngridr-1) / (gi->Ngriddf-1);
    for (i = (gi->Ngriddf-1), j = (gi->Ngridr-1); i >= 0; i--, j = j-dj) {
	fE = integraldf(j,gi,si);
	if (fE < 0) {
	    fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	    fprintf(stderr,"You have chosen parameters that lead to a negative and hence unphysical distribution function.\n");
	    usage();
	    }
	griddf->r[i] = gridr->r[j];
	griddf->logr[i] = gridr->logr[j];
	griddf->E[i] = gridr->Pot[j];
	griddf->logE[i] = gridr->logPot[j];
	griddf->fE[i] = fE;
	griddf->logfE[i] = log(fE);
	}
    }

/*
** Routine for initialising shell
*/

void initialise_shell(const GI *gi, SI *si) {

    INT i, Nshell;
    DOUBLE dlogr, logr;
    DOUBLE mass0;
    SHELL *shell;

    Nshell = si->Nshell;
    shell = malloc((Nshell+3)*sizeof(SHELL));
    assert(shell != NULL);
    if (Nshell > 0) {
	dlogr = (log(si->rso)-log(si->rsi))/Nshell;
	}
    else {
	dlogr = 0;
	}
    shell[0].rinner = exp(gi->gridr->logr[0]);
    for (i = 1; i < (Nshell+2); i++) {
	logr = log(si->rsi) + (i-1)*dlogr;
	shell[i].rinner = exp(logr);
	shell[i-1].router = shell[i].rinner;
	}
    shell[Nshell+2].rinner = exp(gi->gridr->logr[gi->Ngridr-1]);
    shell[Nshell+1].router = shell[Nshell+2].rinner;
    shell[Nshell+2].router = shell[Nshell+2].rinner;
    for (i = 0; i < (Nshell+3); i++) {
	shell[i].Menc = exp(lininterpolate(gi->Ngridr,gi->gridr->logr,si->logMenc,log(shell[i].rinner)));
	}
    mass0 = (shell[1].Menc-shell[0].Menc)/(2.0*(si->N0/2));
    for (i = 0; i < (Nshell+2); i++) {
	shell[i].massfac = pow(si->DRMmax,i);
	shell[i].mass = mass0*shell[i].massfac;
	shell[i].Mp = 0;
	shell[i].Ninitial = (INT)(((shell[i+1].Menc-shell[i].Menc)/(2.0*shell[i].mass))+0.5);
	shell[i].Nnosplit = shell[i].Ninitial;
	shell[i].Nnew = 0;
	shell[i].N = shell[i].Ninitial;
	shell[i].soft = si->soft0*pow(shell[i].massfac,1.0/(3.0-si->sp->gamma));
	shell[i].p = malloc(shell[i].N*sizeof(PARTICLE));
	assert(shell[i].p != NULL);
	}
    shell[Nshell+2].mass = 0;
    shell[Nshell+2].Mp = 0;
    shell[Nshell+2].Ninitial = 0;
    shell[Nshell+2].Nnosplit = 0;
    shell[Nshell+2].Nnew = 0;
    shell[Nshell+2].N = 0;
    shell[Nshell+2].soft = 0;
    shell[Nshell+2].p = malloc(0*sizeof(PARTICLE));
    assert(shell[Nshell+2].p != NULL);
    si->shell = shell;
    }
