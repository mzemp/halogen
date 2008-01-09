/* 
** routines.c 
**
** General routines for HALOGEN
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
#include "routines.h"
#include "write.h"
#include "usage.h"

/*
** Routine for calculatein parameters for general info structure
*/

void calculate_parameters_general_info(GI *gi) {

    gi->rhocritz = 3/(8*M_PI*G)*(pow((100*gi->h0*Ecosmo(gi)*VelConvertFac/1000),2));
    if (gi->Deltavirz == -1) {
	gi->OmegaMz = gi->OmegaM0*pow((1+gi->z),3)/pow(Ecosmo(gi),2);
	if (gi->OmegaK0 == 0) {
	    gi->Deltavirz = 178*(pow(gi->OmegaMz,0.45));
	    }
	else if (gi->OmegaL0 == 0) {
	    gi->Deltavirz = 178*(pow(gi->OmegaMz,0.3));
	    }
	else {
	    fprintf(stderr,"HALOGEN can't calculate a value for Deltavirz for that choice of cosmological parameters.\n");
	    fprintf(stderr,"Set a value vor Deltavirz by hand or choose a different cosmology.\n");
	    usage();
	    }
	}
    }

/*
** Routine for calculating parameters for system
*/

void calculate_parameters_system(const GI *gi, SI *si) {

    DOUBLE I_M;

    if (si->sp->beta > 3) {
	/*
	** Finite mass models
	*/
	if (si->sp->rs == -1) {
	    fprintf(stderr,"Missing or bad parameter for %s.\n",si->systemname);
	    fprintf(stderr,"For finite mass models you have to set a value for the scale radius rs.\n");
	    usage();
	    }
	if (si->sp->cvir != -1) {
	    fprintf(stderr,"Warning for %s!\n",si->systemname);
	    fprintf(stderr,"For finite mass models the virial concentration cvir is calculated selfconsistently.\n");
	    fprintf(stderr,"Hence, your input for the virial concentration cvir (= "OFD1") was ignored.\n",si->sp->cvir);
	    }
	if (si->sp->rcutoff != -1) {
	    fprintf(stderr,"Warning for %s!\n",si->systemname);
	    fprintf(stderr,"For finite mass models the cutoff radius rcutoff is not needed!\n");
	    fprintf(stderr,"Hence, your input for the cutoff radius rcutoff (= "OFD1" kpc) was ignored.\n",si->sp->rcutoff);
	    }
	I_M = exp(lgamma((si->sp->beta-3)/si->sp->alpha));
	I_M *= exp(lgamma((3-si->sp->gamma)/si->sp->alpha));
	I_M /= (si->sp->alpha*exp(lgamma((si->sp->beta-si->sp->gamma)/si->sp->alpha)));
	si->sp->rho0 = si->sp->M/(4*M_PI*(si->sp->rs*si->sp->rs*si->sp->rs)*I_M);
	si->sp->rcutoff = SBI;
	si->sp->rdecay = 0;
	si->sp->delta = 0;
	}
    else {
	/*
	** Cutoff models
	*/
	if (((si->sp->cvir != -1) && (si->sp->rs != -1) && (si->sp->rcutoff != -1)) ||
	    ((si->sp->cvir == -1) && ((si->sp->rs == -1) || (si->sp->rcutoff == -1)))) {
	    fprintf(stderr,"HALOGEN is confused!\n");
	    fprintf(stderr,"Missing or bad parameter for %s.\n",si->systemname);
	    fprintf(stderr,"Specify either just a value for the virial concentration cvir or give values\n");
	    fprintf(stderr,"for the scale radius rs and cutoff radius rcutoff.\n");
	    usage();
	    }
	if (si->sp->cvir != -1) {
	    /*
	    ** Models with virial overdensity
	    */
	    if (si->sp->rs != -1) {
		fprintf(stderr,"Warning for %s!\n",si->systemname);
		fprintf(stderr,"If you set a virial concentration cvir the scale radius rs is calculated selfconsistently.\n");
		fprintf(stderr,"Hence, your input for the scale radius rs (= "OFD1" kpc) was ignored.\n",si->sp->rs);
		}
	    if (si->sp->rcutoff != -1) {
		fprintf(stderr,"Warning for %s!\n",si->systemname);
		fprintf(stderr,"If you set a virial concentration cvir the cutoff radius rcutoff is set to the virial radius rvir.\n");
		fprintf(stderr,"Hence, your input for the cutoff radius rcutoff (= "OFD1" kpc) was ignored.\n",si->sp->rcutoff);
		}
	    si->sp->rvir = pow(3*si->sp->M/(4*M_PI*gi->rhocritz*gi->Deltavirz),1.0/3.0);
	    si->sp->vvir = sqrt(G*si->sp->M/si->sp->rvir);
	    si->sp->rs = si->sp->rvir/si->sp->cvir;
	    si->sp->rcutoff = si->sp->rvir;
	    }
	I_M = pow(1e-6,3-si->sp->gamma)/(3-si->sp->gamma); /* approximate inner integral */
	I_M += integral(integrandIM,1e-6,si->sp->rcutoff/si->sp->rs,si);
	si->sp->rho0 = si->sp->M/(4*M_PI*(si->sp->rs*si->sp->rs*si->sp->rs)*I_M);
	si->sp->rdecay = CutoffFac*si->sp->rcutoff;
	si->sp->delta = si->sp->rcutoff/si->sp->rdecay + dlrhodlr(si->sp->rcutoff,si);
	}
    }

/*
** Routine for calculating virial stuff for finite mass models and cutoff models
*/

void calculate_virial_stuff(const GI *gi, SI *si) {

    if (si->sp->beta > 3) {
	/*
	** Finite mass models
	*/
	si->sp->rvir = exp(lininterpolate(gi->Ngridr,gi->gridr->logrhoencHalo,gi->gridr->logr,log(gi->Deltavirz*gi->rhocritz)));
	si->sp->cvir = si->sp->rvir/si->sp->rs;
	si->sp->vvir = sqrt(G*si->sp->M/si->sp->rvir);
	}
    else {
	/*
	** Cutoff models
	*/
	if (si->sp->rvir == -1) {
	    si->sp->rvir = exp(lininterpolate(gi->Ngridr,gi->gridr->logrhoencHalo,gi->gridr->logr,log(gi->Deltavirz*gi->rhocritz)));
	    si->sp->cvir = si->sp->rvir/si->sp->rs;
	    si->sp->vvir = sqrt(G*si->sp->M/si->sp->rvir);
	    }
	}
    }

/*
** Routine for setting remaining parameters
*/

void set_remaining_parameters(const GI *gi, SI *si) {

    /*
    ** Set rsi
    */
    if (si->rsi_in_rs_units == 1) {
	si->rsi *= si->sp->rs;
	}
    else if (si->rsi_in_rvir_units == 1) {
	si->rsi *= si->sp->rvir;
	}
    else if (si->rsi_in_rcutoff_units == 1) {
	si->rsi *= si->sp->rcutoff;
	}
    if (si->rsi == -1) {
	/*
	** Default values
	*/
	if (si->sp->beta > 3) {
	    si->rsi = gi->router;
	    }
	else {
	    si->rsi = si->sp->rcutoff;
	    }
	}
    /*
    ** Set rso
    */
    if (si->rso_in_rs_units == 1) {
	si->rso *= si->sp->rs;
	}
    else if (si->rso_in_rvir_units == 1) {
	si->rso *= si->sp->rvir;
	}
    else if (si->rso_in_rcutoff_units == 1) {
	si->rso *= si->sp->rcutoff;
	}
    if (si->rso == -1) {
	/*
	** Default values
	*/
	if (si->sp->beta > 3) {
	    si->rso = gi->router;
	    }
	else {
	    si->rso = si->sp->rcutoff;
	    }
	}
    /*
    ** Ser rmor
    */
    if (si->rmor_in_rs_units == 1) {
	si->rmor *= si->sp->rs;
	}
    else if (si->rmor_in_rvir_units == 1) {
	si->rmor *= si->sp->rvir;
	}
    else if (si->rmor_in_rcutoff_units == 1) {
	si->rmor *= si->sp->rcutoff;
	}
    if (si->rmor == -1) {
	/*
	** Default value
	*/
	si->rmor = 0;
	}
    /*
    ** Set Ismor
    */
    if (si->Ismor == -1) {
	si->Ismor = si->Nshell+1;
	}
    }

/*
** Routine for setting position of particles
*/

void set_positions(const GI *gi, SI *si) {
    
    INT i, j, N;
    DOUBLE Mrand, logMrand, Mmin, Mmax;
    DOUBLE rrand, logrrand;
    DOUBLE costheta, sintheta, phi, cosphi, sinphi;
    PARTICLE *p;

    for (j = 0; j < (si->Nshell+2); j++) {
	N = si->shell[j].N;
	p = si->shell[j].p;
	Mmin = si->shell[j].Menc;
	Mmax = si->shell[j+1].Menc;
	for (i = 0; i < N; i++) {
	    Mrand = Mmin + rand01()*(Mmax - Mmin);
	    logMrand = log(Mrand);
	    logrrand = lininterpolate(gi->Ngridr,gi->gridr->logMencHalo,gi->gridr->logr,logMrand);
	    rrand = exp(logrrand);
	    costheta = 2.0*rand01() - 1.0;
	    sintheta = sqrt(1-costheta*costheta);
	    phi = rand01()*2.0*M_PI;
	    cosphi = cos(phi);
	    sinphi = sin(phi);
	    p[i].r[0] = rrand;
	    p[i].r[1] = rrand*sintheta*cosphi;
	    p[i].r[2] = rrand*sintheta*sinphi;
	    p[i].r[3] = rrand*costheta;
	    }
	}
    }

/*
** Routine for setting velocities of particles 
*/

void set_velocities(const GI *gi, SI *si) {
    
    INT i, j, N, isplit;
    DOUBLE r, Potr, Esplit, Erand;
    DOUBLE fEmax, fEsplit, fErand, fEcheck;
    DOUBLE Compmax, Compsplit, Comprand;
    DOUBLE vesc, vsplit, vrand;
    DOUBLE costheta, sintheta, phi, cosphi, sinphi;
    GRIDR *gridr;
    GRIDDF *griddf;
    PARTICLE *p;

    gridr = gi->gridr;
    griddf = si->griddf;
    for (j = 0; j < (si->Nshell+2); j++) {
	N = si->shell[j].N;
	p = si->shell[j].p;
	for (i = 0; i < N; i++) {
	    r = p[i].r[0];
	    Potr = Pot(r,gi);
	    vesc = vescape(r,gi);
	    fEmax = f2(r,gi,si);
	    if (si->dfsf == 0 || si->dfsf == 1) {
		/*
		** No splitting
		*/
		isplit = gi->Ngriddf-1;
		}
	    else {
		isplit = locate(gi->Ngriddf,griddf->fE,si->dfsf*fEmax);
		}
	    if (griddf->fE[isplit] > fEmax) {
		isplit += 1;
		}
	    assert(isplit < gi->Ngriddf);
	    Esplit = griddf->E[isplit];
	    vsplit = sqrt(2*(Esplit-Potr));
	    fEsplit = griddf->fE[isplit];
	    Compsplit = fEmax*vsplit*vsplit*vsplit/3.0;
	    Compmax = Compsplit + (fEsplit/3.0)*(vesc*vesc*vesc-vsplit*vsplit*vsplit);
	    vrand = 0;
	    Erand = 0;
	    fErand = 0;
	    fEcheck = 1;
	    while (fEcheck > fErand) {
		/*
		** Get random value of comparison function
		** and calculate the according velocity vrand
		*/
		Comprand = rand01()*Compmax;
		if (Comprand <= Compsplit) {
		    vrand = (3.0/fEmax)*Comprand;
		    vrand = pow(vrand,1.0/3.0);
		    }
		else {
		    vrand = Comprand - Compsplit + (fEsplit/3.0)*vsplit*vsplit*vsplit;
		    vrand *= 3.0/fEsplit;
		    vrand = pow(vrand,1.0/3.0);
		    }
		/*
		** Calculate the value of the distribution
		** function at vrand
		*/
		Erand = 0.5*vrand*vrand + Potr;
		fErand = f1(Erand,gi,si);
		/*
		** Generate random value fEcheck for
		** acceptance/rejection criterion
		*/
		if (Comprand <= Compsplit) {
		    fEcheck = rand01()*fEmax;
		    }
		else {
		    fEcheck = rand01()*fEsplit;
		    }
		}
	    costheta = 2.0*rand01() - 1.0;
	    sintheta = sqrt(1-costheta*costheta);
	    phi = rand01()*2.0*M_PI;
	    cosphi = cos(phi);
	    sinphi = sin(phi);
	    p[i].v[0] = vrand;
	    p[i].v[1] = vrand*sintheta*cosphi;
	    p[i].v[2] = vrand*sintheta*sinphi;
	    p[i].v[3] = vrand*costheta;
	    }
	}
    }

/*
** Routine for setting velocities of particles to zero 
*/

void set_velocities_zero(SI *si) {
    
    INT i, j, k, N;
    PARTICLE *p;

    for (j = 0; j < (si->Nshell+2); j++) {
	N = si->shell[j].N;
	p = si->shell[j].p;
	for (i = 0; i < N; i++) {
	    for (k = 0; k < 4; k++) {
		p[i].v[k] = 0;
		}
	    }
	}
    }

/*
** Routine for setting the remaining attributes
*/

void set_attributes(const GI *gi, SI *si) {

    INT i, j, N;
    DOUBLE Ekin, Epot, Etot;
    DOUBLE mass, soft;
    GRIDR *gridr;
    PARTICLE *p;
    
    gridr = gi->gridr;
    for (j = 0; j < (si->Nshell+2); j++) {
	N = si->shell[j].N;
	p = si->shell[j].p;
	mass = si->shell[j].mass;
	soft = si->shell[j].soft;
	for (i = 0; i < N; i++) {
	    Ekin = p[i].v[0]*p[i].v[0]/2.0;
	    Epot = Pot(p[i].r[0],gi);
	    Etot = Ekin + Epot;
	    p[i].L[1] = p[i].r[2]*p[i].v[3] - p[i].r[3]*p[i].v[2];
	    p[i].L[2] = p[i].r[3]*p[i].v[1] - p[i].r[1]*p[i].v[3];
	    p[i].L[3] = p[i].r[1]*p[i].v[2] - p[i].r[2]*p[i].v[1];
	    p[i].L[0] = sqrt(p[i].L[1]*p[i].L[1]+p[i].L[2]*p[i].L[2]+p[i].L[3]*p[i].L[3]);
	    p[i].mass = mass;
	    p[i].soft = soft;
	    p[i].Ekin = Ekin;
	    p[i].Epot = Epot;
	    p[i].Etot = Etot;
	    }
	}
    }

/*
** Routines for refining
*/

void refine(const GI *gi, SI *si) {

    INT i, j, k, N, Nnew, Nrealloc, splitfac;
    INT index[2];
    DOUBLE Etot, L2;
    DOUBLE dgx, dgy, dx, m;
    DOUBLE massnew, softnew;
    DOUBLE costheta, sintheta, phi, cosphi, sinphi;
    DOUBLE ScalarProd;
    DOUBLE erad1, erad2, erad3;
    DOUBLE ephi1, ephi2, ephi3;
    DOUBLE etheta1, etheta2, etheta3;
    DOUBLE v1, v2, v3;
    DOUBLE vrad, vtan;
    DOUBLE rootfunc[gi->Ngridr];
    GRIDR *gridr;
    PARTICLE *p;

    gridr = gi->gridr;
    if((si->DRMmax > 1) && (si->rmor > 0)) {
	for (k = 1; k < (si->Ismor+1); k++) {
	    N = si->shell[k].N;
	    p = si->shell[k].p;
	    Nnew = 0;
	    Nrealloc = 2*N;
	    p = realloc(p,Nrealloc*sizeof(PARTICLE));
	    assert(p != NULL);
	    for (i = 0; i < N; i++) {
		Etot = p[i].Etot;
		L2 = p[i].L[0]*p[i].L[0];
		for(j = 0; j < gi->Ngridr; j++) {
		    rootfunc[j] = 1.0/(gridr->r[j]*gridr->r[j]) + 2.0*(gridr->Pot[j]-Etot)/L2;
		    }
		if(rootfunc[0] < 0) {
		    p[i].rperi = 0;
		    index[0] = gi->Ngridr;
		    index[1] = gi->Ngridr;
		    searchroot(gi,index,rootfunc);
		    if(index[0] == gi->Ngridr) {
			fprintf(stderr,"Something strange happened! You should never get here!\n");
			fprintf(stderr,"N = "OFI1" i = "OFI1" r = "OFD3" E = "OFD3" L2 = "OFD3"\n",N,i,p[i].r[0],Etot,L2);
			}
		    dgx = rootfunc[index[0]+1] - rootfunc[index[0]];
		    dgy = gridr->logr[index[0]+1] - gridr->logr[index[0]];
		    m = dgy/dgx;
		    dx = 0.0 - rootfunc[index[0]];
		    p[i].rapo = exp(gridr->logr[index[0]] + dx*m);
		    }
		else {
		    index[0] = gi->Ngridr;
		    index[1] = gi->Ngridr;
		    searchroot(gi,index,rootfunc);
		    if(index[0] == gi->Ngridr || index[1] == gi->Ngridr) {
			searchmin(gi,index,rootfunc);
			p[i].rperi = exp(gridr->logr[index[0]]);
			p[i].rapo = p[i].rperi;
			}
		    else {
			dgx = rootfunc[index[0]+1] - rootfunc[index[0]];
			dgy = gridr->logr[index[0]+1] - gridr->logr[index[0]];
			m = dgy/dgx;
			dx = 0.0 - rootfunc[index[0]];
			p[i].rperi = exp(gridr->logr[index[0]] + dx*m);
			dgx = rootfunc[index[1]+1] - rootfunc[index[1]];
			dgy = gridr->logr[index[1]+1] - gridr->logr[index[1]];
			m = dgy/dgx;
			dx = 0.0 - rootfunc[index[1]];
			p[i].rapo = exp(gridr->logr[index[1]] + dx*m);
			}
		    }
		p[i].ecc = (p[i].rapo - p[i].rperi)/(p[i].rapo + p[i].rperi);
		splitfac = split(k,p[i].rperi,si);
		if (splitfac > 1) {
		    si->shell[k].Nnosplit--;
		    Nnew = Nnew + (splitfac - 1);
		    if(N+Nnew > Nrealloc) {
			Nrealloc = 2*(N+Nnew);
			p = realloc(p,Nrealloc*sizeof(PARTICLE));
			assert(p != NULL);
			}
		    massnew = p[i].mass/splitfac;
		    softnew = p[i].soft/pow(splitfac,1.0/(3.0-si->sp->gamma));
		    p[i].mass = massnew;
		    p[i].soft = softnew;
		    ScalarProd = p[i].r[1]*p[i].v[1] + p[i].r[2]*p[i].v[2] + p[i].r[3]*p[i].v[3];
		    vrad = ScalarProd/p[i].r[0];
		    vtan = p[i].L[0]/p[i].r[0];
		    for(j = N + Nnew - (splitfac - 1); j < (N + Nnew); j++) {
			/*
			** Create orthonormal basis in spherical coordinates with random direction
			*/
			costheta = 2.0*rand01() - 1.0;
			sintheta = sqrt(1-costheta*costheta);
			phi = rand01()*2.0*M_PI;
			cosphi = cos(phi);
			sinphi = sin(phi);
			erad1 = sintheta*cosphi;
			erad2 = sintheta*sinphi;
			erad3 = costheta;
			ephi1 = -sinphi;
			ephi2 = cosphi;
			ephi3 = 0;
			etheta1 = -costheta*cosphi;
			etheta2 = -costheta*sinphi;
			etheta3 = sintheta;
			/*
			** Calculate new position
			*/
			p[j].r[0] = p[i].r[0];
			p[j].r[1] = p[i].r[0]*erad1;
			p[j].r[2] = p[i].r[0]*erad2;
			p[j].r[3] = p[i].r[0]*erad3;
			/*
			** Generate random direction in tangential plane
			*/
			phi = rand01()*2.0*M_PI;
			cosphi = cos(phi);
			sinphi = sin(phi);
			v1 = cosphi*ephi1 + sinphi*etheta1;
			v2 = cosphi*ephi2 + sinphi*etheta2;
			v3 = cosphi*ephi3 + sinphi*etheta3;
			/*
			** Calculate new velocity
			*/
			p[j].v[0] = p[i].v[0];
			p[j].v[1] = vrad*erad1 + vtan*v1;
			p[j].v[2] = vrad*erad2 + vtan*v2;
			p[j].v[3] = vrad*erad3 + vtan*v3;
			/*
			** Calculate angular momentum
			*/
			p[j].L[1] = p[j].r[2]*p[j].v[3] - p[j].r[3]*p[j].v[2];
			p[j].L[2] = p[j].r[3]*p[j].v[1] - p[j].r[1]*p[j].v[3];
			p[j].L[3] = p[j].r[1]*p[j].v[2] - p[j].r[2]*p[j].v[1];
			p[j].L[0] = sqrt(p[j].L[1]*p[j].L[1]+p[j].L[2]*p[j].L[2]+p[j].L[3]*p[j].L[3]);
			/*
			** Copy rest of characteristics
			*/
			p[j].mass = massnew;
			p[j].soft = softnew;
			p[j].Ekin = p[i].Ekin;
			p[j].Epot = p[i].Epot;
			p[j].Etot = p[i].Etot;
			p[j].rperi = p[i].rperi;
			p[j].rapo = p[i].rapo;
			}
		    }
		}
	    si->shell[k].Nnew = Nnew;
	    si->shell[k].N = N+Nnew;
	    si->shell[k].p = p;
	    }
	}
    }

/*
** Routine for searching roots
*/

void searchroot(const GI *gi, INT *index, DOUBLE *grid) {

    INT i, j;

    j = 0;
    for(i = 0; i < (gi->Ngridr-2); i++) {
	if(grid[i]*grid[i+1] < 0) {
	    index[j] = i;
	    j ++;
	    if(j == 2) {
		break;
		}
	    }
	}
    }

/*
** Routine for searching minimum in case no root is found
*/

void searchmin(const GI *gi, INT *index, DOUBLE *grid) {

    INT i;
    DOUBLE min;

    min = grid[0];
    for(i = 1; i < (gi->Ngridr-1); i++) {
	if(grid[i] < min) {
	    min = grid[i];
	    index[0] = i;
	    }
	}
    }

/*
** Routine for doubling particles with mirror halo
*/

void double_particles(SI *si) {

    INT i, j, N;
    PARTICLE *p;

    for(j = 0; j < (si->Nshell+2); j++) {
	N = si->shell[j].N;
	p = si->shell[j].p;
	p = realloc(p,2*N*sizeof(PARTICLE));
	if (N > 0) {
	    assert(p != NULL);
	    }
	for (i = 0; i < N; i++) {
	    p[N+i].r[0] = p[i].r[0];
	    p[N+i].r[1] = -p[i].r[1];
	    p[N+i].r[2] = -p[i].r[2];
	    p[N+i].r[3] = -p[i].r[3];
	    p[N+i].v[0] = p[i].v[0];
	    p[N+i].v[1] = -p[i].v[1];
	    p[N+i].v[2] = -p[i].v[2];
	    p[N+i].v[3] = -p[i].v[3];
	    p[N+i].L[0] = p[i].L[0];
	    p[N+i].L[1] = p[i].L[1];
	    p[N+i].L[2] = p[i].L[2];
	    p[N+i].L[3] = p[i].L[3];
	    p[N+i].mass = p[i].mass; 
	    p[N+i].soft = p[i].soft;
	    p[N+i].Ekin = p[i].Ekin;
	    p[N+i].Epot = p[i].Epot;
	    p[N+i].Etot = p[i].Etot;
	    p[N+i].rperi = p[i].rperi;
	    p[N+i].rapo = p[i].rapo;
	    p[N+i].ecc = p[i].ecc;
	    }
	si->shell[j].Ninitial = 2*si->shell[j].Ninitial;
	si->shell[j].Nnosplit = 2*si->shell[j].Nnosplit;
	si->shell[j].Nnew = 2*si->shell[j].Nnew;
	si->shell[j].N = 2*si->shell[j].N;
	si->shell[j].p = p;
	}
    }

/*
** Routine for calculating center of mass position and velocity, and angular momentum
*/

void calculate_stuff(GI *gi, PARTICLE *bh, SI *halo) {

    INT i, j, N;
    DOUBLE temp;
    STUFF *stuff;
    PARTICLE *p;

    stuff = gi->stuff;
    stuff->Ntot = 0;
    stuff->Ninitialtot = 0;
    stuff->Nnosplittot = 0;
    stuff->Nnewtot = 0;
    stuff->Mp = 0;
    stuff->Ekin = 0;
    stuff->Epot = 0;
    stuff->Nfemm = 0;
    stuff->Nfesm = 0;
    for(i = 0; i < 4; i++) {
	stuff->Cr[i] = 0;
	stuff->Cv[i] = 0;
	stuff->Ltot[i] = 0;
	}
    if (bh->mass > 0) {
	stuff->Ntot++;
	stuff->Mp += bh->mass;
	}
    for(j = 0; j < (halo->Nshell+2); j++) {
	N = halo->shell[j].N;
	p = halo->shell[j].p;
	stuff->Ntot += halo->shell[j].N;
	stuff->Ninitialtot += halo->shell[j].Ninitial;
	stuff->Nnosplittot += halo->shell[j].Nnosplit;
	stuff->Nnewtot += halo->shell[j].Nnew;
	for (i = 0; i < N; i++) {
	    if((j == 0) && (p[i].r[0] < halo->rimp)) {
		halo->rimp = p[i].r[0];
		}
	    temp = p[i].mass;
	    halo->shell[j].Mp += temp;
	    stuff->Mp += temp;
	    stuff->Cr[1] += temp*p[i].r[1];
	    stuff->Cr[2] += temp*p[i].r[2];
	    stuff->Cr[3] += temp*p[i].r[3];
	    stuff->Cv[1] += temp*p[i].v[1];
	    stuff->Cv[2] += temp*p[i].v[2];
	    stuff->Cv[3] += temp*p[i].v[3];
	    stuff->Ltot[1] += temp*p[i].L[1];
	    stuff->Ltot[2] += temp*p[i].L[2];
	    stuff->Ltot[3] += temp*p[i].L[3];
	    stuff->Ekin += temp*p[i].Ekin;
	    stuff->Epot += temp*p[i].Epot;
	    temp = (2*M_PI)/(0.03*Tdyn(p[i].r[0],gi));
	    stuff->Nfemm += temp;
	    stuff->Nfesm += temp*p[i].mass/halo->shell[0].mass;
	    }
	}
    for(i = 1; i < 4; i++) {
	stuff->Cr[i] = stuff->Cr[i]/stuff->Mp;
	stuff->Cv[i] = stuff->Cv[i]/stuff->Mp;
	}
    stuff->Cr[0] = sqrt(stuff->Cr[1]*stuff->Cr[1]+stuff->Cr[2]*stuff->Cr[2]+stuff->Cr[3]*stuff->Cr[3]);
    stuff->Cv[0] = sqrt(stuff->Cv[1]*stuff->Cv[1]+stuff->Cv[2]*stuff->Cv[2]+stuff->Cv[3]*stuff->Cv[3]);
    stuff->Ltot[0] = sqrt(stuff->Ltot[1]*stuff->Ltot[1]+stuff->Ltot[2]*stuff->Ltot[2]+stuff->Ltot[3]*stuff->Ltot[3]);
    stuff->Epot = stuff->Epot/2.0;
    stuff->Etot = stuff->Ekin + stuff->Epot;
    if (gi->gridr->eqrvcmax[0] < 0) {
	halo->sp->rvcmax = exp(lininterpolate(gi->Ngridr,gi->gridr->eqrvcmax,gi->gridr->logr,0.0));
	halo->sp->vcmax = sqrt(G*Menc(halo->sp->rvcmax,gi)/halo->sp->rvcmax);
	}
    halo->sp->rhalf = exp(lininterpolate(gi->Ngridr,gi->gridr->Menc,gi->gridr->logr,halo->sp->M/2.0));
    halo->r1 = pow(((3.0-halo->sp->gamma)*halo->shell[0].mass)/(4*M_PI*halo->sp->rho0*pow(halo->sp->rs,halo->sp->gamma)),1.0/(3.0-halo->sp->gamma));
    halo->r100 = pow(((3.0-halo->sp->gamma)*100*halo->shell[0].mass)/(4*M_PI*halo->sp->rho0*pow(halo->sp->rs,halo->sp->gamma)),1.0/(3.0-halo->sp->gamma));
    }

