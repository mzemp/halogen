/*
** check.c
**
** Checking routines for HALOGEN
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "definitions.h"
#include "usage.h"

void check_main_parameters_general(GI *gi) {

    if (strcmp(gi->outputname,"none") == 0) {
	fprintf(stderr,"You have not set a name for the output model.\n");
	usage();
	}
    if ((gi->Ngridr-1) % (gi->Ngriddf-1) != 0) {
	fprintf(stderr,"Bad choice of Ngridr (= %d) and Ngriddf (= %d)!\n",gi->Ngridr,gi->Ngriddf);
	fprintf(stderr,"These numbers have to fulfill the condition (Ngridr-1) mod (Ngriddf-1) == 0.\n");
	usage();
	}
    if (gi->coordinates == 1 && gi->positionsonly == 0) {
	gi->positionsonly = 1;
	fprintf(stderr,"You use ellipsoidal coordinates. Only the initialisation of positions is supported so far. Reset positionsonly value to 1.\n");
	}
    }

void check_main_parameters_system(const SI *si) {

    if (si->sp->alpha == -1) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"You have not set a value for alpha.\n");
	usage();
	}
    if (si->sp->beta == -1) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"You have not set a value for beta.\n");
	usage();
	}
    if (si->sp->gamma == -1) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"You have not set a value for gamma.\n");
	usage();
	}
    if (si->sp->gamma >= 3) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"You have chosen gamma = "OFD1".\n",si->sp->gamma);
	fprintf(stderr,"This means your cumulative mass function is diverging at the centre.\n");
	fprintf(stderr,"Use a smaller value for gamma.\n");
	usage();
	}
    if (si->sp->M == -1) {	
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"You have not set a value for the mass M.\n");
	usage();
	}
    if (si->sp->rba > 1) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"You have set rba > 1.\n");
	usage();
	}
    if (si->sp->rca > 1) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"You have set rca > 1.\n");
	usage();
	}
    if (si->sp->rca > si->sp->rba) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"You have set rca > rba.\n");
	usage();
	}
    if (si->N0 == -1) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"You have not set a value for the number of paricles N0.\n");
	usage();
	}
    if (si->soft0 == -1) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"You have not set a value for the softening of the particles soft0.\n");
	usage();
	}
    }

void check_more_parameters_system(const GI *gi, const SI *si) {

    if (si->rsi < 0) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"Your choice of rsi (= "OFD1" LU) is negative!\n",si->rsi);
	fprintf(stderr,"Please choose rsi positive.\n");
	usage();
	}
    if (si->rso < 0) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"Your choice of rso (= "OFD1" LU) is negative!\n",si->rso);
	fprintf(stderr,"Please choose rso positive.\n");
	usage();
	}
    if (si->rmor < 0) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"Your choice of rmor (= "OFD1" LU) is negative!\n",si->rmor);
	fprintf(stderr,"Please choose rmor positive.\n");
	usage();
	}
    if (si->rsi < 10*gi->rinner) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"Your choice of rsi (= "OFD1" LU) is smaller than 10*rinner (= "OFD1" LU).\n",si->rsi,10*gi->rinner);
	fprintf(stderr,"Please choose rsi larger.\n");
	usage();
	}
    if (si->rso > gi->router) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"Your choice of rso (= "OFD1" LU) is larger than router (= "OFD1" LU).\n",si->rso,gi->router);
	fprintf(stderr,"Please choose rso smaller.\n");
	usage();
	}
    if (si->rso < si->rsi) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"Your choice of rso (= "OFD1" LU) is smaller than rsi (= "OFD1" LU).\n",si->rso,si->rsi);
	fprintf(stderr,"Please choose different values.\n");
	usage();
	}
    if (si->Nshell < 0) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"You have chosen Nshell = "OFI1".\n",si->Nshell);
	fprintf(stderr,"Please choose a positive value for Nshell.\n");
	usage();
	}
    if (si->Nshell == 0 && si->rso != si->rsi) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"You have chosen Nshell = "OFI1".\n",si->Nshell);
	fprintf(stderr,"You have chosen rsi (= "OFD1" LU) not equal to rso (= "OFD1" LU)!\n",si->rsi,si->rso);
	fprintf(stderr,"Please set them equal.\n");
	usage();
	}
    if (si->Nshell > 0 && si->rso == si->rsi) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"You have chosen Nshell = "OFI1".\n",si->Nshell);
	fprintf(stderr,"You have chosen rsi (= "OFD1" LU) equal to rso (= "OFD1" LU)!\n",si->rsi,si->rso);
	fprintf(stderr,"Please set them unequal.\n");
	usage();
	}
    if (si->DRMmax == 0) {
	fprintf(stderr,"Missing or bad parameter for the %s.\n",si->systemname);
	fprintf(stderr,"You have chosen DRMmax = "OFI1".\n",si->DRMmax);
	fprintf(stderr,"Please choose a different value.\n");
	usage();
	}
    if ((si->DRMmax != 1) && (si->rmor < si->rsi)) {
	fprintf(stderr,"Warning for the %s!\n",si->systemname);
	fprintf(stderr,"You have chosen rmor (= "OFD1" LU) smaller than rsi (= "OFD1" LU)!\n",si->rmor,si->rsi);
	fprintf(stderr,"Are you sure about that? All particles with rperi < rsi will be\n");
	fprintf(stderr,"refined by default if 0 LU < rmor < rsi.\n");
	fprintf(stderr,"Info: if you choose to set rmor = 0 LU then no refinement will be done.\n");
	}
    if ((si->sp->beta <= 3) && (si->sp->rvir < si->sp->rcutoff)) {
	fprintf(stderr,"Warning for the %s!\n",si->systemname);
	fprintf(stderr,"You have chosen rvir (= "OFD1" LU) smaller than rcutoff (= "OFD1" LU)!\n",si->sp->rvir,si->sp->rcutoff);
	fprintf(stderr,"Are you sure about that?\n");
	}
    }
