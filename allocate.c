/*
** allocate.c
**
** Some allocation routines for HALOGEN
*/

#include <string.h>
#include <malloc.h>
#include <assert.h>
#include "definitions.h"

/*
** Routine for allocating memory for general info structure
*/

void allocate_general(GI *gi) {

	gi->gridr = malloc(sizeof(GRIDR));
	assert(gi->gridr != NULL);
	gi->gridr->r = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->r != NULL);
	gi->gridr->logr = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->logr != NULL);
	gi->gridr->rho = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->rho != NULL);
	gi->gridr->logrho = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->logrho != NULL);
	gi->gridr->rhoBulge = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->rhoBulge != NULL);
	gi->gridr->logrhoBulge = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->logrhoBulge != NULL);
	gi->gridr->rhoHalo = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->rhoHalo != NULL);
	gi->gridr->logrhoHalo = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->logrhoHalo != NULL);
	gi->gridr->rhoenc = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->rhoenc != NULL);
	gi->gridr->logrhoenc = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->logrhoenc != NULL);
	gi->gridr->rhoencBulge = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->rhoencBulge != NULL);
	gi->gridr->logrhoencBulge = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->logrhoencBulge != NULL);
	gi->gridr->rhoencHalo = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->rhoencHalo != NULL);
	gi->gridr->logrhoencHalo = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->logrhoencHalo != NULL);
	gi->gridr->Menc = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->Menc != NULL);
	gi->gridr->logMenc = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->logMenc != NULL);
	gi->gridr->MencBulge = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->MencBulge != NULL);
	gi->gridr->logMencBulge = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->logMencBulge != NULL);
	gi->gridr->MencHalo = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->MencHalo != NULL);
	gi->gridr->logMencHalo = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->logMencHalo != NULL);
	gi->gridr->Pot = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->Pot != NULL);
	gi->gridr->logPot = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->logPot != NULL);
	gi->gridr->Potoutr = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->Potoutr != NULL);
	gi->gridr->sigma = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->sigma != NULL);
	gi->gridr->logsigma = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->logsigma != NULL);
	gi->gridr->sigmaintoutr = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->sigmaintoutr != NULL);
	gi->gridr->eqrvcmax = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->eqrvcmax != NULL);
	gi->gridr->eqrvcmaxBulge = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->eqrvcmaxBulge != NULL);
	gi->gridr->eqrvcmaxHalo = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(gi->gridr->eqrvcmaxHalo != NULL);
	}

/*
** Routine for allocating memory for a system
*/

void allocate_system(const GI *gi, SI *si) {

	si->shell = malloc(sizeof(SHELL));
	assert(si->shell != NULL);
	si->griddf = malloc(sizeof(GRIDDF));
	assert(si->griddf != NULL);
	si->griddf->r = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(si->griddf->r != NULL);
	si->griddf->logr = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(si->griddf->logr != NULL);
	si->griddf->E = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(si->griddf->E != NULL);
	si->griddf->logE = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(si->griddf->logE != NULL);
	si->griddf->fE = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(si->griddf->fE != NULL);
	si->griddf->logfE = malloc(gi->Ngridr*sizeof(DOUBLE));
	assert(si->griddf->logfE != NULL);
	if (strcmp(si->systemname,"bulge") == 0) {
		si->logrhoenc = gi->gridr->logrhoencBulge;
		si->logMenc = gi->gridr->logMencBulge;
		si->eqrvcmax = gi->gridr->eqrvcmaxBulge;
		}
	else if (strcmp(si->systemname,"halo") == 0) {
		si->logrhoenc = gi->gridr->logrhoencHalo;
		si->logMenc = gi->gridr->logMencHalo;
		si->eqrvcmax = gi->gridr->eqrvcmaxHalo;
		}
	}
