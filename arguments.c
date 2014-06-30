/*
** arguments.c
**
** Routine for processing input arguments
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "definitions.h"
#include "usage.h"

/*
** Routine for processing input arguments
*/

void process_arguments(int argc, char **argv, GI *gi, PARTICLE *bh, SI *bulge, SI *halo) {

	INT i;
	SI *si;

	si = NULL;
	i = 1;
	while (i < argc) {
		/*
		** Help
		*/
		if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
			usage();
			}
		/*
		** Version
		*/
		if (strcmp(argv[i],"-version") == 0) {
			fprintf(stderr,"HALOGEN (%s) by Marcel Zemp\n",VERSION);
			exit(1);
			}
		/*
		** Model parameters
		*/
		if (strcmp(argv[i],"-bulge") == 0) {
			si = bulge;
			gi->do_bulge = 1;
			i++;
			}
		else if (strcmp(argv[i],"-halo") == 0) {
			si = halo;
			gi->do_halo = 1;
			i++;
			}
		if (si == NULL) {
			fprintf(stderr,"You have no system specified!\n");
			usage();
			}
		if (strcmp(argv[i],"-a") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->alpha = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-b") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->beta = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-c") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->gamma = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-M") == 0) {
			i++;
			if (i >= argc) usage();
			if (strcmp("Mo",argv[i]) == 0) {
				si->sp->M = 1;
				}
			else if (strcmp("MU",argv[i]) == 0) {
				si->sp->M = MU;
				}
			else {
				si->sp->M = atof(argv[i]);
				}
			i++;
			if (i < argc) {
				if (strcmp("Mo",argv[i]) == 0) {
					si->sp->M /= MU;
					i++;
					}
				else if (strcmp("MU",argv[i]) == 0) {
					i++;
					}
				else {
					si->sp->M /= MU;
					}
				}
			}
		else if (strcmp(argv[i],"-rs") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->rs = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-rcutoff") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->rcutoff = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-cvir") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->cvir = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-N0") == 0) {
			i++;
			if (i >= argc) usage();
			si->N0 = (INT) (atof(argv[i]));
			i++;
			}
		else if (strcmp(argv[i],"-soft0") == 0) {
			i++;
			if (i >= argc) usage();
			if (strcmp("LU",argv[i]) == 0) {
				si->soft0 = 1;
				}
			else if (strcmp("rs",argv[i]) == 0) {
				si->soft0 = 1;
				si->soft0_in_rs_units = 1;
				}
			else if (strcmp("rvir",argv[i]) == 0) {
				si->soft0 = 1;
				si->soft0_in_rvir_units = 1;
				}
			else if (strcmp("rcutoff",argv[i]) == 0) {
				si->soft0 = 1;
				si->soft0_in_rcutoff_units = 1;
				}
			else {
				si->soft0 = atof(argv[i]);
				}
			i++;
			if (i < argc) {
				if (strcmp("LU",argv[i]) == 0) {
					i++;
					}
				else if (strcmp("rs",argv[i]) == 0) {
					si->soft0_in_rs_units = 1;
					i++;
					}
				else if (strcmp("rvir",argv[i]) == 0) {
					si->soft0_in_rvir_units = 1;
					i++;
					}
				else if (strcmp("rcutoff",argv[i]) == 0) {
					si->soft0_in_rcutoff_units = 1;
					i++;
					}
				}
			}
		else if (strcmp(argv[i],"-rba_at_r0") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->rba_at_r0 = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-rba_r0") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->rba_r0 = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-rba_slope") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->rba_slope = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-rba_min") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->rba_min = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-rba_max") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->rba_max = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-rca_at_r0") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->rca_at_r0 = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-rca_r0") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->rca_r0 = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-rca_slope") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->rca_slope = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-rca_min") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->rca_min = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-rca_max") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->rca_max = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-alpha_at_r0") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->alpha_at_r0 = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-alpha_r0") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->alpha_r0 = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-alpha_slope") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->alpha_slope = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-alpha_min") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->alpha_min = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-alpha_max") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->alpha_max = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-beta_at_r0") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->beta_at_r0 = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-beta_r0") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->beta_r0 = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-beta_slope") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->beta_slope = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-beta_min") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->beta_min = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-beta_max") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->beta_max = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-gamma_at_r0") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->gamma_at_r0 = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-gamma_r0") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->gamma_r0 = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-gamma_slope") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->gamma_slope = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-gamma_min") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->gamma_min = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-gamma_max") == 0) {
			i++;
			if (i >= argc) usage();
			si->sp->gamma_max = atof(argv[i]);
			i++;
			}
		/*
		** Multi-mass parameters
		*/
		else if (strcmp(argv[i],"-rsi") == 0) {
			i++;
			if (i >= argc) usage();
			if (strcmp("LU",argv[i]) == 0) {
				si->rsi = 1;
				}
			else if (strcmp("rs",argv[i]) == 0) {
				si->rsi = 1;
				si->rsi_in_rs_units = 1;
				}
			else if (strcmp("rvir",argv[i]) == 0) {
				si->rsi = 1;
				si->rsi_in_rvir_units = 1;
				}
			else if (strcmp("rcutoff",argv[i]) == 0) {
				si->rsi = 1;
				si->rsi_in_rcutoff_units = 1;
				}
			else {
				si->rsi = atof(argv[i]);
				}
			i++;
			if (i < argc) {
				if (strcmp("LU",argv[i]) == 0) {
					i++;
					}
				else if (strcmp("rs",argv[i]) == 0) {
					si->rsi_in_rs_units = 1;
					i++;
					}
				else if (strcmp("rvir",argv[i]) == 0) {
					si->rsi_in_rvir_units = 1;
					i++;
					}
				else if (strcmp("rcutoff",argv[i]) == 0) {
					si->rsi_in_rcutoff_units = 1;
					i++;
					}
				}
			}
		else if (strcmp(argv[i],"-rso") == 0) {
			i++;
			if (i >= argc) usage();
			if (strcmp("LU",argv[i]) == 0) {
				si->rso = 1;
				}
			else if (strcmp("rs",argv[i]) == 0) {
				si->rso = 1;
				si->rso_in_rs_units = 1;
				}
			else if (strcmp("rvir",argv[i]) == 0) {
				si->rso = 1;
				si->rso_in_rvir_units = 1;
				}
			else if (strcmp("rcutoff",argv[i]) == 0) {
				si->rso = 1;
				si->rso_in_rcutoff_units = 1;
				}
			else {
				si->rso = atof(argv[i]);
				}
			i++;
			if (i < argc) {
				if (strcmp("LU",argv[i]) == 0) {
					i++;
					}
				else if (strcmp("rs",argv[i]) == 0) {
					si->rso_in_rs_units = 1;
					i++;
					}
				else if (strcmp("rvir",argv[i]) == 0) {
					si->rso_in_rvir_units = 1;
					i++;
					}
				else if (strcmp("rcutoff",argv[i]) == 0) {
					si->rso_in_rcutoff_units = 1;
					i++;
					}
				}
			}
		else if (strcmp(argv[i],"-Nshell") == 0) {
			i++;
			if (i >= argc) usage();
			si->Nshell = (INT) atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-Ismor") == 0) {
			i++;
			if (i >= argc) usage();
			si->Ismor = (INT) atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-DRMmax") == 0) {
			i++;
			if (i >= argc) usage();
			si->DRMmax = (INT) atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-rmor") == 0) {
			i++;
			if (i >= argc) usage();
			if (strcmp("LU",argv[i]) == 0) {
				si->rmor = 1;
				}
			else if (strcmp("rs",argv[i]) == 0) {
				si->rmor = 1;
				si->rmor_in_rs_units = 1;
				}
			else if (strcmp("rvir",argv[i]) == 0) {
				si->rmor = 1;
				si->rmor_in_rvir_units = 1;
				}
			else if (strcmp("rcutoff",argv[i]) == 0) {
				si->rmor = 1;
				si->rmor_in_rcutoff_units = 1;
				}
			else {
				si->rmor = atof(argv[i]);
				}
			i++;
			if (i < argc) {
				if (strcmp("LU",argv[i]) == 0) {
					i++;
					}
				else if (strcmp("rs",argv[i]) == 0) {
					si->rmor_in_rs_units = 1;
					i++;
					}
				else if (strcmp("rvir",argv[i]) == 0) {
					si->rmor_in_rvir_units = 1;
					i++;
					}
				else if (strcmp("rcutoff",argv[i]) == 0) {
					si->rmor_in_rcutoff_units = 1;
					i++;
					}
				}
			}
		/*
		** Special parameters
		*/
		else if (strcmp(argv[i],"-dfsf") == 0) {
			i++;
			if (i >= argc) usage();
			si->dfsf = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-Ngridr") == 0) {
			i++;
			if (i >= argc) usage();
			gi->Ngridr = atoi(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-Ngriddf") == 0) {
			i++;
			if (i >= argc) usage();
			gi->Ngriddf = atoi(argv[i]);
			i++;
			}
		/*
		** Black hole parameters
		*/
		else if (strcmp(argv[i],"-MBH") == 0) {
			gi->do_bh = 1;
			i++;
			if (i >= argc) usage();
			if (strcmp("Mo",argv[i]) == 0) {
				bh->mass = 1;
				}
			else if (strcmp("MU",argv[i]) == 0) {
				bh->mass = MU;
				}
			else {
				bh->mass = atof(argv[i]);
				}
			i++;
			if (i < argc) {
				if (strcmp("Mo",argv[i]) == 0) {
					bh->mass /= MU;
					i++;
					}
				else if (strcmp("MU",argv[i]) == 0) {
					i++;
					}
				else {
					bh->mass /= MU;
					}
				}
			}
		else if (strcmp(argv[i],"-softBH") == 0) {
			gi->do_bh = 1;
			i++;
			if (i >= argc) usage();
			bh->soft = atof(argv[i]);
			i++;
			}
		/*
		** Model name
		*/
		else if (strcmp(argv[i],"-name") == 0) {
			i++;
			if (i >= argc) usage();
			sprintf(gi->outputname,"%s",argv[i]);
			i++;
			}
		/*
		** Coordinates
		*/
		else if (strcmp(argv[i],"-spherical") == 0) {
			gi->coordinates = 0;
			i++;
			}
		else if (strcmp(argv[i],"-ellipsoidal") == 0) {
			gi->coordinates = 1;
			i++;
			}
		/*
		** Output parameters
		*/
		else if (strcmp(argv[i],"-ogr") == 0) {
			gi->output_gridr = 1;
			i++;
			}
		else if (strcmp(argv[i],"-ogdf") == 0) {
			gi->output_griddf = 1;
			i++;
			}
		else if (strcmp(argv[i],"-ots") == 0) {
			gi->output_tipsy_standard = 1;
			i++;
			}
		else if (strcmp(argv[i],"-otsdpp") == 0) {
			gi->output_tipsy_standard_dpp = 1;
			i++;
			}
		else if (strcmp(argv[i],"-po") == 0) {
			gi->positionsonly = 1;
			i++;
			}
		/*
		** Cosmological parameters
		*/
		else if (strcmp(argv[i],"-OmegaM0") == 0) {
			i++;
			if (i >= argc) usage();
			gi->OmegaM0 = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-OmegaK0") == 0) {
			i++;
			if (i >= argc) usage();
			gi->OmegaK0 = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-OmegaL0") == 0) {
			i++;
			if (i >= argc) usage();
			gi->OmegaL0 = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-h0") == 0) {
			i++;
			if (i >= argc) usage();
			gi->h0 = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-z") == 0) {
			i++;
			if (i >= argc) usage();
			gi->z = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-Deltavirz") == 0) {
			i++;
			if (i >= argc) usage();
			gi->Deltavirz = atof(argv[i]);
			i++;
			}
		/*
		** Special parameters
		*/
		else if (strcmp(argv[i],"-f_rinner") == 0) {
			i++;
			if (i >= argc) usage();
			gi->f_rinner = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-f_router") == 0) {
			i++;
			if (i >= argc) usage();
			gi->f_router = atof(argv[i]);
			i++;
			}
		else if (strcmp(argv[i],"-f_cutoff") == 0) {
			i++;
			if (i >= argc) usage();
			gi->f_cutoff = atof(argv[i]);
			i++;
			}
		/*
		** Failure
		*/
		else {
			usage();
			}
		}
	}
