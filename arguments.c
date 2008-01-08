/*
** arguments.c
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "definitions.h"
#include "usage.h"

/*
** Routine for processing input arguments
*/

void process_arguments(int argc, char **argv, GI *gi, PARTICLE *bh, SI *halo) {

    INT i;

    i = 1;
    while (i < argc) {
	/*
	** Model parameters
	*/
	if (strcmp(argv[i],"-a") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->sp->alpha = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-b") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->sp->beta = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-c") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->sp->gamma = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-M") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->sp->M = atof(argv[i]);
	    i++;
	    if (strcmp("Mo",argv[i]) == 0) {
		halo->sp->M /= MU;
		i++;
		}
	    else if (strcmp("MU",argv[i]) == 0) {
		i++;
		}
	    else {
		halo->sp->M /= MU;
		}
	    }
	else if (strcmp(argv[i],"-rs") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->sp->rs = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-rcutoff") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->sp->rcutoff = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-cvir") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->sp->cvir = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-N0") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->N0 = (INT) (atof(argv[i]));
	    i++;
	    }
	else if (strcmp(argv[i],"-soft0") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->soft0 = atof(argv[i]);
	    i++;
	    }
	/*
	** Multi-mass parameters
	*/
	else if (strcmp(argv[i],"-rsi") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    if (strcmp("rs",argv[i]) == 0) {
		halo->set_rsi_to_rs = 1;
		}
	    else if (strcmp("rvir",argv[i]) == 0) {
		halo->set_rsi_to_rvir = 1;
		}
	    else if (strcmp("rcutoff",argv[i]) == 0) {
		halo->set_rsi_to_rcutoff = 1;
		}
	    else {
		halo->rsi = atof(argv[i]);
		}
	    i++;
	    }
	else if (strcmp(argv[i],"-rso") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    if (strcmp("rs",argv[i]) == 0) {
		halo->set_rso_to_rs = 1;
		}
	    else if (strcmp("rvir",argv[i]) == 0) {
		halo->set_rso_to_rvir = 1;
		}
	    else if (strcmp("rcutoff",argv[i]) == 0) {
		halo->set_rso_to_rcutoff = 1;
		}
	    else {
		halo->rso = atof(argv[i]);
		}
	    i++;
	    }
	else if (strcmp(argv[i],"-Nshell") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->Nshell = (INT) atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-Ismor") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->Ismor = (INT) atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-DRMmax") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->DRMmax = (INT) atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-rmor") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    if (strcmp("rs",argv[i]) == 0) {
		halo->set_rmor_to_rs = 1;
		}
	    else if (strcmp("rvir",argv[i]) == 0) {
		halo->set_rmor_to_rvir = 1;
		}
	    else if (strcmp("rcutoff",argv[i]) == 0) {
		halo->set_rmor_to_rcutoff = 1;
		}
	    else {
		halo->rmor = atof(argv[i]);
		}
	    i++;
	    }
	/*
	** Special parameters
	*/
	else if (strcmp(argv[i],"-dfsf") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    halo->dfsf = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-Ngridr") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    gi->Ngridr = atoi(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-Ngriddf") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    gi->Ngriddf = atoi(argv[i]);
	    i++;
	    }
	/*
	** Black hole parameters
	*/
	else if (strcmp(argv[i],"-MBH") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    bh->mass = atof(argv[i]);
	    i++;
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
	else if (strcmp(argv[i],"-softBH") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    bh->soft = atof(argv[i]);
	    i++;
	    }
	/*
	** Model name
	*/
	else if (strcmp(argv[i],"-name") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    sprintf(gi->inputname,"%s",argv[i]);
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
	else if (strcmp(argv[i],"-ogb") == 0) {
	    gi->output_gadget_binary = 1;
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
	    if (i >= argc) {
		usage();
		}
	    gi->OmegaM0 = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-OmegaK0") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    gi->OmegaK0 = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-OmegaL0") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    gi->OmegaL0 = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-h0") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    gi->h0 = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-z") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    gi->z = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-Deltavirz") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    gi->Deltavirz = atof(argv[i]);
	    i++;
	    }
	/*
	** Special parameters
	*/
	else if (strcmp(argv[i],"-randomseed") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    gi->randomseed = atof(argv[i]);
	    i++;
	    }
	/*
	** Help or failure
	*/
	else if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
	    usage();
	    }
	else {
	    usage();
	    }
	}
    }
