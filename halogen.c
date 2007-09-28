/* halogen.c
**
** Program written in order to generate multi-mass spherical structures
**
** written by Marcel Zemp (mzemp@ucolick.org)
**
** This program works in units where
** 
** [G] = 1 [L]^3 [T]^-2 [M]^-1
** [L] = kpc
** [T] = Gyr 
** [V] = kpc Gyr^-1
** [M] = MU
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <time.h>
#include "IOfunctions.h"
#include "definitions.h"
#include "functions.h"
#include "routines.h"

int main(int argc, char **argv) {

    /*
    ** Variables
    */

    INT i, k;
    INT output_gridr, output_griddf;
    INT output_tipsy_ascii, output_tipsy_binary, output_tipsy_standard;
    INT output_gadget_binary;
    INT positionsonly;
    DOUBLE randomseed, OmegaMz;
    DOUBLE t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;
    PARTICLE *bh;
    SI *halo;
    GI *gi;
    TIPSY_STRUCTURE *ts;
    CHAR FILENAME[STRINGSIZE], INPUTNAME[STRINGSIZE];
    FILE *file;

    randomseed = time(NULL);

    t0 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);

    /*
    ** Initialise fixed structures
    */

    bh = malloc(sizeof(PARTICLE));
    assert(bh != NULL);

    halo = malloc(sizeof(SI));
    assert(halo != NULL);
    halo->sp = malloc(sizeof(SP));
    assert(halo->sp != NULL);
    halo->griddf = malloc(sizeof(GRIDDF));
    assert(halo->griddf != NULL);

    gi = malloc(sizeof(GI));
    assert(gi != NULL);
    gi->stuff = malloc(sizeof(STUFF));
    assert(gi->stuff != NULL);
    gi->gridr = malloc(sizeof(GRIDR));
    assert(gi->gridr != NULL);

    ts = malloc(sizeof(TIPSY_STRUCTURE));
    assert(ts != NULL);

    /*
    ** Set standard values for parameters
    */

    gi->OmegaM0 = 0.3;
    gi->OmegaK0 = 0;
    gi->OmegaL0 = 0.7;
    gi->h0 = 0.7;
    gi->z = 0;
    gi->rhocritz = -1;
    gi->Deltavirz = -1;

    output_gridr = 0;
    output_griddf = 0;
    output_tipsy_ascii = 0;
    output_tipsy_binary = 0;
    output_tipsy_standard = 0;
    output_gadget_binary = 0;
    positionsonly = 0;

    sprintf(INPUTNAME,"none");
    sprintf(halo->systemname,"halo");

    initialise_black_hole(bh);
    initialise_parameters(halo);

    /*
    ** Read in and calculate model parameters
    */

    i = 1;
    while (i < argc) {
	/*
	** Halo parameters
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
	    sprintf(INPUTNAME,"%s",argv[i]);
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
	** Output parameters
	*/
	else if (strcmp(argv[i],"-ogr") == 0) {
	    output_gridr = 1;
	    i++;
	    }
	else if (strcmp(argv[i],"-ogdf") == 0) {
	    output_griddf = 1;
	    i++;
	    }
	else if (strcmp(argv[i],"-ota") == 0) {
	    output_tipsy_ascii = 1;
	    i++;
	    }
	else if (strcmp(argv[i],"-otb") == 0) {
	    output_tipsy_binary = 1;
	    i++;
	    }
	else if (strcmp(argv[i],"-ots") == 0) {
	    output_tipsy_standard = 1;
	    i++;
	    }
	else if (strcmp(argv[i],"-ogb") == 0) {
	    output_gadget_binary = 1;
	    i++;
	    }
	else if (strcmp(argv[i],"-po") == 0) {
	    positionsonly = 1;
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
	    randomseed = atof(argv[i]);
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

    fprintf(stderr,"Checking parameters, calculating halo properties and initialising grid in r... \n");
    srand(randomseed);

    /*
    ** Check main input parameters
    */

    if (strcmp(INPUTNAME,"none") == 0) {
	fprintf(stderr,"You have not set a name for the output model.\n");
	usage();
	}
    if ((NGRIDR-1) % (NGRIDDF-1) != 0) {
	fprintf(stderr,"Bad choice of NGRIDR and NGRIDDF!\n");
	fprintf(stderr,"These numbers have to fulfill the condition (NGRIDR-1) mod (NGRIDDF-1) == 0.\n");
	usage();
	}

    check_main_parameters(halo);

    /*
    ** Derived cosmological parameters
    */

    gi->rhocritz = 3/(8*M_PI*G)*(pow((100*gi->h0*Ecosmo(gi)*VelConvertFac/1000),2));
    if (gi->Deltavirz == -1) {
	OmegaMz = gi->OmegaM0*pow((1+gi->z),3)/pow(Ecosmo(gi),2);
	if (gi->OmegaK0 == 0) {
	    gi->Deltavirz = 178*(pow(OmegaMz,0.45));
	    }
	else if (gi->OmegaL0 == 0) {
	    gi->Deltavirz = 178*(pow(OmegaMz,0.3));
	    }
	else {
	    fprintf(stderr,"HALOGEN can't calculate a value for Deltavirz for that choice of cosmological parameters.\n");
	    fprintf(stderr,"Set a value vor Deltavirz by hand or choose a different cosmology.\n");
	    usage();
	    }
	}

    calculate_parameters(gi,halo);

    /*
    ** Initialise gridr
    */

    gi->rinner = FACTORRINNER*halo->sp->rs;
    gi->router = halo->sp->rs;
    while (rho(halo->sp->rs,halo)/rho(gi->router,halo) < FACTORROUTER) {
	gi->router = gi->router*10;
	}

    initialise_gridr(gi,bh,halo);
    
    if (output_gridr == 1) {
	sprintf(FILENAME,"%s.gridr.dat",INPUTNAME);
	file = fopen(FILENAME,"w");
	assert(file != NULL);
	write_gridr(gi->gridr,file);
	fclose(file);
	}

    /*
    ** Calculate virial stuff for finite mass models and cutoff models
    */

    calculate_virial_stuff(gi,halo);

    /*
    ** Set remaining parameters
    */

    set_remaining_parameters(gi,halo);

    /*
    ** Check some more things
    */

    check_more_parameters(gi,halo);

    /*
    ** Initialise griddf
    */

    t1 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds.\nInitialising grid for distribution function... \n",t1-t0);

    if (positionsonly == 0) {
	initialise_griddf(gi,halo);
	if (output_griddf == 1) {
	    sprintf(FILENAME,"%s.griddf.halo.dat",INPUTNAME);
	    file = fopen(FILENAME,"w");
	    assert(file != NULL);
	    write_griddf(halo,file);
	    fclose(file);
	    }
	}

    /*
    ** Initialise shell
    */

    t2 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds\nInitialising shells... \n",t2-t1);

    initialise_shell(halo);

    /*
    ** Set particle positions
    */

    t3 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds.\nSetting particle positions... \n",t3-t2);

    set_positions(halo);

    /*
    ** Set particle velocities
    */

    t4 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds.\nSetting particle velocities... \n",t4-t3);

    if (positionsonly == 0) {
	set_velocities(gi,halo);
	}
    else {
	set_velocities_zero(halo);
	}

    /*
    ** Set remaining attributes
    */

    t5 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds.\nSetting remaining particle attributes... \n",t5-t4);

    set_attributes(gi,halo);

    /*
    ** Do orbit dependent refining
    */

    t6 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds.\nDoing orbit dependent refining... \n",t6-t5);

    refine(gi,halo);

    /*
    ** Calculate a few things and do center of mass correction
    */

    t7 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds\nCalculating a few things and correct center of mass position and velocity... \n",t7-t6);

    double_particles(halo);
    calculate_stuff(gi,bh,halo);

    /*
    ** Write Output
    */

    t8 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds\nWriting Output... \n",t8-t7);

    transfer_particles(bh,halo,ts);
    
    if (output_tipsy_ascii == 1) {
	sprintf(FILENAME,"%s.tipsy.ascii",INPUTNAME);
	file = fopen(FILENAME,"w");
	assert(file != NULL);
	write_tipsy_ascii(file,ts);
	fclose(file);
	}
    if (output_tipsy_binary == 1) {
	sprintf(FILENAME,"%s.tipsy.bin",INPUTNAME);
	file = fopen(FILENAME,"w");
	assert(file != NULL);
	write_tipsy_binary(file,ts);
	fclose(file);
	}
    if (output_tipsy_standard == 1) {
	sprintf(FILENAME,"%s.tipsy.std",INPUTNAME);
	file = fopen(FILENAME,"w");
	assert(file != NULL);
	write_tipsy_standard(file,ts);
	fclose(file);
	}
    if (output_gadget_binary == 1) {
	sprintf(FILENAME,"%s.gadget.bin",INPUTNAME);
	file = fopen(FILENAME,"w");
	assert(file != NULL);
	write_gadget_binary(file,ts,1,0,0,0,3,1,1000.0/VelConvertFac);
	fclose(file);
	}

    /* 
    ** Print some output in file
    */

    sprintf(FILENAME,"%s.out",INPUTNAME);
    file = fopen(FILENAME,"w");
    assert(file != NULL);
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
	fprintf(file,"Neff within rcutoff = "OFD1"\n",Menc(halo->sp->rcutoff,gi)/halo->shell[0].mass);
	}
    fprintf(file,"Ninitialtot         = "OFD3" = "OFI1"\n",(DOUBLE)gi->stuff->Ninitialtot,gi->stuff->Ninitialtot);
    fprintf(file,"Nnosplittot         = "OFD3" = "OFI1"\n",(DOUBLE)gi->stuff->Nnosplittot,gi->stuff->Nnosplittot);
    fprintf(file,"Nnewtot             = "OFD3" = "OFI1"\n",(DOUBLE)gi->stuff->Nnewtot,gi->stuff->Nnewtot);
    fprintf(file,"rimp                = "OFD3" kpc\n",halo->rimp);
    fprintf(file,"r1                  = "OFD3" kpc\n",halo->r1);
    fprintf(file,"r100                = "OFD3" kpc\n",halo->r100);
    fprintf(file,"Tdyn(rimp)          = "OFD3" Gyr\n",Tdyn(halo->rimp,gi));
    fprintf(file,"Tdyn(r1)            = "OFD3" Gyr\n",Tdyn(halo->r1,gi));
    fprintf(file,"Tdyn(r100)          = "OFD3" Gyr\n",Tdyn(halo->r100,gi));
    fprintf(file,"Mtheo               = "OFD3" Mo = "OFD3" MU\n",Menc(gi->router,gi)*MU,Menc(gi->router,gi));
    fprintf(file,"Msamp               = "OFD3" Mo = "OFD3" MU\n",gi->stuff->Mp*MU,gi->stuff->Mp);
    fprintf(file,"(Msamp-Mtheo)/Mtheo = "OFD3"\n",gi->stuff->Mp/Menc(gi->router,gi)-1.0);
    fprintf(file,"Random seed         = "OFD3"\n",randomseed);
    fprintf(file,"\n");
    fprintf(file,"Times for individual steps\n\n");
    fprintf(file,"Calculation of halo properties and initialisation of grid in r: "OFD1" seconds.\n",t1-t0);
    fprintf(file,"Initialisation of grid for distribution function: "OFD1" seconds.\n",t2-t1);
    fprintf(file,"Initialisation of shells: "OFD1" seconds.\n",t3-t2);
    fprintf(file,"Setting particle positions: "OFD1" seconds\n",t4-t3);
    fprintf(file,"Setting particle velocities: "OFD1" seconds\n",t5-t4);
    fprintf(file,"Setting remaining particle attributes: "OFD1" seconds\n",t6-t5);
    fprintf(file,"Orbit refinement: "OFD1" seconds\n",t7-t6);
    fprintf(file,"Calculating a few things and correct center of mass: "OFD1" seconds\n",t8-t7);
    t9 = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(file,"Writing output: "OFD1" seconds\n",t9-t8);
    fprintf(file,"Total time: "OFD1" seconds\n",t9-t0);
    fclose(file);
   
    fprintf(stderr,"Done in "OFD1" seconds\nTotal time needed was "OFD1" seconds\n",t9-t8,t9-t0);

    free(bh);
    free(halo);
    free(gi);
    exit(0);

    } /* end of main function */


