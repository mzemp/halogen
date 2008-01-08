/* 
** halogen.c
**
** Program written in order to generate multi-mass spherical structures
**
** written by Marcel Zemp
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
#include "definitions.h"
#include "arguments.h"
#include "functions.h"
#include "routines.h"
#include "write.h"
#include "check.h"
#include "usage.h"

int main(int argc, char **argv) {

    /*
    ** Variables
    */

    GI *gi;
    PARTICLE *bh;
    SI *halo;
    SI *bulge;
    CHAR FILENAME[STRINGSIZE];
    FILE *file;

    /*
    ** Initialise structures for reading parameters and start clock
    */

    bh = malloc(sizeof(PARTICLE));
    assert(bh != NULL);

    halo = malloc(sizeof(SI));
    assert(halo != NULL);
    halo->sp = malloc(sizeof(SP));
    assert(halo->sp != NULL);

    bulge = malloc(sizeof(SI));
    assert(bulge != NULL);
    bulge->sp = malloc(sizeof(SP));
    assert(bulge->sp != NULL);

    gi = malloc(sizeof(GI));
    assert(gi != NULL);

    gi->t[0] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);

    /*
    ** Set standard values for parameters
    */

    sprintf(gi->inputname,"none");
    sprintf(halo->systemname,"halo");
    sprintf(bulge->systemname,"bulge");

    initialise_general_info(gi);
    initialise_particle(bh);
    initialise_system(halo);

    /*
    ** Read in and process arguments
    */

    process_arguments(argc,argv,gi,bh,halo);

    fprintf(stderr,"Checking parameters, calculating halo properties and initialising grid in r... \n");

    /*
    ** Initialise random number generator
    */

    srand(gi->randomseed);

    /*
    ** Check main input parameters
    */

    check_main_parameters_general_info(gi);
    check_main_parameters_system(halo);

    /*
    ** Allocate memory for the structures
    */

    allocate_general_info(gi);
    allocate_system(gi,halo);

    /*
    ** Derived cosmological parameters
    */

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
    
    if (gi->output_gridr == 1) {
	sprintf(FILENAME,"%s.gridr.dat",gi->inputname);
	file = fopen(FILENAME,"w");
	assert(file != NULL);
	write_gridr(file,gi);
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

    check_more_parameters_system(gi,halo);

    /*
    ** Initialise griddf
    */

    gi->t[1] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds.\nInitialising grid for distribution function... \n",gi->t[1]-gi->t[0]);

    if (gi->positionsonly == 0) {
	initialise_griddf(gi,halo);
	if (gi->output_griddf == 1) {
	    sprintf(FILENAME,"%s.griddf.halo.dat",gi->inputname);
	    file = fopen(FILENAME,"w");
	    assert(file != NULL);
	    write_griddf(file,gi,halo);
	    fclose(file);
	    }
	}

    /*
    ** Initialise shell
    */

    gi->t[2] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds\nInitialising shells... \n",gi->t[2]-gi->t[1]);

    initialise_shell(gi,halo);

    /*
    ** Set particle positions
    */

    gi->t[3] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds.\nSetting particle positions... \n",gi->t[3]-gi->t[2]);

    set_positions(gi,halo);

    /*
    ** Set particle velocities
    */

    gi->t[4] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds.\nSetting particle velocities... \n",gi->t[4]-gi->t[3]);

    if (gi->positionsonly == 0) {
	set_velocities(gi,halo);
	}
    else {
	set_velocities_zero(halo);
	}

    /*
    ** Set remaining attributes
    */

    gi->t[5] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds.\nSetting remaining particle attributes... \n",gi->t[5]-gi->t[4]);

    set_attributes(gi,halo);

    /*
    ** Do orbit dependent refining
    */

    gi->t[6] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds.\nDoing orbit dependent refining... \n",gi->t[6]-gi->t[5]);

    refine(gi,halo);

    /*
    ** Calculate a few things and do center of mass correction
    */

    gi->t[7] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds\nCalculating a few things and correct center of mass position and velocity... \n",gi->t[7]-gi->t[6]);

    double_particles(halo);
    calculate_stuff(gi,bh,halo);

    /*
    ** Write Output
    */

    gi->t[8] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
    fprintf(stderr,"Done in "OFD1" seconds\nWriting output... \n",gi->t[8]-gi->t[7]);

    if (gi->output_tipsy_standard == 1) {
	sprintf(FILENAME,"%s.tipsy.std",gi->inputname);
	file = fopen(FILENAME,"w");
	assert(file != NULL);
	write_tipsy_standard_2(file,bh,halo);
	fclose(file);
	}

    /* 
    ** Print some output in file
    */

    sprintf(FILENAME,"%s.out",gi->inputname);
    file = fopen(FILENAME,"w");
    assert(file != NULL);
    write_general_output(file,argc,argv,gi,bh,halo);
    fclose(file);

    fprintf(stderr,"Done in "OFD1" seconds\nTotal time needed was "OFD1" seconds\n",gi->t[9]-gi->t[8],gi->t[9]-gi->t[0]);

    free(bh);
    free(halo);
    free(gi);

    exit(0);
    } /* end of main function */


