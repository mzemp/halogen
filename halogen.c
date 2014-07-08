/* 
** halogen.c
**
** Program written in order to generate multi-mass spherical structures
**
** Written by Marcel Zemp
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
#include "arguments.h"
#include "allocate.h"
#include "check.h"
#include "initialise.h"
#include "usage.h"
#include "write.h"

int main(int argc, char **argv) {

	/*
	** Variables
	*/

	GI *gi;
	PARTICLE *bh;
	SI *bulge;
	SI *halo;
	CHAR FILENAME[STRINGSIZE];
	FILE *file;

	/*
	** Initialise structures for reading parameters and start clock
	*/

	gi = malloc(sizeof(GI));
	assert(gi != NULL);
	bh = malloc(sizeof(PARTICLE));
	assert(bh != NULL);
	bulge = malloc(sizeof(SI));
	assert(bulge != NULL);
	halo = malloc(sizeof(SI));
	assert(halo != NULL);

	gi->t[0] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);

	/*
	** Set standard values for parameters
	*/

	sprintf(gi->outputname,"none");
	sprintf(bulge->systemname,"bulge");
	sprintf(halo->systemname,"halo");

	initialise_general(gi);
	initialise_particle(bh);
	initialise_system(bulge);
	initialise_system(halo);

	/*
	** Read in and process arguments
	*/

	process_arguments(argc,argv,gi,bh,bulge,halo);

	fprintf(stderr,"Checking parameters, calculating halo properties and initialising grid in r... \n");

	/*
	** Initialise random number generator
	*/

	srand(gi->randomseed);

	/*
	** Check main input parameters
	*/

	check_main_parameters_general(gi);
	if (gi->do_bulge == 1) {
		check_main_parameters_system(bulge);
		}
	if (gi->do_halo == 1) {
		check_main_parameters_system(halo);
		}

	/*
	** Allocate memory for the structures
	*/

	allocate_general(gi);
	if (gi->do_bulge == 1) {
		allocate_system(gi,bulge);
		}
	if (gi->do_halo == 1) {
		allocate_system(gi,halo);
		}

	/*
	** Calculate parameters
	*/

	calculate_parameters_general(gi);
	if (gi->do_bulge == 1) {
		calculate_parameters_system(gi,bulge);
		}
	if (gi->do_halo == 1) {
		calculate_parameters_system(gi,halo);
		}

	/*
	** Initialise gridr
	*/

	initialise_gridr(gi,bh,bulge,halo);

	if (gi->output_gridr == 1) {
		sprintf(FILENAME,"%s.gridr.total.dat",gi->outputname);
		file = fopen(FILENAME,"w");
		assert(file != NULL);
		write_gridr_total(file,gi);
		fclose(file);
		if (gi->do_bulge == 1) {
			sprintf(FILENAME,"%s.gridr.bulge.dat",gi->outputname);
			file = fopen(FILENAME,"w");
			write_gridr_system(file,gi,bulge);
			fclose(file);
			}
		if (gi->do_halo == 1) {
			sprintf(FILENAME,"%s.gridr.halo.dat",gi->outputname);
			file = fopen(FILENAME,"w");
			write_gridr_system(file,gi,halo);
			fclose(file);
			}
		}

	/*
	** Calculate virial stuff for finite mass models and cutoff models
	*/

	if (gi->do_bulge == 1) {
		calculate_virial_stuff(gi,bulge);
		}
	if (gi->do_halo == 1) {
		calculate_virial_stuff(gi,halo);
		}

	/*
	** Set remaining parameters
	*/

	if (gi->do_bulge == 1) {
		set_remaining_parameters(gi,bulge);
		}
	if (gi->do_halo == 1) {
		set_remaining_parameters(gi,halo);
		}

	/*
	** Check some more things
	*/

	if (gi->do_bulge == 1) {
		check_more_parameters_system(gi,bulge);
		}
	if (gi->do_halo == 1) {
		check_more_parameters_system(gi,halo);
		}

	/*
	** Initialise griddf
	*/

	gi->t[1] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
	fprintf(stderr,"Done in "OFD1" seconds.\nInitialising grid for distribution function... \n",gi->t[1]-gi->t[0]);

	if (gi->positionsonly == 0) {
		if (gi->do_bulge == 1) {
			initialise_griddf(gi,bulge);
			if (gi->output_griddf == 1) {
				sprintf(FILENAME,"%s.griddf.bulge.dat",gi->outputname);
				file = fopen(FILENAME,"w");
				assert(file != NULL);
				write_griddf_system(file,gi,bulge);
				fclose(file);
				}
			}
		if (gi->do_halo == 1) {
			initialise_griddf(gi,halo);
			if (gi->output_griddf == 1) {
				sprintf(FILENAME,"%s.griddf.halo.dat",gi->outputname);
				file = fopen(FILENAME,"w");
				assert(file != NULL);
				write_griddf_system(file,gi,halo);
				fclose(file);
				}
			}
		}

	/*
	** Initialise shell
	*/

	gi->t[2] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
	fprintf(stderr,"Done in "OFD1" seconds\nInitialising shells... \n",gi->t[2]-gi->t[1]);

	if (gi->do_bulge == 1) {
		initialise_shell(gi,bulge);
		}
	if (gi->do_halo == 1) {
		initialise_shell(gi,halo);
		}

	/*
	** Set particle positions
	*/

	gi->t[3] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
	fprintf(stderr,"Done in "OFD1" seconds.\nSetting particle positions... \n",gi->t[3]-gi->t[2]);

	if (gi->do_bulge == 1) {
		set_positions(gi,bulge);
		}
	if (gi->do_halo == 1) {
		set_positions(gi,halo);
		}

	/*
	** Set particle velocities
	*/

	gi->t[4] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
	fprintf(stderr,"Done in "OFD1" seconds.\nSetting particle velocities... \n",gi->t[4]-gi->t[3]);

	if (gi->positionsonly == 0) {
		if (gi->do_bulge == 1) {
			set_velocities(gi,bulge);
			}
		if (gi->do_halo == 1) {
			set_velocities(gi,halo);
			}
		}
	else {
		if (gi->do_bulge == 1) {
			set_velocities_zero(bulge);
			}
		if (gi->do_halo == 1) {
			set_velocities_zero(halo);
			}
		}

	/*
	** Set remaining attributes
	*/

	gi->t[5] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
	fprintf(stderr,"Done in "OFD1" seconds.\nSetting remaining particle attributes... \n",gi->t[5]-gi->t[4]);

	if (gi->do_bulge == 1) {
		set_attributes(gi,bulge);
		}
	if (gi->do_halo == 1) {
		set_attributes(gi,halo);
		}

	/*
	** Do orbit dependent refining
	*/

	gi->t[6] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
	fprintf(stderr,"Done in "OFD1" seconds.\nDoing orbit dependent refining... \n",gi->t[6]-gi->t[5]);

	if (gi->do_bulge == 1) {
		refine(gi,bulge);
		}
	if (gi->do_halo == 1) {
		refine(gi,halo);
		}

	/*
	** Calculate a few things and do center of mass correction
	*/

	gi->t[7] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
	fprintf(stderr,"Done in "OFD1" seconds\nCalculating a few things and correct center of mass position and velocity... \n",gi->t[7]-gi->t[6]);

	if (gi->do_bulge == 1) {
		double_particles(bulge);
		calculate_samplinginfo_system(gi,bulge);
		}
	if (gi->do_halo == 1) {
		double_particles(halo);
		calculate_samplinginfo_system(gi,halo);
		}
	calculate_samplinginfo_general(gi,bh);

	/*
	** Write Output
	*/

	gi->t[8] = ((DOUBLE) clock())/((DOUBLE) CLOCKS_PER_SEC);
	fprintf(stderr,"Done in "OFD1" seconds\nWriting output... \n",gi->t[8]-gi->t[7]);

	if (gi->output_tipsy_standard == 1) {
		sprintf(FILENAME,"%s.tipsy.std",gi->outputname);
		file = fopen(FILENAME,"w");
		assert(file != NULL);
		write_tipsy_xdr_halogen(file,gi,bh,bulge,halo);
		fclose(file);
		}
	if (gi->output_tipsy_standard_dpp == 1) {
		sprintf(FILENAME,"%s.tipsy.dpp.std",gi->outputname);
		file = fopen(FILENAME,"w");
		assert(file != NULL);
		write_tipsy_xdr_dpp_halogen(file,gi,bh,bulge,halo);
		fclose(file);
		}

	/* 
	** Print some output in file
	*/

	sprintf(FILENAME,"%s.info.dat",gi->outputname);
	file = fopen(FILENAME,"w");
	assert(file != NULL);
	write_general_output(file,argc,argv,gi,bh,bulge,halo);
	fclose(file);

	fprintf(stderr,"Done in "OFD1" seconds\nTotal time needed was "OFD1" seconds\n",gi->t[9]-gi->t[8],gi->t[9]-gi->t[0]);

	free(gi);
	free(bh);
	free(halo);
	free(bulge);

	exit(0);
	} /* end of main function */
