/*
** usage.c
*/

#include <stdio.h>
#include <stdlib.h>

/*
** Usage description
*/

void usage() {

    fprintf(stderr,"\n");
    fprintf(stderr,"You can specify the following arguments.\n\n");
    fprintf(stderr,"-a <value>          : alpha parameter in density profile\n");
    fprintf(stderr,"-b <value>          : beta parameter in density profile\n");
    fprintf(stderr,"-c <value>          : gamma parameter in density profile\n");
    fprintf(stderr,"-M <value> <unit>   : mass within rcutoff or total mass for finite mass models, <unit> = Mo or MU (optional - default: Mo)\n");
    fprintf(stderr,"-rs <value>         : scale radius [kpc]\n");
    fprintf(stderr,"-rcutoff <value>    : cutoff radius [kpc]\n");
    fprintf(stderr,"-cvir <value>       : virial concentration parameter\n");
    fprintf(stderr,"-N0 <value>         : number of particles within rsi\n");
    fprintf(stderr,"-soft0 <value>      : softening of particles within rsi [kpc]\n");
    fprintf(stderr,"-rsi <value>        : inner shell radius - rcutoff, rvir, rs or numerical value [kpc] (default: rcutoff / router)\n");
    fprintf(stderr,"-rso <value>        : outer shell radius - rcutoff, rvir, rs or numerical value [kpc] (default: rcutoff / router)\n");
    fprintf(stderr,"-Nshell <value>     : number of shells between rsi and rso (default: 0)\n");
    fprintf(stderr,"-DRMmax <value>     : maximum mass ratio between two neighbouring shells (default: 1)\n");
    fprintf(stderr,"-rmor <value>       : maximum orbital refinement radius - rcutoff, rvir, rs or numerical value [kpc] (default: 0 kpc)\n");
    fprintf(stderr,"-Ismor <value>      : maximum index of shell where orbital refinement is done (default: Nshell+1)\n");
    fprintf(stderr,"-dfsf <value>       : distribution function split factor (default: 0.01)\n");
    fprintf(stderr,"-MBH <value> <unit> : mass of black hole, <unit> = Mo or MU (optional - default: Mo)\n");
    fprintf(stderr,"-softBH <value>     : softening of black hole [kpc]\n");
    fprintf(stderr,"-name <value>       : name of the output file\n");
    fprintf(stderr,"-OmegaM0 <value>    : OmegaM0 (default: 0.3)\n");
    fprintf(stderr,"-OmegaK0 <value>    : OmegaK0 (default: 0.0)\n");
    fprintf(stderr,"-OmegaL0 <value>    : OmegaL0 (default: 0.7)\n");
    fprintf(stderr,"-h0 <value>         : Hubble parameter h0 (default: 0.7)\n");
    fprintf(stderr,"-z <value>          : Redshift z (default: 0.0)\n");
    fprintf(stderr,"-Deltavirz <value>  : virial overdensity (default: 178*[OmegaMz^0.45] if OmegaK0 = 0 and 178*[OmegaMz^0.30] if OmegaL0 = 0)\n");
    fprintf(stderr,"-ogr                : set this flag for outputting grid in r\n");
    fprintf(stderr,"-ogdf               : set this flag for outputting grid for distribution function\n");
    fprintf(stderr,"-ota                : set this flag for writing particles in tipsy ascii format\n");
    fprintf(stderr,"-otb                : set this flag for writing particles in tipsy binary format\n");
    fprintf(stderr,"-ots                : set this flag for writing particles in tipsy standard binary format\n");
    fprintf(stderr,"-ogb                : set this flag for writing particles in gadget binary format\n");
    fprintf(stderr,"-po                 : set this flag for initialising positions only (velocities are 0)\n");
    fprintf(stderr,"-randomseed <value> : set this flag for setting a value for a random seed (default: random value)\n");
    fprintf(stderr,"\n");
    exit(1);
    }
