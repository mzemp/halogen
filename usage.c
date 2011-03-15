/*
** usage.c
**
** Usage for HALOGEN
*/

#include <stdio.h>
#include <stdlib.h>

/*
** Usage description
*/

void usage() {

    fprintf(stderr,"\n");
    fprintf(stderr,"You can specify the following arguments:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"System specification:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-bulge                : all system parameters after this flag are set for the bulge system\n");
    fprintf(stderr,"-halo                 : all system parameters after this flag are set for the halo system\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"System parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-a <value>            : alpha parameter in density profile\n");
    fprintf(stderr,"-b <value>            : beta parameter in density profile\n");
    fprintf(stderr,"-c <value>            : gamma parameter in density profile\n");
    fprintf(stderr,"-M <value> <unit>     : mass within rcutoff or total mass for finite mass models, <unit>: Mo or MU (optional - default: Mo)\n");
    fprintf(stderr,"-rs <value>           : scale radius [LU]\n");
    fprintf(stderr,"-rcutoff <value>      : cutoff radius [LU]\n");
    fprintf(stderr,"-cvir <value>         : virial concentration parameter\n");
    fprintf(stderr,"-N0 <value>           : number of particles within rsi\n");
    fprintf(stderr,"-soft0 <value> <unit> : softening of particles within rsi, <unit>: LU, rs, rvir or rcutoff (optional - default: LU)\n");
    fprintf(stderr,"-rba <value>          : axis ratio b/a (default: 1)\n");
    fprintf(stderr,"-rca <value>          : axis ratio c/a (default: 1)\n");
    fprintf(stderr,"-rsi <value> <unit>   : inner shell radius (default: rcutoff / router), <unit>: LU, rs, rvir or rcutoff (optional - default: LU)\n");
    fprintf(stderr,"-rso <value> <unit>   : outer shell radius (default: rcutoff / router), <unit>: LU, rs, rvir or rcutoff (optional - default: LU)\n");
    fprintf(stderr,"-Nshell <value>       : number of shells between rsi and rso (default: 0)\n");
    fprintf(stderr,"-DRMmax <value>       : maximum mass ratio between two neighbouring shells (default: 1)\n");
    fprintf(stderr,"-rmor <value> <unit>  : maximum orbital refinement radius (default: 0 LU), <unit>: LU, rs, rvir or rcutoff (optional - default: LU)\n");
    fprintf(stderr,"-Ismor <value>        : maximum index of shell where orbital refinement is done (default: Nshell+1)\n");
    fprintf(stderr,"-dfsf <value>         : distribution function split factor (default: 0.01)\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Black hole parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-MBH <value> <unit>   : mass of black hole, <unit>: Mo or MU (optional - default: Mo)\n");
    fprintf(stderr,"-softBH <value>       : softening of black hole [LU]\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"General parameters:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-name <value>         : name of the output file\n");
    fprintf(stderr,"-Ngridr <value>       : number of grid points for grid in r (default: 2001)\n");
    fprintf(stderr,"-Ngriddf <value>      : number of grid points for grid for distribution function for each system (default: 101)\n");
    fprintf(stderr,"-OmegaM0 <value>      : OmegaM0 (default: 0.3)\n");
    fprintf(stderr,"-OmegaK0 <value>      : OmegaK0 (default: 0.0)\n");
    fprintf(stderr,"-OmegaL0 <value>      : OmegaL0 (default: 0.7)\n");
    fprintf(stderr,"-h0 <value>           : Hubble parameter h0 (default: 0.7)\n");
    fprintf(stderr,"-z <value>            : Redshift z (default: 0)\n");
    fprintf(stderr,"-Deltavirz <value>    : virial overdensity (default: 178*[OmegaMz^0.45] if OmegaK0 = 0 and 178*[OmegaMz^0.30] if OmegaL0 = 0)\n");
    fprintf(stderr,"-ogr                  : set this flag for outputting grid in r\n");
    fprintf(stderr,"-ogdf                 : set this flag for outputting grid for distribution function for each system\n");
    fprintf(stderr,"-ots                  : set this flag for writing particles in tipsy standard binary format\n");
    fprintf(stderr,"-otsdpp               : set this flag for writing particles in tipsy standard binary format with double precision positions\n");
    fprintf(stderr,"-spherical            : set this flag for initialising positions in spherical coordinates (default)\n");
    fprintf(stderr,"-ellipsoidal          : set this flag for initialising positions in ellipsoidal coordinates\n");
    fprintf(stderr,"-po                   : set this flag for initialising positions only (velocities are 0)\n");
    fprintf(stderr,"-randomseed <value>   : value for a random seed (default: random value)\n");
    fprintf(stderr,"\n");
    exit(1);
    }
