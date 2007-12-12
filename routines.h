/* 
** routines.h
**
** Header file for routines.c
*/

void initialise_parameters(SI (*));
void initialise_black_hole(PARTICLE (*));
void check_main_parameters(const SI (*));
void calculate_parameters(const GI (*), SI (*)); 
void initialise_gridr(GI (*), PARTICLE (*), SI (*));
void calculate_virial_stuff(const GI (*), SI (*));
void set_remaining_parameters(const GI (*), SI (*));
void check_more_parameters(const GI (*), const SI (*));
void initialise_griddf(const GI (*), SI (*));
void initialise_shell(SI (*));
void set_positions(SI (*)); 
void set_velocities(const GI (*), SI (*));
void set_velocities_zero(SI (*));
void set_attributes(const GI (*), SI (*));
void refine(const GI (*), SI (*));
void searchroot(INT (*), DOUBLE (*));
void searchmin(INT (*), DOUBLE (*));
void double_particles(SI (*));
void correct_cm_velocity(SI (*));
void calculate_stuff(GI (*), PARTICLE (*), SI (*));
void transfer_particles(const PARTICLE (*), const SI (*), TIPSY_STRUCTURE (*));
void write_tipsy_standard_2(FILE (*), const PARTICLE (*), const SI (*));
void write_griddf(SI (*), FILE (*));
void write_gridr(GRIDR (*), FILE (*));
void usage(void);
