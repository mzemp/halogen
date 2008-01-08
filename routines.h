/* 
** routines.h
*/

void calculate_parameters(const GI (*), SI (*)); 
void calculate_virial_stuff(const GI (*), SI (*));
void set_remaining_parameters(const GI (*), SI (*));
void set_positions(const GI (*), SI (*)); 
void set_velocities(const GI (*), SI (*));
void set_velocities_zero(SI (*));
void set_attributes(const GI (*), SI (*));
void refine(const GI (*), SI (*));
void searchroot(const GI (*), INT (*), DOUBLE (*));
void searchmin(const GI (*), INT (*), DOUBLE (*));
void double_particles(SI (*));
void correct_cm_velocity(SI (*));
void calculate_stuff(GI (*), PARTICLE (*), SI (*));
