/*
** initialise.h
*/

void initialise_general_info(GI (*));
void initialise_system(SI (*));
void initialise_particle(PARTICLE (*));
void initialise_gridr(GI (*), PARTICLE (*), SI (*), SI (*));
void initialise_griddf(const GI (*), SI (*));
void initialise_shell(const GI (*), SI (*));
