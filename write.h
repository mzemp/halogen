/*
** write.h
*/

void write_gridr_total(FILE (*), const GI (*));
void write_gridr_system(FILE (*), const GI (*), const SI (*));
void write_griddf_system(FILE (*), const GI (*), const SI (*));
void write_tipsy_standard_halogen(FILE (*), const GI (*), const PARTICLE (*), const SI (*), const SI (*));
void write_tipsy_standard_dpp_halogen(FILE (*), const GI (*), const PARTICLE (*), const SI (*), const SI (*));
void write_general_output(FILE (*), int, char (**), GI (*), const PARTICLE (*), const SI (*), const SI (*));
void write_output_system(FILE (*), const GI (*), const SI (*));
