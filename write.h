/*
** write.h
*/

void write_tipsy_standard_2(FILE (*), const PARTICLE (*), const SI (*));
void write_griddf(FILE (*), const GI (*), const SI (*));
void write_gridr(FILE (*), const GI (*));
void write_general_output(FILE (*), int, char (**), GI (*), const PARTICLE (*), const SI(*));
