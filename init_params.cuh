/**
*\fn void read_param(void)
*\brief function used to read the parameter file (named params.in by default)
*\param void
*\return void
*/

#include <stdio.h>
#include <stdlib.h>

#ifndef CUDAUBUNTU_INIT_PARAMS_CUH
#define CUDAUBUNTU_INIT_PARAMS_CUH
#include "def_types.h"
void read_table_mat(material_data *tab_mat_part);
void read_geom_param(geom_struct *geom_para);
int get_line_char(char *detect,char **tab_char);
void adi_params(material_data *prop_mat_part,geom_struct *geom);
#endif //CUDAUBUNTU_INIT_PARAMS_CUH
