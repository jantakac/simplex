#ifndef SIMPLEX_TABLE_H
#define SIMPLEX_TABLE_H

#include <stdio.h>
#include <locale.h>
#include <stdint.h>
#include <stdlib.h>
#include <float.h>

typedef struct SimplexTable SimplexTable;
SimplexTable *simplex_table_create(void);
void simplex_table_solve(SimplexTable *self);
float simplex_table_elem_val(SimplexTable *self, size_t i, size_t j);
float *simplex_table_elem_ptr(SimplexTable *self, size_t row, size_t column);
size_t simplex_table_size(SimplexTable *self);
float simplex_table_obj_func(SimplexTable *self);
void simplex_table_print(SimplexTable *self);
void simplex_table_print_solution(SimplexTable *self);

void simplex_table_destroy(SimplexTable *self);

#endif