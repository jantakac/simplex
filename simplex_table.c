#include "simplex_table.h"

typedef struct SimplexTable
{
    size_t rows;
    size_t cols;
    size_t *basic_vars_idxs;
    float *flat_table;
} SimplexTable;

static bool simplex_table_is_var_in_basis(SimplexTable *self, size_t var_idx);

SimplexTable *simplex_table_create(void)
{
    SimplexTable *self = malloc(sizeof(SimplexTable));

    FILE *fhandler = fopen("priklad1.txt", "r");
    if (!fhandler)
        return nullptr;

    setlocale(LC_NUMERIC, "sk_SK.utf8");

    fscanf(fhandler, "%zd %zd", &(self->rows), &(self->cols));
    ++self->rows;
    self->cols += 2;
    self->flat_table = malloc(sizeof(float) * self->rows * self->cols);
    self->basic_vars_idxs = malloc(sizeof(size_t) * (self->rows - 1));

    // loading the flat simplex table
    self->flat_table[0] = 1;
    for (size_t i = 1; i < self->rows * self->cols; ++i)
    {
        if (i % self->cols == 0)
            self->flat_table[i] = 0;
        else
        {
            fscanf(fhandler, "%f", &(self->flat_table[i]));
            if (i < self->cols)
                self->flat_table[i] = -self->flat_table[i];
        }
    }

    for (size_t i = 0; i < self->rows - 1; ++i)
        fscanf(fhandler, "%zd", &(self->basic_vars_idxs[i]));

    fclose(fhandler);
    return self;
}

void simplex_table_solve(SimplexTable *self)
{
    while (true)
    {
        size_t pivot_col;
        float smallest_obj_func_num = FLT_MAX;

        // calculating the column of the pivot
        for (size_t i = 0; i < self->cols - 1; ++i)
        {
            if (self->flat_table[i] < 0.0f && self->flat_table[i] < smallest_obj_func_num)
            {
                smallest_obj_func_num = self->flat_table[i];
                pivot_col = i;
            }
        }

        // if pivot_col was not found, end
        if (smallest_obj_func_num == FLT_MAX)
        {
            return;
        }

        size_t pivot_row;
        float best_pivot_row_div = FLT_MAX;
        // calculating the row of the pivot
        for (size_t i = 1; i < self->rows; ++i)
        {
            if (simplex_table_elem_val(self, i, pivot_col) <= 0.0f)
                continue;
            float curr_div = simplex_table_elem_val(self, i, self->cols - 1) / simplex_table_elem_val(self, i, pivot_col);
            if (curr_div > 0 && curr_div < best_pivot_row_div)
            {
                best_pivot_row_div = curr_div;
                pivot_row = i;
            }
        }

        // if pivot_row was not found, end
        if (best_pivot_row_div == FLT_MAX)
            return;

        simplex_table_print(self);
        // change basic var
        self->basic_vars_idxs[pivot_row - 1] = pivot_col;

        float pivot_val = simplex_table_elem_val(self, pivot_row, pivot_col);
        for (size_t i = 0; i < self->cols; ++i)
        {
            *simplex_table_elem_ptr(self, pivot_row, i) /= pivot_val;
        }

        for (size_t i = 0; i < self->rows; ++i)
        {
            if (i == pivot_row)
                continue;
            float mul_by = -simplex_table_elem_val(self, i, pivot_col);
            for (size_t j = 0; j < self->cols; ++j)
            {
                *simplex_table_elem_ptr(self, i, j) += simplex_table_elem_val(self, pivot_row, j) * mul_by;
            }
        }
        simplex_table_print(self);
    }
}

float simplex_table_elem_val(SimplexTable *self, size_t row, size_t column)
{
    return self->flat_table[row * self->cols + column];
}

float *simplex_table_elem_ptr(SimplexTable *self, size_t row, size_t column)
{
    return &(self->flat_table[row * self->cols + column]);
}

size_t simplex_table_size(SimplexTable *self)
{
    return self->rows * self->cols;
}

void simplex_table_print(SimplexTable *self)
{
    for (size_t i = 0; i < self->rows; ++i)
    {
        for (size_t j = 0; j < self->cols; ++j)
        {
            printf("%7.2f ", simplex_table_elem_val(self, i, j));
        }
        puts("");
    }
    puts("Basic variables:");
    for (size_t i = 0; i < self->rows - 1; ++i)
    {
        printf("%zd ", self->basic_vars_idxs[i]);
    }
    puts("");
    puts("---------------------");
}

void simplex_table_print_solution(SimplexTable *self)
{
    float *solution = calloc(0.0f, sizeof(float) * (self->cols - 1));
    for (size_t i = 0; i < self->rows - 1; ++i)
        solution[self->basic_vars_idxs[i] - 1] = self->flat_table[self->basic_vars_idxs[i]];

    printf("Optimal solution: <");
    for (size_t i = 0; i < self->cols - 1; ++i)
        printf("%f, ", solution[i]);
    printf(">, f(x)=%.2f", self->flat_table[self->cols - 1]);

    free(solution);
}

void simplex_table_destroy(SimplexTable *self)
{
    free(self->basic_vars_idxs);
    free(self->flat_table);
    free(self);
}

static bool simplex_table_is_var_in_basis(SimplexTable *self, size_t var_idx)
{
    for (size_t i = 0; i < self->rows - 1; ++i)
        if (self->basic_vars_idxs[i] == var_idx)
            return true;
    return false;
}