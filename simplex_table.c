#include "simplex_table.h"

typedef struct SimplexTable
{
    size_t rows_st;
    size_t cols_st;
    size_t *basic_vars;
    float *st_vars;
    float *costs;
} SimplexTable;

SimplexTable *simplex_table_create(void)
{
    SimplexTable *self = malloc(sizeof(SimplexTable));

    FILE *fhandler = fopen("priklad1.txt", "r");
    if (!fhandler)
        return nullptr;

    setlocale(LC_NUMERIC, "sk_SK.utf8");

    fscanf(fhandler, "%zd %zd %*s ", &(self->rows_st), &(self->cols_st));

    for (size_t i = 0; i < self->rows_st; ++i)
        fscanf(fhandler, "%*c ");

    ++self->cols_st;
    self->costs = malloc(sizeof(float) * (self->cols_st - 1));
    self->st_vars = malloc(sizeof(float) * self->rows_st * self->cols_st);
    self->basic_vars = malloc(sizeof(size_t) * self->rows_st);

    for (size_t i = 0; i < self->cols_st - 1; ++i)
        fscanf(fhandler, "%f", &(self->costs[i]));

    // loading the st vars flat table
    for (size_t i = 0; i < self->rows_st * self->cols_st; ++i)
        fscanf(fhandler, "%f", &(self->st_vars[i]));

    // loading the basic vars
    for (size_t i = 0; i < self->rows_st; ++i)
    {
        fscanf(fhandler, "%zd", &(self->basic_vars[i]));
        --self->basic_vars[i];
    }

    fclose(fhandler);
    return self;
}

void simplex_table_solve(SimplexTable *self)
{
    while (true)
    {
        size_t pivot_col;
        float smallest_reduced_cost = FLT_MAX;
        float reduced_cost;
        bool found_piv = false;
        // calculating the column of the pivot
        for (size_t i = 0; i < self->cols_st - 1; ++i)
        {
            reduced_cost = 0.0f;
            for (size_t j = 0; j < self->rows_st; ++j)
                reduced_cost += self->costs[self->basic_vars[j]] * simplex_table_elem_val(self, j, i);
            reduced_cost = self->costs[i] - reduced_cost;
            if (reduced_cost < 0 && reduced_cost < smallest_reduced_cost)
            {
                smallest_reduced_cost = reduced_cost;
                pivot_col = i;
                found_piv = true;
            }
        }

        // if pivot_col was not found, end
        if (!found_piv)
            return;

        size_t pivot_row;
        float smallest_ratio = FLT_MAX;
        found_piv = false;
        // calculating the row of the pivot
        for (size_t i = 0; i < self->rows_st; ++i)
        {
            if (simplex_table_elem_val(self, i, pivot_col) <= 0.0f)
                continue;
            float curr_div = simplex_table_elem_val(self, i, self->cols_st - 1) / simplex_table_elem_val(self, i, pivot_col);
            if (curr_div > 0 && curr_div < smallest_ratio)
            {
                smallest_ratio = curr_div;
                pivot_row = i;
                found_piv = true;
            }
        }

        // if pivot_row was not found, end
        if (!found_piv)
            return;

        simplex_table_print(self);
        // change basic var
        self->basic_vars[pivot_row] = pivot_col;

        float pivot_val = simplex_table_elem_val(self, pivot_row, pivot_col);
        for (size_t i = 0; i < self->cols_st; ++i)
        {
            *simplex_table_elem_ptr(self, pivot_row, i) /= pivot_val;
        }

        for (size_t i = 0; i < self->rows_st; ++i)
        {
            if (i == pivot_row)
                continue;
            float mul_by = -simplex_table_elem_val(self, i, pivot_col);
            for (size_t j = 0; j < self->cols_st; ++j)
            {
                *simplex_table_elem_ptr(self, i, j) += simplex_table_elem_val(self, pivot_row, j) * mul_by;
            }
        }
        simplex_table_print(self);
    }
}

float simplex_table_elem_val(SimplexTable *self, size_t row, size_t column)
{
    return self->st_vars[row * self->cols_st + column];
}

float *simplex_table_elem_ptr(SimplexTable *self, size_t row, size_t column)
{
    return &(self->st_vars[row * self->cols_st + column]);
}

size_t simplex_table_size(SimplexTable *self)
{
    return self->rows_st * self->cols_st;
}

float simplex_table_obj_func(SimplexTable *self)
{
    float obj_func = 0.0f;
    for (size_t i = 0; i < self->rows_st; ++i)
        obj_func += simplex_table_elem_val(self, i, self->cols_st - 1) * self->costs[self->basic_vars[i]];

    return obj_func;
}

void simplex_table_print(SimplexTable *self)
{
    for (size_t i = 0; i < self->rows_st; ++i)
    {
        for (size_t j = 0; j < self->cols_st; ++j)
        {
            printf("%7.2f ", simplex_table_elem_val(self, i, j));
        }
        puts("");
    }
    puts("Basic variables:");
    for (size_t i = 0; i < self->rows_st; ++i)
    {
        printf("%zd ", self->basic_vars[i]);
    }
    puts("");
    puts("---------------------");
}

void simplex_table_print_solution(SimplexTable *self)
{
    float *solution = calloc((self->cols_st - 1), sizeof(float));
    for (size_t i = 0; i < self->rows_st; ++i)
        solution[self->basic_vars[i]] = simplex_table_elem_val(self, i, self->cols_st - 1);

    printf("Optimal solution: <");
    for (size_t i = 0; i < self->cols_st - 1; ++i)
        printf("%.2f, ", solution[i]);
    printf(">, f(x)=%.2f", simplex_table_obj_func(self));

    free(solution);
}

void simplex_table_destroy(SimplexTable *self)
{
    free(self->basic_vars);
    free(self->st_vars);
    free(self);
}