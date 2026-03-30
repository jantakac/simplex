#include "simplex_table.h"

typedef struct SimplexTable
{
    size_t rows_st;
    size_t cols_st;
    size_t *basic_vars;
    float *st_vars;
    float *costs;
    bool minimizing;
} SimplexTable;

static inline size_t calc_pivot_column(SimplexTable *self, bool *out_found);
static inline size_t calc_pivot_row(SimplexTable *self, size_t pivot_col, bool *out_found);
static inline void apply_pivot_transform(SimplexTable *self, size_t pivot_row, size_t pivot_col);

SimplexTable *simplex_table_create(void)
{
    SimplexTable *self = malloc(sizeof(SimplexTable));

    FILE *fhandler = fopen("priklad1.txt", "r");
    if (!fhandler)
        return nullptr;

    setlocale(LC_NUMERIC, "sk_SK.utf8");

    char optimization_direction[4];
    fscanf(fhandler, "%zd %zd %3s ", &(self->rows_st), &(self->cols_st), optimization_direction);

    // mi'n' vs ma'x'
    if (optimization_direction[2] == 'n')
        self->minimizing = true;
    else
        self->minimizing = false;

    for (size_t i = 0; i < self->rows_st; ++i)
        fscanf(fhandler, "%*c ");

    ++self->cols_st;
    self->costs = malloc(sizeof(float) * (self->cols_st - 1));
    self->st_vars = malloc(sizeof(float) * self->rows_st * self->cols_st);
    self->basic_vars = malloc(sizeof(size_t) * self->rows_st);

    // loading the objective function costs
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
        simplex_table_print(self);

        bool found_piv;
        size_t pivot_col = calc_pivot_column(self, &found_piv);
        // if pivot_col was not found, end
        if (!found_piv)
            return;

        size_t pivot_row = calc_pivot_row(self, pivot_col, &found_piv);
        // if pivot_row was not found, end
        if (!found_piv)
            return;

        apply_pivot_transform(self, pivot_row, pivot_col);
        simplex_table_print(self);
    }
}

void simplex_table_dual_solve(SimplexTable *self)
{
    while (true)
    {
        simplex_table_print(self);

        bool found_piv = false;
        float best_pivot_row_val = FLT_MAX;
        size_t pivot_row;

        for (size_t i = 0; i < self->rows_st; ++i)
        {
            float curr_val = simplex_table_elem_val(self, i, self->cols_st - 1);
            if (curr_val < 0 && curr_val < best_pivot_row_val)
            {
                found_piv = true;
                best_pivot_row_val = curr_val;
                pivot_row = i;
            }
        }

        if (!found_piv)
            return;

        found_piv = false;
        float best_pivot_col_val = FLT_MAX;
        size_t pivot_col;

        for (size_t i = 0; i < self->cols_st; ++i)
        {
            if (simplex_table_elem_val(self, pivot_row, i) >= 0)
                continue;

            float curr_val = abs(self->costs[i] / simplex_table_elem_val(self, pivot_row, i));
            if (curr_val < best_pivot_col_val)
            {
                found_piv = true;
                best_pivot_col_val = curr_val;
                pivot_col = i;
            }
        }

        if (!found_piv)
            return;

        apply_pivot_transform(self, pivot_row, pivot_col);
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
            printf("%7.2f ", simplex_table_elem_val(self, i, j));
        puts("");
    }
    puts("Basic variables:");
    for (size_t i = 0; i < self->rows_st; ++i)
        printf("%zd ", self->basic_vars[i] + 1);
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
    free(self->costs);
    free(self->basic_vars);
    free(self->st_vars);
    free(self);
}

static inline size_t calc_pivot_column(SimplexTable *self, bool *out_found)
{
    *out_found = false;
    size_t pivot_col;
    float best_reduced_cost;
    float reduced_cost;

    if (self->minimizing)
        best_reduced_cost = FLT_MAX;
    else
        best_reduced_cost = -FLT_MAX;

    for (size_t i = 0; i < self->cols_st - 1; ++i)
    {
        reduced_cost = 0.0f;
        for (size_t j = 0; j < self->rows_st; ++j)
            reduced_cost += self->costs[self->basic_vars[j]] * simplex_table_elem_val(self, j, i);
        reduced_cost = self->costs[i] - reduced_cost;
        bool is_better = false;
        if (self->minimizing)
        {
            if (reduced_cost < 0 && reduced_cost < best_reduced_cost)
                is_better = true;
        }
        else
        {
            if (reduced_cost > 0 && reduced_cost > best_reduced_cost)
                is_better = true;
        }
        if (is_better)
        {
            best_reduced_cost = reduced_cost;
            pivot_col = i;
            *out_found = true;
        }
    }
    return pivot_col;
}

static inline size_t calc_pivot_row(SimplexTable *self, size_t pivot_col, bool *out_found)
{
    size_t pivot_row;
    float smallest_ratio = FLT_MAX;
    *out_found = false;
    // calculating the row of the pivot
    for (size_t i = 0; i < self->rows_st; ++i)
    {
        if (simplex_table_elem_val(self, i, pivot_col) <= 0.0f)
            continue;
        float curr_div = simplex_table_elem_val(self, i, self->cols_st - 1) / simplex_table_elem_val(self, i, pivot_col);
        if (curr_div >= 0 && curr_div < smallest_ratio)
        {
            smallest_ratio = curr_div;
            pivot_row = i;
            *out_found = true;
        }
    }
    return pivot_row;
}

static inline void apply_pivot_transform(SimplexTable *self, size_t pivot_row, size_t pivot_col)
{
    // change basic var
    self->basic_vars[pivot_row] = pivot_col;

    float pivot_val = simplex_table_elem_val(self, pivot_row, pivot_col);
    for (size_t i = 0; i < self->cols_st; ++i)
        *simplex_table_elem_ptr(self, pivot_row, i) /= pivot_val;

    for (size_t i = 0; i < self->rows_st; ++i)
    {
        if (i == pivot_row)
            continue;
        float mul_by = -simplex_table_elem_val(self, i, pivot_col);
        for (size_t j = 0; j < self->cols_st; ++j)
            *simplex_table_elem_ptr(self, i, j) += simplex_table_elem_val(self, pivot_row, j) * mul_by;
    }
}