#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32
#include <conio.h>
#elif defined(__LINUX__)
#include <curses.h>
#endif

#include "simplex_table.h"

int main(int argc, char **argv)
{
    SimplexTable *simplex_table = simplex_table_create();
    simplex_table_solve(simplex_table);
    simplex_table_print_solution(simplex_table);
    simplex_table_destroy(simplex_table);
    getch();
    return EXIT_SUCCESS;
}