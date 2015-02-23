#include <stdio.h>
#ifndef _create_hybrid_maps
#define _create_hybrid_maps

void create_hybrid_maps(FILE *out_f, int num_procs, int my_rank, int nn, int nb, int ntet, int npyr, int npri, int nhex, int nply,
                        int **tet_n, int **pyr_n, int **pri_n, int **hex_n, int ***poly_n,
                        int *nt, int ***t_n, int *nq, int ***q_n, int *ngon, int ***ngon_n,
                        int **nmap, int **tri_map, int **quad_map, int **ngon_map,
                        int *tet_map, int *pyr_map, int *pri_map, int *hex_map, int *poly_map);

#endif
