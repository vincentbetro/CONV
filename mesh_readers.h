#include <stdio.h>
#include "StarCD.h"

int Ensight(char filename[], int mode, int &nn, Point **node, int &nb, char ***b_name, int **nt, int ****t_n, 
             int **nq, int ****q_n, int **ngon, int ****ngon_n, int &ntet, int ***tet_n, int &npyr, int ***pyr_n,
             int &npri, int ***pri_n, int &nhex, int ***hex_n, int &nply, int ****poly_n);
int Generic_mesh(char filename[], int mode, int &nn, Point **node, int &nb, char ***b_name, int **nt, int ****t_n, 
             int **nq, int ****q_n, int **ngon, int ****ngon_n, int &ntet, int ***tet_n, int &npyr, int ***pyr_n,
             int &npri, int ***pri_n, int &nhex, int ***hex_n, int &nply, int ****poly_n);
int Plotfile(char filename[], int mode, int &nn, Point **node, int &nb, char ***b_name, int **nt, int ****t_n, 
             int **nq, int ****q_n, int **ngon, int ****ngon_n, int &ntet, int ***tet_n, int &npyr, int ***pyr_n,
             int &npri, int ***pri_n, int &nhex, int ***hex_n, int &nply, int ****poly_n);
//int StarCD_read(char *fname, int &nnodes, Point **nodes, int &nb, char ***b_name,
//                int **ntri, int ****tri_conn,
//                int **nquad, int ****quad_conn,
//                int **ngon, int ****gon_conn,
//                int &nvol, char ***vol_name, 
//                int &ntet, int ***tet_conn, int **tet_vc,
//                int &npyr, int ***pyr_conn, int **pyr_vc,
//                int &npri, int ***pri_conn, int **pri_vc,
//                int &nhex, int ***hex_conn, int **hex_vc,
//                int &nply, int ****poly_conn, int **poly_vc);
//int StarCD_write(char *fname, const int nnodes, Point *nodes,
//                 const int nb, char **b_name,
//                 const int * const ntri, int ***tri_conn,
//                 const int * const nquad, int ***quad_conn,
//                 const int * const ngon, int ***gon_conn,
//                 const int nvol, char **vol_name, 
//                 const int ntet, int **tet_conn, const int * const tet_vc,
//                 const int npyr, int **pyr_conn, const int * const pyr_vc,
//                 const int npri, int **pri_conn, const int * const pri_vc,
//                 const int nhex, int **hex_conn, const int * const hex_vc,
