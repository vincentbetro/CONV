#include <stdio.h>

#ifndef _mesh_io
#define _mesh_io
int Read_Mesh(char sname[], char **&b_name, int &nn, double *&x, double *&y, double *&z,
               int &nb, int *&nt, int *&nq, int *&ngon, int ***&t_n, int ***&q_n, int ***&ngon_n,
               int &ntet, int **&tet_n, int &npyr, int **&pyr_n,
               int &npri, int **&pri_n, int &nhex, int **&hex_n, int &nply, int ***&poly_n);

int Read_Mesh(char sname[], char **&b_name, int &nn, double *&x, double *&y, double *&z,
               int &nb, int *&nt, int *&nq, int *&ngon, int ***&t_n, int ***&q_n, int ***&ngon_n,
               int &ntet, int **&tet_n, int &npyr, int **&pyr_n,
               int &npri, int **&pri_n, int &nhex, int **&hex_n, int &nply, int ***&poly_n,
               int &nvol, char **&vol_name, int *&tet_vc, int *&pyr_vc, int *&pri_vc, int *&hex_vc, int *&poly_vc);

int Write_Mesh(char sname[], char **&b_name, int &nn, double *&x, double *&y, double *&z,
               int &nb, int *&nt, int *&nq, int *&ngon, int ***&t_n, int ***&q_n, int ***&ngon_n,
               int &ntet, int **&tet_n, int &npyr, int **&pyr_n,
               int &npri, int **&pri_n, int &nhex, int **&hex_n, int &nply, int ***&poly_n);

int Write_Mesh(char sname[], char **&b_name, int &nn, double *&x, double *&y, double *&z,
               int &nb, int *&nt, int *&nq, int *&ngon, int ***&t_n, int ***&q_n, int ***&ngon_n,
               int &ntet, int **&tet_n, int &npyr, int **&pyr_n,
               int &npri, int **&pri_n, int &nhex, int **&hex_n, int &nply, int ***&poly_n,
               int &nvol, char **&vol_name, int *&tet_vc, int *&pyr_vc, int *&pri_vc, int *&hex_vc, int *&poly_vc);

int displaced(char dname[], int *&dnode, double *&dxn, double *&dyn, double *&dzn,
              int nn, double *x, double *y, double*z,
              int nb, int *nt, int *nq, int *ngon, int ***t_n, int ***q_n, int ***ngon_n);

#endif
