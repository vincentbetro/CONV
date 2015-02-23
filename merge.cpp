#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "Point.h"
#include "Vector.h"
#include "Util.h"
#include "Linked_List.h"
#include "List.h"
#include "merge.h"
#include "convert.h"
#include "PList.h"
#include "Octree_Storage.h"
#include "mesh_io.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a,b) ((b) < 0.0 ? -ABS(a) : ABS(a))

#define CHUNK 100000

int merge(char fname1[], char fname2[], char mname[], double tol)
{
  int b, c, f, i, j, k, m, n, nf, n0, n1, n2, n3, m0, m1, m2, m3;
  FILE *fp;
  int error = 0;
 
  // 1st mesh
  int nn, nb, ntet, npyr, npri, nhex, nply, nvc;
  int *nt, *nq, *ngon, ***t_n, ***q_n, ***ngon_n;
  int **tet_n, **pyr_n, **pri_n, **hex_n, ***poly_n;
  int *tet_vc, *pyr_vc, *pri_vc, *hex_vc, *poly_vc;
  double *x, *y, *z;
  Point *node;
  char **b_name, **vc_name;
  nn = nb = nvc = ntet = npyr = npri = nhex = nply = 0;
  // set all pointers to null
  nt = 0;
  nq = 0;
  ngon = 0;
  t_n = 0;
  q_n = 0;
  ngon_n = 0;
  tet_n = 0;
  pyr_n = 0;
  pri_n = 0;
  hex_n = 0;
  poly_n = 0;
  tet_vc = 0;
  pyr_vc = 0;
  pri_vc = 0;
  hex_vc = 0;
  poly_vc = 0;
  node = 0;

  //error = Read_Mesh(fname1,b_name,nn,x,y,z,nb,nt,nq,ngon,t_n,q_n,ngon_n,ntet,tet_n,npyr,pyr_n,npri,pri_n,nhex,hex_n,nply,poly_n);
  error = Read_Mesh(fname1,b_name,nn,x,y,z,nb,nt,nq,ngon,t_n,q_n,ngon_n,ntet,tet_n,npyr,pyr_n,npri,pri_n,nhex,hex_n,nply,poly_n,
            nvc,vc_name,tet_vc,pyr_vc,pri_vc,hex_vc,poly_vc);

  node = (Point*)malloc(nn*sizeof(Point));
  for (n=0; n < nn; n++)
    node[n] = Point(x[n],y[n],z[n]);

  free(x);
  free(y);
  free(z);

  printf("\n\n1st mesh statistics:");
  printf("\n # of nodes = %d",nn);
  if (ntet > 0) printf("\n # of tetrahedra = %d",ntet);
  if (npyr > 0) printf("\n # of pyramids   = %d",npyr);
  if (npri > 0) printf("\n # of prisms     = %d",npri);
  if (nhex > 0) printf("\n # of hexahedra  = %d",nhex);
  if (nply > 0) printf("\n # of polyhedra  = %d",nply);
  printf("\n # of volume conditions = %i",nvc);
  for (i=0; i < nvc; i++)
    printf("\nVolume %i name = %s",i+1,vc_name[i]);
  printf("\n # of boundaries = %d",nb);
  int nttot = 0;
  int nqtot = 0;
  int ngtot = 0;
  for (b=0; b < nb; b++)
  {
    printf("\nBoundary %d name = %s",b+1,b_name[b]);
    if (  nt[b] > 0) printf("\n  # of triangles      = %d",nt[b]);
    if (  nq[b] > 0) printf("\n  # of quadrilaterals = %d",nq[b]);
    if (ngon[b] > 0) printf("\n  # of polygons       = %d",ngon[b]);
    nttot += nt[b];
    nqtot += nq[b];
    ngtot += ngon[b];
  }
  if (nttot > 0) printf("\nTotal number of surface triangles      = %d",nttot);
  if (nqtot > 0) printf("\nTotal number of surface quadrilaterals = %d",nqtot);
  if (ngtot > 0) printf("\nTotal number of surface polygons       = %d",ngtot);
  printf("\n\n");
  
  List *vlist = new List();
  for (i=0; i < ntet; i++)
    vlist->Check_List(tet_vc[i]);
  for (i=0; i < npyr; i++)
    vlist->Check_List(pyr_vc[i]);
  for (i=0; i < npri; i++)
    vlist->Check_List(pri_vc[i]);
  for (i=0; i < nhex; i++)
    vlist->Check_List(hex_vc[i]);
  for (i=0; i < nply; i++)
    vlist->Check_List(poly_vc[i]);

  for (i=0; i < ntet; i++)
    tet_vc[i] = vlist->Index(tet_vc[i]);
  for (i=0; i < npyr; i++)
    pyr_vc[i] = vlist->Index(pyr_vc[i]);
  for (i=0; i < npri; i++)
    pri_vc[i] = vlist->Index(pri_vc[i]);
  for (i=0; i < nhex; i++)
    hex_vc[i] = vlist->Index(hex_vc[i]);
  for (i=0; i < nply; i++)
    poly_vc[i] = vlist->Index(poly_vc[i]);
  delete vlist;

  // 2nd mesh
  int nn2, nb2, nvc2, ntet2, npyr2, npri2, nhex2, nply2;
  int *nt2, *nq2, *ngon2, ***t_n2, ***q_n2, ***ngon_n2;
  int **tet_n2, **pyr_n2, **pri_n2, **hex_n2, ***poly_n2;
  int *tet_vc2, *pyr_vc2, *pri_vc2, *hex_vc2, *poly_vc2;
  Point *node2;
  char **b_name2, **vc_name2;
  nn2 = nb2 = nvc2 = ntet2 = npyr2 = npri2 = nhex2 = nply2 = 0;
  // set all pointers to null
  nt2 = 0;
  nq2 = 0;
  ngon2 = 0;
  t_n2 = 0;
  q_n2 = 0;
  ngon_n2 = 0;
  tet_n2 = 0;
  pyr_n2 = 0;
  pri_n2 = 0;
  hex_n2 = 0;
  poly_n2 = 0;
  tet_vc2 = 0;
  pyr_vc2 = 0;
  pri_vc2 = 0;
  hex_vc2 = 0;
  poly_vc2 = 0;
  node2 = 0;

  //error = Read_Mesh(fname2,b_name2,nn2,x,y,z,nb2,nt2,nq2,ngon2,t_n2,q_n2,ngon_n2,ntet2,tet_n2,npyr2,pyr_n2,npri2,pri_n2,nhex2,hex_n2,nply2,poly_n2);
  error = Read_Mesh(fname2,b_name2,nn2,x,y,z,nb2,nt2,nq2,ngon2,t_n2,q_n2,ngon_n2,ntet2,tet_n2,npyr2,pyr_n2,npri2,pri_n2,nhex2,hex_n2,nply2,poly_n2,
            nvc2,vc_name2,tet_vc2,pyr_vc2,pri_vc2,hex_vc2,poly_vc2);

  node2 = (Point*)malloc(nn2*sizeof(Point));
  for (n=0; n < nn2; n++)
    node2[n] = Point(x[n],y[n],z[n]);

  free(x);
  free(y);
  free(z);

  printf("\n\n2nd mesh statistics:");
  printf("\n # of nodes = %d",nn2);
  if (ntet > 0) printf("\n # of tetrahedra = %d",ntet2);
  if (npyr > 0) printf("\n # of pyramids   = %d",npyr2);
  if (npri > 0) printf("\n # of prisms     = %d",npri2);
  if (nhex > 0) printf("\n # of hexahedra  = %d",nhex2);
  if (nply > 0) printf("\n # of polyhedra  = %d",nply2);
  printf("\n # of volume conditions = %i",nvc2);
  for (i=0; i < nvc2; i++)
    printf("\nVolume %i name = %s",i+1,vc_name2[i]);
  printf("\n # of boundaries = %d",nb2);
  nttot = 0;
  nqtot = 0;
  ngtot = 0;
  for (b=0; b < nb2; b++)
  {
    printf("\nBoundary %d name = %s",b+1,b_name2[b]);
    if (  nt2[b] > 0) printf("\n  # of triangles      = %d",nt2[b]);
    if (  nq2[b] > 0) printf("\n  # of quadrilaterals = %d",nq2[b]);
    if (ngon2[b] > 0) printf("\n  # of polygons       = %d",ngon2[b]);
    nttot += nt2[b];
    nqtot += nq2[b];
    ngtot += ngon2[b];
  }
  if (nttot > 0) printf("\nTotal number of surface triangles      = %d",nttot);
  if (nqtot > 0) printf("\nTotal number of surface quadrilaterals = %d",nqtot);
  if (ngtot > 0) printf("\nTotal number of surface polygons       = %d",ngtot);
  printf("\n\n");

  List *vlist2 = new List();
  for (i=0; i < ntet2; i++)
    vlist2->Check_List(tet_vc2[i]);
  for (i=0; i < npyr2; i++)
    vlist2->Check_List(pyr_vc2[i]);
  for (i=0; i < npri2; i++)
    vlist2->Check_List(pri_vc2[i]);
  for (i=0; i < nhex2; i++)
    vlist2->Check_List(hex_vc2[i]);
  for (i=0; i < nply2; i++)
    vlist2->Check_List(poly_vc2[i]);

  for (i=0; i < ntet2; i++)
    tet_vc2[i] = vlist2->Index(tet_vc2[i]);
  for (i=0; i < npyr2; i++)
    pyr_vc2[i] = vlist2->Index(pyr_vc2[i]);
  for (i=0; i < npri2; i++)
    pri_vc2[i] = vlist2->Index(pri_vc2[i]);
  for (i=0; i < nhex2; i++)
    hex_vc2[i] = vlist2->Index(hex_vc2[i]);
  for (i=0; i < nply2; i++)
    poly_vc2[i] = vlist2->Index(poly_vc2[i]);
  delete vlist2;

  // add mesh 2 to mesh 1
  vc_name = (char**)realloc((void*)vc_name,(nvc+nvc2)*sizeof(char*));
  for (i=0; i < nvc2; i++)
  {
    vc_name[i+nvc] = (char*)malloc(33*sizeof(char));
    sprintf(vc_name[i+nvc],"%-32s",vc_name2[i]);
    free(vc_name2[i]);
  }
  free(vc_name2);

  node = (Point*)realloc((void*)node,(nn+nn2)*sizeof(Point));
  for (n=0; n < nn2; n++)
    node[n+nn] = node2[n];
  free(node2); node2=0;
  if (ntet2 > 0)
  {
    tet_vc = (int*)realloc((void*)tet_vc,(ntet+ntet2)*sizeof(int));
    tet_n = (int**)realloc((void*)tet_n,(ntet+ntet2)*sizeof(int*));
    for (c=0; c < ntet2; c++)
    {
      tet_n[c+ntet] = (int*)malloc(4*sizeof(int));
      tet_n[c+ntet][0] = tet_n2[c][0]+nn;
      tet_n[c+ntet][1] = tet_n2[c][1]+nn;
      tet_n[c+ntet][2] = tet_n2[c][2]+nn;
      tet_n[c+ntet][3] = tet_n2[c][3]+nn;
      tet_vc[c+ntet] = tet_vc2[c]+nvc;
      free(tet_n2[c]);
    }
    free(tet_vc2); tet_vc2=0;
    free(tet_n2); tet_n2=0;
    ntet += ntet2;
  }
  if (npyr2 > 0)
  {
    pyr_vc = (int*)realloc((void*)pyr_vc,(npyr+npyr2)*sizeof(int));
    pyr_n = (int**)realloc((void*)pyr_n,(npyr+npyr2)*sizeof(int*));
    for (c=0; c < npyr2; c++)
    {
      pyr_n[c+npyr] = (int*)malloc(5*sizeof(int));
      pyr_n[c+npyr][0] = pyr_n2[c][0]+nn;
      pyr_n[c+npyr][1] = pyr_n2[c][1]+nn;
      pyr_n[c+npyr][2] = pyr_n2[c][2]+nn;
      pyr_n[c+npyr][3] = pyr_n2[c][3]+nn;
      pyr_n[c+npyr][4] = pyr_n2[c][4]+nn;
      pyr_vc[c+npyr] = pyr_vc2[c]+nvc;
      free(pyr_n2[c]);
    }
    free(pyr_vc2); pyr_vc2=0;
    free(pyr_n2); pyr_n2=0;
    npyr += npyr2;
  }
  if (npri2 > 0)
  {
    pri_vc = (int*)realloc((void*)pri_vc,(npri+npri2)*sizeof(int));
    pri_n = (int**)realloc((void*)pri_n,(npri+npri2)*sizeof(int*));
    for (c=0; c < npri2; c++)
    {
      pri_n[c+npri] = (int*)malloc(6*sizeof(int));
      pri_n[c+npri][0] = pri_n2[c][0]+nn;
      pri_n[c+npri][1] = pri_n2[c][1]+nn;
      pri_n[c+npri][2] = pri_n2[c][2]+nn;
      pri_n[c+npri][3] = pri_n2[c][3]+nn;
      pri_n[c+npri][4] = pri_n2[c][4]+nn;
      pri_n[c+npri][5] = pri_n2[c][5]+nn;
      pri_vc[c+npri] = pri_vc2[c]+nvc;
      free(pri_n2[c]);
    }
    free(pri_vc2); pri_vc2=0;
    free(pri_n2); pri_n2=0;
    npri += npri2;
  }
  if (nhex2 > 0)
  {
    hex_vc = (int*)realloc((void*)hex_vc,(nhex+nhex2)*sizeof(int));
    hex_n = (int**)realloc((void*)hex_n,(nhex+nhex2)*sizeof(int*));
    for (c=0; c < nhex2; c++)
    {
      hex_n[c+nhex] = (int*)malloc(8*sizeof(int));
      hex_n[c+nhex][0] = hex_n2[c][0]+nn;
      hex_n[c+nhex][1] = hex_n2[c][1]+nn;
      hex_n[c+nhex][2] = hex_n2[c][2]+nn;
      hex_n[c+nhex][3] = hex_n2[c][3]+nn;
      hex_n[c+nhex][4] = hex_n2[c][4]+nn;
      hex_n[c+nhex][5] = hex_n2[c][5]+nn;
      hex_n[c+nhex][6] = hex_n2[c][6]+nn;
      hex_n[c+nhex][7] = hex_n2[c][7]+nn;
      hex_vc[c+nhex] = hex_vc2[c]+nvc;
      free(hex_n2[c]);
    }
    free(hex_vc2); hex_vc2=0;
    free(hex_n2); hex_n2=0;
    nhex += nhex2;
  }
  if (nply2 > 0)
  {
    poly_vc = (int*)realloc((void*)poly_vc,(nply+nply2)*sizeof(int));
    poly_n = (int***)realloc((void*)poly_n,(nply+nply2)*sizeof(int**));
    for (c=0; c < nply2; c++)
    {
      nf = poly_n2[c][0][0];
      poly_n[c+nply] = (int**)malloc((nf+1)*sizeof(int*));
      poly_n[c+nply][0] = (int*)malloc(sizeof(int));
      poly_n[c+nply][0][0] = nf;
      for (f=1; f <= nf; f++)
      {
        i=poly_n2[c][f][0];
        poly_n[c+nply][f] = (int*)malloc((i+1)*sizeof(int));
        poly_n[c+nply][f][0] = i;
        for (j=1; j <= i; j++)
          poly_n[c+nply][f][j] = poly_n2[c][f][j]+nn;
      }
      poly_vc[c+nply] = poly_vc2[c]+nvc;
        
      for (j=poly_n2[f][0][0]; j >= 0; j--)
        if (poly_n2[f][j] != 0) free(poly_n2[f][j]);
      if (poly_n2[f] != 0) free(poly_n2[f]);
    }
    free(poly_vc2); poly_vc2=0;
    free(poly_n2); poly_n2 = 0;
    nply += nply2;
  }
  nvc += nvc2;
  
  nt = (int*)realloc((void*)nt,(nb+nb2)*sizeof(int));
  nq = (int*)realloc((void*)nq,(nb+nb2)*sizeof(int));
  ngon = (int*)realloc((void*)ngon,(nb+nb2)*sizeof(int));
  t_n = (int***)realloc((void*)t_n,(nb+nb2)*sizeof(int**));
  q_n = (int***)realloc((void*)q_n,(nb+nb2)*sizeof(int**));
  ngon_n = (int***)realloc((void*)ngon_n,(nb+nb2)*sizeof(int**));
  b_name = (char**)realloc((void*)b_name,(nb+nb2)*sizeof(char*));
  for (b=0; b < nb2; b++)
  {
    nt[b+nb] = 0;
    nq[b+nb] = 0;
    ngon[b+nb] = 0;
    t_n[b+nb] = 0;
    q_n[b+nb] = 0;
    ngon_n[b+nb] = 0;
    b_name[b+nb] = (char*)malloc(33*sizeof(char));
    sprintf(b_name[b+nb],"%-32s",b_name2[b]);
    free(b_name2[b]);
    if (nt2[b] > 0)
    {
      t_n[b+nb] = (int**)malloc(nt2[b]*sizeof(int*));
      for (c=0; c < nt2[b]; c++)
      {
        t_n[b+nb][c] = (int*)malloc(3*sizeof(int));
        t_n[b+nb][c][0] = t_n2[b][c][0]+nn;
        t_n[b+nb][c][1] = t_n2[b][c][1]+nn;
        t_n[b+nb][c][2] = t_n2[b][c][2]+nn;
        free(t_n2[b][c]);
      }
      free(t_n2[b]);
      nt[b+nb] = nt2[b];
    }
    if (nq2[b] > 0)
    {
      q_n[b+nb] = (int**)malloc(nq2[b]*sizeof(int*));
      for (c=0; c < nq2[b]; c++)
      {
        q_n[b+nb][c] = (int*)malloc(4*sizeof(int));
        q_n[b+nb][c][0] = q_n2[b][c][0]+nn;
        q_n[b+nb][c][1] = q_n2[b][c][1]+nn;
        q_n[b+nb][c][2] = q_n2[b][c][2]+nn;
        q_n[b+nb][c][3] = q_n2[b][c][3]+nn;
        free(q_n2[b][c]);
      }
      free(q_n2[b]);
      nq[b+nb] = nq2[b];
    }
    if (ngon2[b] > 0)
    {
      ngon_n[b+nb] = (int**)malloc(ngon2[b]*sizeof(int*));
      for (c=0; c < ngon2[b]; c++)
      {
        ngon_n[b+nb][c] = (int*)malloc((ngon_n2[b][c][0]+1)*sizeof(int));
        ngon_n[b+nb][c][0] = ngon_n2[b][c][0];
        for (i=1; i <= ngon_n2[b][c][0]; i++)
          ngon_n[b+nb][c][i] = ngon_n2[b][c][i]+nn;
        free(ngon_n2[b][c]);
      }
      free(ngon_n2[b]);
      ngon[b+nb] = ngon2[b];
    }
  }
  free(b_name2);
  free(nt2); nt2=0;
  free(nq2); nq2=0;
  free(ngon2); ngon2=0;
  free(t_n2); t_n2=0;
  free(q_n2); q_n2=0;
  free(ngon_n2); ngon_n2=0;

  nb += nb2;
  nn += nn2;

  // tag boundary nodes
  int *map = new int[nn];
  for (n=0; n < nn; n++)
    map[n] = -1;

  double dsmn = 1.0e20;
  for (b=0; b < nb; b++)
  {
    for (c=0; c < nt[b]; c++)
      for (i=0; i < 3; i++)
      {
        n0 = t_n[b][c][i];
        n1 = t_n[b][c][(i+1)%3];
        dsmn = MIN(dsmn,distance(node[n0],node[n1]));
        map[t_n[b][c][i]] = 1;
      }
    for (c=0; c < nq[b]; c++)
      for (i=0; i < 4; i++)
      {
        n0 = q_n[b][c][i];
        n1 = q_n[b][c][(i+1)%4];
        dsmn = MIN(dsmn,distance(node[n0],node[n1]));
        map[q_n[b][c][i]] = 1;
      }
    for (c=0; c < ngon[b]; c++)
      for (i=1; i <= ngon_n[b][c][0]; i++)
      {
        n0 = ngon_n[b][c][i];
        if (i == ngon_n[b][c][0])
          n1 = 1;
        else
          n1 = ngon_n[b][c][i+1];
        dsmn = MIN(dsmn,distance(node[n0],node[n1]));
        map[ngon_n[b][c][i]] = 1;
      }
  }

  printf("\nComputed minimum spacing between boundary nodes = %lg",dsmn);
  tol = MIN(tol,dsmn*0.1);
  printf("\nTolerance used for node comparison = %lg",tol);

  double olo[3], ohi[3];
  List *bnode = new List();
  for (n=0; n < nn; n++)
    if (map[n] == 1) 
    {
      map[n] = n;
      bnode->Add_To_List(n);
      olo[0] = MIN(olo[0],node[n][0]);
      olo[1] = MIN(olo[1],node[n][1]);
      olo[2] = MIN(olo[2],node[n][2]);
      ohi[0] = MAX(ohi[0],node[n][0]);
      ohi[1] = MAX(ohi[1],node[n][1]);
      ohi[2] = MAX(ohi[2],node[n][2]);
    }

  printf("\nNumber of boundary nodes = %i",bnode->max);

  olo[0] -= dsmn;
  olo[1] -= dsmn;
  olo[2] -= dsmn;
  ohi[0] += dsmn;
  ohi[1] += dsmn;
  ohi[2] += dsmn;

  Octree_Storage *dummy=0;
  Octree_Storage *Octree_root = new Octree_Storage((Octree_Storage*)dummy,olo,ohi,25,dsmn*0.1);

  //typedef struct {
  //  int n;
  //  Point p;
  //} BNODE;

  //BNODE *bnds = new BNODE[bnode->max];
  
  void *ptr;
  dsmn *= 0.25;
  for (i=0; i < bnode->max; i++)
  {
    n = bnode->list[i];

    olo[0] = node[n][0]-dsmn;
    olo[1] = node[n][1]-dsmn;
    olo[2] = node[n][2]-dsmn;
    ohi[0] = node[n][0]+dsmn;
    ohi[1] = node[n][1]+dsmn;
    ohi[2] = node[n][2]+dsmn;

    //bnds[i].n = n;
    //bnds[i].p = node[n];

    //ptr = (void*)&bnds[i];
    ptr = (void*)&map[n];
    Octree_root->Store_In_Octree(ptr,olo,ohi);
  }

  PList *plist = new PList();
  plist->Redimension(bnode->max);

  double pt[3], ptol[3];
  ptol[0] = dsmn;
  ptol[1] = dsmn;
  ptol[2] = dsmn;

  int tenth = bnode->max/10;
  int hundredth = bnode->max/100;
  printf("\nComparing boundary nodes. '.' every %i node\n",hundredth);

  //for (i=0; i < bnode->max-1; i++)
  //{
  //  n = bnode->list[i];
  //  if (((i+1)%hundredth) == 0) { printf("."); fflush(stdout); }
  //  if (((i+1)%tenth) == 0) { printf("\n"); fflush(stdout); }
  //  if (map[n] < 0 || map[n] != n) continue;

  //  for (j=i+1; j < bnode->max; j++)
  //  {
  //    m = bnode->list[j];
  //    if (map[m] < 0 || map[m] != m) continue;
  //    if (distance(node[n],node[m]) < tol)
  //      map[m] = n;
  //  }
  //}

  //BNODE *bptr;
  int *iptr;
  for (i=0; i < bnode->max; i++)
  {
    n = bnode->list[i];
    if (((i+1)%hundredth) == 0) { printf("."); fflush(stdout); }
    if (((i+1)%tenth) == 0) { printf("\n"); fflush(stdout); }
    if (map[n] < 0 || map[n] != n) continue;

    pt[0] = node[n][0];
    pt[1] = node[n][1];
    pt[2] = node[n][2];

    plist->max = 0;
    Octree_root->retrieve_list(pt,ptol,0,plist);

    for (j=0; j < plist->max; j++)
    {
      //bptr=(BNODE*)plist->list[j];
      //m = bptr->n;
      //if (m <= n || map[m] < 0 || map[m] != m) continue;
      //if (distance(node[n],node[m]) < tol)
      //  map[m] = n;

      iptr=(int*)plist->list[j];
      if (*iptr < 0 || *iptr <= n) continue;
      if (distance(node[n],node[*iptr]) < tol)
        *iptr = n;
    }
  }

  //delete[] bnds;
  delete plist;
  delete Octree_root;

  m = 0;
  for (n=0; n < nn; n++)
  {
    if (map[n] == -1)
      map[n] = m++;
    else if (map[n] == n)
      map[n] = m++;
    else
      map[n] = map[map[n]];
  }

  printf("\nOld number of nodes = %i",nn);
  printf("\nNew number of nodes = %i",m);
  fflush(stdout);

  // modify connectivities
  for (c=0; c < ntet; c++)
    for (i=0; i < 4; i++)
      tet_n[c][i] = map[tet_n[c][i]];
  for (c=0; c < npyr; c++)
    for (i=0; i < 5; i++)
      pyr_n[c][i] = map[pyr_n[c][i]];
  for (c=0; c < npri; c++)
    for (i=0; i < 6; i++)
      pri_n[c][i] = map[pri_n[c][i]];
  for (c=0; c < nhex; c++)
    for (i=0; i < 8; i++)
      hex_n[c][i] = map[hex_n[c][i]];
  for (c=0; c < nply; c++)
    for (f=1; f <= poly_n[c][0][0]; f++)
      for (i=1; i <= poly_n[c][f][0]; i++)
        poly_n[c][f][i] = map[poly_n[c][f][i]];
  for (b=0; b < nb; b++)
  {
    for (c=0; c < nt[b]; c++)
      for (i=0; i < 3; i++)
        t_n[b][c][i] = map[t_n[b][c][i]];
    for (c=0; c < nq[b]; c++)
      for (i=0; i < 4; i++)
        q_n[b][c][i] = map[q_n[b][c][i]];
    for (c=0; c < ngon[b]; c++)
      for (i=1; i <= ngon_n[b][c][0]; i++)
        ngon_n[b][c][i] = map[ngon_n[b][c][i]];
  }

  // now delete duplicate nodes
  for (n=0; n < nn; n++)
    if (map[n] != n)
      node[map[n]] = node[n];
  nn = m;
  node = (Point*)realloc((void*)node,nn*sizeof(Point));

  delete[] map;
  delete bnode;

  // now identify duplicate boundary faces

  List **hash;
  hash = new List*[nn];
  for (n=0; n < nn; n++)
    hash[n] = new List();

  bool match;

  // process triangles
  int total = 0;
  int mb, mc;
  for (b=0; b < nb; b++)
  {
    for (c=0; c < nt[b]; c++)
    {
      n0 = t_n[b][c][0];
      n1 = t_n[b][c][1];
      n2 = t_n[b][c][2];
      match = false;
      for (i=0; i < hash[n0]->max && !match; i++)
      {
        m = hash[n0]->list[i];
        if (m != (c+total) && hash[n1]->Is_In_List(m) && hash[n2]->Is_In_List(m))
        {
          match=true;
          t_n[b][c][0] = -1;
          t_n[b][c][1] = -1;
          t_n[b][c][2] = -1;
          for (k=mb=0; mb < nb; mb++)
          {
            if (k+nt[mb] > m)
              break;
            k+=nt[mb];
          }
          mc = m-k;
          hash[n0]->Delete_From_List(m);
          hash[n1]->Delete_From_List(m);
          hash[n2]->Delete_From_List(m);
          t_n[mb][mc][0] = -1;
          t_n[mb][mc][1] = -1;
          t_n[mb][mc][2] = -1;
        }
      }
      if (!match)
      {
        m = c+total;
        hash[n0]->Add_To_List(m);
        hash[n1]->Add_To_List(m);
        hash[n2]->Add_To_List(m);
      }
    }
    total += nt[b];
  }

  // process quads
  for (n=0; n < nn; n++)
    hash[n]->Redimension(0);

  total = 0;
  for (b=0; b < nb; b++)
  {
    for (c=0; c < nq[b]; c++)
    {
      n0 = q_n[b][c][0];
      n1 = q_n[b][c][1];
      n2 = q_n[b][c][2];
      n3 = q_n[b][c][3];
      match = false;
      for (i=0; i < hash[n0]->max && !match; i++)
      {
        m = hash[n0]->list[i];
        if (m != (c+total) && hash[n1]->Is_In_List(m) && hash[n2]->Is_In_List(m) && hash[n3]->Is_In_List(m))
        {
          match=true;
          q_n[b][c][0] = -1;
          q_n[b][c][1] = -1;
          q_n[b][c][2] = -1;
          q_n[b][c][3] = -1;
          for (k=mb=0; mb < nb; mb++)
          {
            if (k+nq[mb] > m)
              break;
            k+=nq[mb];
          }
          mc = m-k;
          hash[n0]->Delete_From_List(m);
          hash[n1]->Delete_From_List(m);
          hash[n2]->Delete_From_List(m);
          hash[n3]->Delete_From_List(m);
          q_n[mb][mc][0] = -1;
          q_n[mb][mc][1] = -1;
          q_n[mb][mc][2] = -1;
          q_n[mb][mc][3] = -1;
        }
      }
      if (!match)
      {
        m = c+total;
        hash[n0]->Add_To_List(m);
        hash[n1]->Add_To_List(m);
        hash[n2]->Add_To_List(m);
        hash[n3]->Add_To_List(m);
      }
    }
    total += nq[b];
  }

  // process polygons
  for (n=0; n < nn; n++)
    hash[n]->Redimension(0);

  for (b=0; b < nb; b++)
  {
    for (c=0; c < ngon[b]; c++)
    {
      n0 = ngon_n[b][c][1];
      match = false;
      for (i=0; i < hash[n0]->max && !match; i++)
      {
        m = hash[n0]->list[i];
        match = true;
        for (j=2; j <= ngon_n[b][c][0] && match; j++)
        {
          n1 = ngon_n[b][c][j];
          if (!hash[n1]->Is_In_List(m))
            match = false;
        }
        if (m != (c+total) && match)
        {
          for (j=1; j <= ngon_n[b][c][0]; j++)
          {
            n1 = ngon_n[b][c][j];
            hash[n1]->Delete_From_List(m);
            ngon_n[b][c][j] = -1;
          }
          for (k=mb=0; mb < nb; mb++)
          {
            if (k+ngon[mb] > m)
              break;
            k+=ngon[mb];
          }
          mc = m-k;
          for (j=1; j <= ngon_n[mb][mc][0]; j++)
            ngon_n[mb][mc][j] = -1;
        }
      }
      if (!match)
      {
        m = c+total;
        for (j=1; j <= ngon_n[b][c][0]; j++)
        {
          n1 = ngon_n[b][c][j];
          hash[n1]->Add_To_List(m);
        }
      }
    }
    total += ngon[b];
  }

  for (n=0; n < nn; n++)
    delete hash[n];
  delete[] hash;

  // now count elements per boundary and compress if necessary
  for (b=0; b < nb; b++)
  {
    m = 0;
    for (c=0; c < nt[b]; c++)
    {
      if (t_n[b][c][0] >= 0)
      {
        if (m != c)
        {
          t_n[b][m][0] = t_n[b][c][0];
          t_n[b][m][1] = t_n[b][c][1];
          t_n[b][m][2] = t_n[b][c][2];
        }
        m++;
      }
    }
    if (m < nt[b])
    {
      for (c=m; c < nt[b]; c++)
        free(t_n[b][c]);
      if (m > 0)
        t_n[b] = (int**)realloc((void*)t_n[b],m*sizeof(int*));
      else
      {
        free(t_n[b]);
        t_n[b] = 0;
      }
      nt[b] = m;
    }

    m = 0;
    for (c=0; c < nq[b]; c++)
    {
      if (q_n[b][c][0] >= 0)
      {
        if (m != c)
        {
          q_n[b][m][0] = q_n[b][c][0];
          q_n[b][m][1] = q_n[b][c][1];
          q_n[b][m][2] = q_n[b][c][2];
          q_n[b][m][3] = q_n[b][c][3];
        }
        m++;
      }
    }
    if (m < nq[b])
    {
      for (c=m; c < nq[b]; c++)
        free(q_n[b][c]);
      if (m > 0)
        q_n[b] = (int**)realloc((void*)q_n[b],m*sizeof(int*));
      else
      {
        free(q_n[b]);
        q_n[b] = 0;
      }
      nq[b] = m;
    }

    m = 0;
    for (c=0; c < ngon[b]; c++)
    {
      if (ngon_n[b][c][1] >= 0)
      {
        if (m != c)
        {
          ngon_n[b][m] = (int*)realloc((void*)ngon_n[b][m],(ngon_n[b][c][0]+1)*sizeof(int));
          for (i=1; i <= ngon_n[b][c][0]; i++)
            ngon_n[b][m][i] = ngon_n[b][c][i];
        }
        m++;
      }
    }
    if (m < ngon[b])
    {
      for (c=m; c < ngon[b]; c++)
        free(ngon_n[b][c]);
      if (m > 0)
        ngon_n[b] = (int**)realloc((void*)ngon_n[b],m*sizeof(int*));
      else
      {
        free(ngon_n[b]);
        ngon_n[b] = 0;
      }
      ngon[b] = m;
    }
  }

  // check for boundaries with zero elements
  m = 0;
  for (b=0; b < nb; b++)
  {
    if (nt[b]+nq[b]+ngon[b] > 0)
    {
      if (m != b)
      {
        if (nt[b] > 0)
        {
          if (nt[m] > 0)
            for (i=0; i < nt[m]; i++)
              free(t_n[m][i]);
          t_n[m] = (int**)realloc((void*)t_n[m],nt[b]*sizeof(int*));
          for (i=0; i < nt[b]; i++)
          {
            t_n[m][i] = (int*)malloc(3*sizeof(int));
            for (j=0; j < 3; j++)
              t_n[m][i][j] = t_n[b][i][j];
            free(t_n[b][i]);
          }
          free(t_n[b]); t_n[b]=0;
          nt[m] = nt[b];
          nt[b] = 0;
        }

        if (nq[b] > 0)
        {
          if (nq[m] > 0)
            for (i=0; i < nq[m]; i++)
              free(q_n[m][i]);
          q_n[m] = (int**)realloc((void*)q_n[m],nq[b]*sizeof(int*));
          for (i=0; i < nq[b]; i++)
          {
            q_n[m][i] = (int*)malloc(4*sizeof(int));
            for (j=0; j < 4; j++)
              q_n[m][i][j] = q_n[b][i][j];
            free(q_n[b][i]);
          }
          free(q_n[b]); q_n[b]=0;
          nq[m] = nq[b];
          nq[b] = 0;
        }

        if (ngon[b] > 0)
        {
          if (ngon[m] > 0)
            for (i=0; i < ngon[m]; i++)
              free(ngon_n[m][i]);
          ngon_n[m] = (int**)realloc((void*)ngon_n[m],ngon[b]*sizeof(int*));
          for (i=0; i < ngon[b]; i++)
          {
            ngon_n[m][i] = (int*)malloc((ngon_n[b][i][0]+1)*sizeof(int));
            ngon_n[m][i][0] = ngon_n[b][i][0];
            for (j=1; j <= ngon_n[b][i][0]; j++)
              ngon_n[m][i][j] = ngon_n[b][i][j];
            free(ngon_n[b][i]);
          }
          free(ngon_n[b]); ngon_n[b] = 0;
          ngon[m] = ngon[b];
          ngon[b] = 0;
        }

        sprintf(b_name[m],"%-32s",b_name[b]);
      }

      printf("\nBoundary %i = %s",m,b_name[m]);
      printf("\n  # of triangles      = %i",nt[m]);
      printf("\n  # of quadrilaterals = %i",nq[m]);
      printf("\n  # of polygons       = %i",ngon[m]);

      m++;
    } else
      printf("\nDeleting boundary %i, %s",b,b_name[b]);
  }
  if (m < nb)
  {
    for (b=m; b < nb; b++)
    {
      if (nt[b] > 0)
      {
        for (c=0; c < nt[b]; c++)
          if (t_n[b][c] != 0) free(t_n[b][c]);
        if (t_n[b] != 0) free(t_n[b]);
        nt[b] = 0;
      }
      if (nq[b] > 0)
      {
        for (c=0; c < nq[b]; c++)
          if (q_n[b][c] != 0) free(q_n[b][c]);
        if (q_n[b] != 0) free(q_n[b]);
        nq[b] = 0;
      }
      if (ngon[b] > 0)
      {
        for (c=0; c < ngon[b]; c++)
          if (ngon_n[b][c] != 0) free(ngon_n[b][c]);
        if (ngon_n[b] != 0) free(ngon_n[b]);
        ngon[b] = 0;
      }
    }
    nb = m;
    nt = (int*)realloc((void*)nt,nb*sizeof(int));
    nq = (int*)realloc((void*)nq,nb*sizeof(int));
    ngon = (int*)realloc((void*)ngon,nb*sizeof(int));
    t_n = (int***)realloc((void*)t_n,nb*sizeof(int**));
    q_n = (int***)realloc((void*)q_n,nb*sizeof(int**));
    ngon_n = (int***)realloc((void*)ngon_n,nb*sizeof(int**));
  }

  x = (double*)malloc(nn*sizeof(double));
  y = (double*)malloc(nn*sizeof(double));
  z = (double*)malloc(nn*sizeof(double));
  for (n=0; n < nn; n++)
  {
    x[n] = node[n][0];
    y[n] = node[n][1];
    z[n] = node[n][2];
  }

  //Write_Mesh(mname,b_name,nn,x,y,z,nb,nt,nq,ngon,t_n,q_n,ngon_n,ntet,tet_n,npyr,pyr_n,npri,pri_n,nhex,hex_n,nply,poly_n);
  Write_Mesh(mname,b_name,nn,x,y,z,nb,nt,nq,ngon,t_n,q_n,ngon_n,ntet,tet_n,npyr,pyr_n,npri,pri_n,nhex,hex_n,nply,poly_n,
            nvc,vc_name,tet_vc,pyr_vc,pri_vc,hex_vc,poly_vc);

  free(x);
  free(y);
  free(z);

  // free up memory
  for (b=0; b < nb; b++)
  {
    free(b_name[b]);
    for (i=0; i < nt[b]; i++)
      if (t_n[b][i] != 0) free(t_n[b][i]);
    for (i=0; i < nq[b]; i++)
      if (q_n[b][i] != 0) free(q_n[b][i]);
    for (i=0; i < ngon[b]; i++)
      if (ngon_n[b][i] != 0) free(ngon_n[b][i]);
    if (t_n[b] != 0) free(t_n[b]);
    if (q_n[b] != 0) free(q_n[b]);
    if (ngon_n[b] != 0) free(ngon_n[b]);
  }
  if (t_n != 0) free(t_n);
  if (q_n != 0) free(q_n);
  if (ngon_n != 0) free(ngon_n);
  if (nt != 0) free(nt);
  if (nq != 0) free(nq);
  if (ngon != 0) free(ngon);
  for (i=0; i < ntet; i++)
    if (tet_n[i] != 0) free(tet_n[i]);
  if (tet_n != 0) free(tet_n);
  if (tet_vc != 0) free(tet_vc);
  for (i=0; i < npyr; i++)
    if (pyr_n[i] != 0) free(pyr_n[i]);
  if (pyr_n != 0) free(pyr_n);
  if (pyr_vc != 0) free(pyr_vc);
  for (i=0; i < npri; i++)
    if (pri_n[i] != 0) free(pri_n[i]);
  if (pri_n != 0) free(pri_n);
  if (pri_vc != 0) free(pri_vc);
  for (i=0; i < nhex; i++)
    if (hex_n[i] != 0) free(hex_n[i]);
  if (hex_n != 0) free(hex_n);
  if (hex_vc != 0) free(hex_vc);
  for (i=0; i < nply; i++)
  {
    for (j=poly_n[i][0][0]; j >= 0; j--)
      if (poly_n[i][j] != 0) free(poly_n[i][j]);
    if (poly_n[i] != 0) free(poly_n[i]);
  }
  if (poly_n != 0) free(poly_n);
  if (poly_vc != 0) free(poly_vc);
  if (node != 0) free(node);

  return(error);
}
