#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "Point.h"
#include "Vector.h"
#include "Util.h"
#include "Linked_List.h"
#include "List.h"
#include "trimesh.h"
#include "convert.h"
#include "mesh_io.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))

int convert(char cname[], char sname[], int wind, int volchk, int pflag, char vname[])
{
  int b, c, i, j, k, m, n, n0, n1, n2, n3, n4, n5, n6, n7;
  int nn, nb, ntet, npyr, npri, nhex, nply, nvol;
  int *nt, *nq, *ngon, ***t_n, ***q_n, ***ngon_n;
  int **tet_n, **pyr_n, **pri_n, **hex_n, ***poly_n;
  int *tet_vc, *pyr_vc, *pri_vc, *hex_vc, *poly_vc;
  double *x, *y, *z;
  Point *node;
  char **b_name, **vol_name;
  FILE *fp;
  int error = 0;
 
  // set all counts to 0
  nvol = nn = nb = ntet = npyr = npri = nhex = nply = 0;
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
  x = 0;
  y = 0;
  z = 0;
  node = 0;
  tet_vc = 0;
  pyr_vc = 0;
  pri_vc = 0;
  hex_vc = 0;
  poly_vc = 0;

  Read_Mesh(cname,b_name,nn,x,y,z,nb,nt,nq,ngon,t_n,q_n,ngon_n,ntet,tet_n,npyr,pyr_n,npri,pri_n,nhex,hex_n,nply,poly_n,
            nvol,vol_name,tet_vc,pyr_vc,pri_vc,hex_vc,poly_vc);

  node = (Point*)malloc(nn*sizeof(Point));
  for (n=0; n < nn; n++)
    node[n] = Point(x[n],y[n],z[n]);

  free(x);
  free(y);
  free(z);

  if (ntet > 0 && tet_vc==0)
  {
    tet_vc = (int*)malloc(ntet*sizeof(int));
    for (c=0; c < ntet; c++)
      tet_vc[c] = 0;
  }
  if (npyr > 0 && pyr_vc==0)
  {
    pyr_vc = (int*)malloc(npyr*sizeof(int));
    for (c=0; c < npyr; c++)
      pyr_vc[c] = 0;
  }
  if (npri > 0 && pri_vc==0)
  {
    pri_vc = (int*)malloc(npri*sizeof(int));
    for (c=0; c < npri; c++)
      pri_vc[c] = 0;
  }
  if (nhex > 0 && hex_vc==0)
  {
    hex_vc = (int*)malloc(nhex*sizeof(int));
    for (c=0; c < nhex; c++)
      hex_vc[c] = 0;
  }
  if (nply > 0 && poly_vc==0)
  {
    poly_vc = (int*)malloc(nply*sizeof(int));
    for (c=0; c < nply; c++)
      poly_vc[c] = 0;
  }

  // check element windings
  if (volchk)
  {
    printf("\nVOLUME CHECK: ONLY CHECKING THE TETRAHEDRAL ELEMENTS!!!");
    fflush(stdout);

    Point p0, p1, p2, p3;
    m = 0;
    for (c=0; c < ntet; c++)
    {
      p0 = node[n0 = tet_n[c][0]];
      p1 = node[n1 = tet_n[c][1]];
      p2 = node[n2 = tet_n[c][2]];
      p3 = node[n3 = tet_n[c][3]];
      if (tetrahedral_volume(p0, p1, p2, p3) < 0.0)
      {
        tet_n[c][1] = n2;
        tet_n[c][2] = n1;
        m++;
      }
    }
    if (m > 0)
      printf("\nNumber of tetrahedral windings corrected = %d",m);
    fflush(stdout);
  }

  // check windings of boundary facets
  if (wind)
  {
    Linked_List **tet_hash, **pyr_hash, **pri_hash, **hex_hash, **poly_hash;
    Linked_Node *hd0, *hd1, *hd2, *hd3;
    Point p0, p1, p2, p3, fcg, cg;
    Vector v1, v2, norm;
    int *tag;
    int  f;

    tag = new int[nn];
    for (n=0; n < nn; n++)
      tag[n] = 0;

    for (b=0; b < nb; b++)
    {
      for (f=0; f < nt[b]; f++)
      {
        for (i=0; i < 3; i++)
        {
          n = t_n[b][f][i];
          tag[n] = 1;
        }
      }
      for (f=0; f < nq[b]; f++)
      {
        for (i=0; i < 4; i++)
        {
          n = q_n[b][f][i];
          tag[n] = 1;
        }
      }
      for (f=0; f < ngon[b]; f++)
      {
        for (i=1; i <= ngon_n[b][f][0]; i++)
        {
          n = ngon_n[b][f][i];
          tag[n] = 1;
        }
      }
    }

    printf("\nBoundary nodes tagged.");
    fflush(stdout);

    tet_hash = (Linked_List**)malloc(nn*sizeof(Linked_List*));
    pyr_hash = (Linked_List**)malloc(nn*sizeof(Linked_List*));
    pri_hash = (Linked_List**)malloc(nn*sizeof(Linked_List*));
    hex_hash = (Linked_List**)malloc(nn*sizeof(Linked_List*));
    poly_hash = (Linked_List**)malloc(nn*sizeof(Linked_List*));

    for (n=0; n < nn; n++)
    {
      tet_hash[n] = new Linked_List();
      pyr_hash[n] = new Linked_List();
      pri_hash[n] = new Linked_List();
      hex_hash[n] = new Linked_List();
      poly_hash[n] = new Linked_List();
    }

    for (c=0; c < ntet; c++)
    {
      for (i=0; i < 4; i++)
      {
        n = tet_n[c][i];
        if (tag[n])
          tet_hash[n]->Insert(c);
      }
    }
    for (c=0; c < npyr; c++)
    {
      for (i=0; i < 5; i++)
      {
        n = pyr_n[c][i];
        if (tag[n])
          pyr_hash[n]->Insert(c);
      }
    }
    for (c=0; c < npri; c++)
    {
      for (i=0; i < 6; i++)
      {
        n = pri_n[c][i];
        if (tag[n])
          pri_hash[n]->Insert(c);
      }
    }
    for (c=0; c < nhex; c++)
    {
      for (i=0; i < 8; i++)
      {
        n = hex_n[c][i];
        if (tag[n])
          hex_hash[n]->Insert(c);
      }
    }
    for (c=0; c < nply; c++)
    {
      for (i=1; i <= poly_n[c][0][0]; i++)
      {
        for (j=1; j <= poly_n[c][i][0]; j++)
        {
          n = poly_n[c][i][j];
          if (tag[n] && !poly_hash[n]->In_list(c))
            poly_hash[n]->Insert(c);
        }
      }
    }

    delete[] tag;

    printf("\nHash tables created.");
    fflush(stdout);

    for (b=0; b < nb; b++)
    {
      printf("\nChecking windings for boundary %d",b+1);
      fflush(stdout);
  
      int tflip, qflip, gflip;
      tflip = 0;
      qflip = 0;
      gflip = 0;
      for (f=0; f < nt[b]; f++)
      {
        n0 = t_n[b][f][0];
        n1 = t_n[b][f][1];
        n2 = t_n[b][f][2];
        p0 = node[n0];
        p1 = node[n1];
        p2 = node[n2];
        fcg = (p0+p1+p2)/3.0;
        v1 = Vector(p0,p1);
        v2 = Vector(p0,p2);
        norm = v1 % v2;
        norm.normalize();

        c = -1;
        if (c < 0)
        {
          // check tets
          hd0 = tet_hash[n0]->head;
          while (hd0 && c < 0)
          {
            n = hd0->data;
            hd1 = tet_hash[n1]->head;
            while (hd1 && c < 0)
            {
              if (hd1->data == n)
              {
                hd2 = tet_hash[n2]->head;
                while (hd2 && c < 0)
                {
                  if (hd2->data == n)
                    c = n;
                  hd2 = hd2->next;
                }
              }
              hd1 = hd1->next;
            }
            hd0 = hd0->next;
          }
          if (c >= 0)
            cg = (node[tet_n[c][0]]+node[tet_n[c][1]]+
                  node[tet_n[c][2]]+node[tet_n[c][3]])/4.0;
        }
        if (c < 0)
        {
          // check pyramid
          hd0 = pyr_hash[n0]->head;
          while (hd0 && c < 0)
          {
            n = hd0->data;
            hd1 = pyr_hash[n1]->head;
            while (hd1 && c < 0)
            {
              if (hd1->data == n)
              {
                hd2 = pyr_hash[n2]->head;
                while (hd2 && c < 0)
                {
                  if (hd2->data == n)
                    c = n;
                  hd2 = hd2->next;
                }
              }
              hd1 = hd1->next;
            }
            hd0 = hd0->next;
          }
          if (c >= 0)
            cg = (node[pyr_n[c][0]]+node[pyr_n[c][1]]+
                  node[pyr_n[c][2]]+node[pyr_n[c][3]]+
                  node[pyr_n[c][4]])/5.0;
        }
        if (c < 0)
        {
          // check prism
          hd0 = pri_hash[n0]->head;
          while (hd0 && c < 0)
          {
            n = hd0->data;
            hd1 = pri_hash[n1]->head;
            while (hd1 && c < 0)
            {
              if (hd1->data == n)
              {
                hd2 = pri_hash[n2]->head;
                while (hd2 && c < 0)
                {
                  if (hd2->data == n)
                    c = n;
                  hd2 = hd2->next;
                }
              }
              hd1 = hd1->next;
            }
            hd0 = hd0->next;
          }
          if (c >= 0)
            cg = (node[pri_n[c][0]]+node[pri_n[c][1]]+
                  node[pri_n[c][2]]+node[pri_n[c][3]]+
                  node[pri_n[c][4]]+node[pri_n[c][5]])/6.0;
        }
        if (c < 0)
        {
          // check polyhedra
          hd0 = poly_hash[n0]->head;
          while (hd0 && c < 0)
          {
            n = hd0->data;
            hd1 = poly_hash[n1]->head;
            while (hd1 && c < 0)
            {
              if (hd1->data == n)
              {
                hd2 = poly_hash[n2]->head;
                while (hd2 && c < 0)
                {
                  if (hd2->data == n)
                    c = n;
                  hd2 = hd2->next;
                }
              }
              hd1 = hd1->next;
            }
            hd0 = hd0->next;
          }
          if (c >= 0)
          {
            cg = Point(0.0,0.0,0.0);
            k=0;
            for (i=1; i <= poly_n[c][0][0]; i++)
            {
              for (j=1; j <= poly_n[c][i][0]; j++)
              {
                cg += node[poly_n[c][i][j]];
                k++;
              }
            }
            cg /= MAX(1,k);
          }
        }
        if (c >= 0)
        {
          v1 = Vector(fcg,cg);
          if (v1 * norm < 0.0)
          {
            t_n[b][f][1] = n2;
            t_n[b][f][2] = n1;
            tflip++;
          }
        } else if (c < 0)
        {
          printf("\nBoundary %d:",b+1);
          printf("\nNo cell match for boundary nodes %d, %d, %d!",n0,n1,n2);
          exit(0);
        }
      }
      for (f=0; f < nq[b]; f++)
      {
        n0 = q_n[b][f][0];
        n1 = q_n[b][f][1];
        n2 = q_n[b][f][2];
        n3 = q_n[b][f][3];
        p0 = node[n0];
        p1 = node[n1];
        p2 = node[n2];
        p3 = node[n3];
        fcg = (p0+p1+p2+p3)/4.0;
        v1 = Vector(p0,p2);
        v2 = Vector(p1,p3);
        norm = v1 % v2;
        norm.normalize();
        c = -1;
        if (c < 0)
        {
          // check pyramid
          hd0 = pyr_hash[n0]->head;
          while (hd0 && c < 0)
          {
            n = hd0->data;
            hd1 = pyr_hash[n1]->head;
            while (hd1 && c < 0)
            {
              if (hd1->data == n)
              {
                hd2 = pyr_hash[n2]->head;
                while (hd2 && c < 0)
                {
                  if (hd2->data == n)
                  {
                    hd3 = pyr_hash[n3]->head;
                    while (hd3 && c < 0)
                    {
                      if (hd3->data == n)
                        c = n;
                      hd3 = hd3->next;
                    }
                  }
                  hd2 = hd2->next;
                }
              }
              hd1 = hd1->next;
            }
            hd0 = hd0->next;
          }
          if (c >= 0)
            cg = (node[pyr_n[c][0]]+node[pyr_n[c][1]]+
                  node[pyr_n[c][2]]+node[pyr_n[c][3]]+
                  node[pyr_n[c][4]])/5.0;
        }
        if (c < 0)
        {
          // check prism
          hd0 = pri_hash[n0]->head;
          while (hd0 && c < 0)
          {
            n = hd0->data;
            hd1 = pri_hash[n1]->head;
            while (hd1 && c < 0)
            {
              if (hd1->data == n)
              {
                hd2 = pri_hash[n2]->head;
                while (hd2 && c < 0)
                {
                  if (hd2->data == n)
                  {
                    hd3 = pri_hash[n3]->head;
                    while (hd3 && c < 0)
                    {
                      if (hd3->data == n)
                        c = n;
                      hd3 = hd3->next;
                    }
                  }
                  hd2 = hd2->next;
                }
              }
              hd1 = hd1->next;
            }
            hd0 = hd0->next;
          }
          if (c >= 0)
            cg = (node[pri_n[c][0]]+node[pri_n[c][1]]+
                  node[pri_n[c][2]]+node[pri_n[c][3]]+
                  node[pri_n[c][4]]+node[pri_n[c][5]])/6.0;
        }
        if (c < 0)
        {
          // check hexahedra
          hd0 = hex_hash[n0]->head;
          while (hd0 && c < 0)
          {
            n = hd0->data;
            hd1 = hex_hash[n1]->head;
            while (hd1 && c < 0)
            {
              if (hd1->data == n)
              {
                hd2 = hex_hash[n2]->head;
                while (hd2 && c < 0)
                {
                  if (hd2->data == n)
                  {
                    hd3 = hex_hash[n3]->head;
                    while (hd3 && c < 0)
                    {
                      if (hd3->data == n)
                        c = n;
                      hd3 = hd3->next;
                    }
                  }
                  hd2 = hd2->next;
                }
              }
              hd1 = hd1->next;
            }
            hd0 = hd0->next;
          }
          if (c >= 0)
            cg = (node[hex_n[c][0]]+node[hex_n[c][1]]+
                  node[hex_n[c][2]]+node[hex_n[c][3]]+
                  node[hex_n[c][4]]+node[hex_n[c][5]]+
                  node[hex_n[c][6]]+node[hex_n[c][7]])/8.0;
        }
        if (c < 0)
        {
          // check polyhedra
          hd0 = poly_hash[n0]->head;
          while (hd0 && c < 0)
          {
            n = hd0->data;
            hd1 = poly_hash[n1]->head;
            while (hd1 && c < 0)
            {
              if (hd1->data == n)
              {
                hd2 = poly_hash[n2]->head;
                while (hd2 && c < 0)
                {
                  if (hd2->data == n)
                  {
                    hd3 = poly_hash[n3]->head;
                    while (hd3 && c < 0)
                    {
                      if (hd3->data == n)
                        c = n;
                      hd3 = hd3->next;
                    }
                  }
                  hd2 = hd2->next;
                }
              }
              hd1 = hd1->next;
            }
            hd0 = hd0->next;
          }
          if (c >= 0)
          {
            cg = Point(0.0,0.0,0.0);
            k=0;
            for (i=1; i <= poly_n[c][0][0]; i++)
            {
              for (j=1; j <= poly_n[c][i][0]; j++)
              {
                cg += node[poly_n[c][i][j]];
                k++;
              }
            }
            cg /= MAX(1,k);
          }
        }
        if (c >= 0)
        {
          v1 = Vector(fcg,cg);
          if (v1 * norm < 0.0)
          {
            q_n[b][f][1] = n3;
            q_n[b][f][3] = n1;
            qflip++;
          }
        } else if (c < 0)
        {
          printf("\nBoundary %d:",b+1);
          printf("\nNo cell match for boundary nodes %d, %d, %d, %d!",n0,n1,n2,n3);
          exit(0);
        }
      }
      for (f=0; f < ngon[b]; f++)
      {
        norm = Vector(0.0,0.0,0.0);
        n0 = ngon_n[b][f][1];
        p0 = node[n0];
        fcg = p0*2.0;
        for (i=2; i < ngon_n[b][f][0]; i++)
        {
          n1 = ngon_n[b][f][i];
          n2 = ngon_n[b][f][i+1];
          p1 = node[n1];
          p2 = node[n2];
          fcg += p1 + p2;
          v1 = Vector(p0,p1);
          v2 = Vector(p0,p2);
          norm += v1 % v2;
        }
        fcg /= (ngon_n[b][f][0]*2);
        norm.normalize();
        c = -1;
        if (c < 0 && ngon_n[b][f][0] == 3)
        {
          // check tets
          hd0 = tet_hash[n0]->head;
          while (hd0 && c < 0)
          {
            n = hd0->data;
            int match = 1;
            for (i=2; i <= ngon_n[b][f][0] && match; i++)
            {
              n1 = ngon_n[b][f][i];
              match = tet_hash[n1]->In_list(n);
            }
            if (match)
              c = n;

            hd0 = hd0->next;
          }
          if (c >= 0)
          {
            cg = (node[tet_n[c][0]]+node[tet_n[c][1]]+
                  node[tet_n[c][2]]+node[tet_n[c][3]])/4.0;
          }
        }
        if (c < 0 && (ngon_n[b][f][0] == 3 || ngon_n[b][f][0] == 4))
        {
          // check pyramids
          hd0 = pyr_hash[n0]->head;
          while (hd0 && c < 0)
          {
            n = hd0->data;
            int match = 1;
            for (i=2; i <= ngon_n[b][f][0] && match; i++)
            {
              n1 = ngon_n[b][f][i];
              match = pyr_hash[n1]->In_list(n);
            }
            if (match)
              c = n;

            hd0 = hd0->next;
          }
          if (c >= 0)
          {
            cg = (node[pyr_n[c][0]]+node[pyr_n[c][1]]+
                  node[pyr_n[c][2]]+node[pyr_n[c][3]]+
                  node[pyr_n[c][4]])/5.0;
          }
        }
        if (c < 0 && (ngon_n[b][f][0] == 3 || ngon_n[b][f][0] == 4))
        {
          // check prisms
          hd0 = pri_hash[n0]->head;
          while (hd0 && c < 0)
          {
            n = hd0->data;
            int match = 1;
            for (i=2; i <= ngon_n[b][f][0] && match; i++)
            {
              n1 = ngon_n[b][f][i];
              match = pri_hash[n1]->In_list(n);
            }
            if (match)
              c = n;

            hd0 = hd0->next;
          }
          if (c >= 0)
          {
            cg = (node[pri_n[c][0]]+node[pri_n[c][1]]+
                  node[pri_n[c][2]]+node[pri_n[c][3]]+
                  node[pri_n[c][4]]+node[pri_n[c][5]])/6.0;
          }
        }
        if (c < 0 && ngon_n[b][f][0] == 4)
        {
          // check hexes
          hd0 = hex_hash[n0]->head;
          while (hd0 && c < 0)
          {
            n = hd0->data;
            int match = 1;
            for (i=2; i <= ngon_n[b][f][0] && match; i++)
            {
              n1 = ngon_n[b][f][i];
              match = hex_hash[n1]->In_list(n);
            }
            if (match)
              c = n;

            hd0 = hd0->next;
          }
          if (c >= 0)
          {
            cg = (node[hex_n[c][0]]+node[hex_n[c][1]]+
                  node[hex_n[c][2]]+node[hex_n[c][3]]+
                  node[hex_n[c][4]]+node[hex_n[c][5]]+
                  node[hex_n[c][6]]+node[hex_n[c][7]])/8.0;
          }
        }
        if (c < 0)
        {
          // check polyhedra
          hd0 = poly_hash[n0]->head;
          while (hd0 && c < 0)
          {
            n = hd0->data;
            int match = 1;
            for (i=2; i <= ngon_n[b][f][0] && match; i++)
            {
              n1 = ngon_n[b][f][i];
              match = poly_hash[n1]->In_list(n);
            }
            if (match)
              c = n;

            hd0 = hd0->next;
          }
          if (c >= 0)
          {
            cg = Point(0.0,0.0,0.0);
            k=0;
            for (i=1; i <= poly_n[c][0][0]; i++)
            {
              for (j=1; j <= poly_n[c][i][0]; j++)
              {
                cg += node[poly_n[c][i][j]];
                k++;
              }
            }
            cg /= MAX(1,k);
          }
        }
        if (c >= 0)
        {
          v1 = Vector(fcg,cg);
          if (v1 * norm < 0.0)
          {
            int *tmp = new int[ngon_n[b][f][0]];
            j=ngon_n[b][f][0]-1;
            for (i=1; i <= ngon_n[b][f][0]; i++)
              tmp[j--] = ngon_n[b][f][i];
            for (i=1; i <= ngon_n[b][f][0]; i++)
              ngon_n[b][f][i] = tmp[i-1];
            delete[] tmp;
            gflip++;
          }
        } else if (c < 0)
        {
          printf("\nBoundary %d:",b+1);
          printf("\nNo cell match for boundary polygon %d!",f);
          exit(0);
        }
      }
      if (tflip > 0)
        printf("\nNumber of triangle windings corrected      = %d",tflip);
      if (qflip > 0)
        printf("\nNumber of quadrilateral windings corrected = %d",qflip);
      if (gflip > 0)
        printf("\nNumber of polygonal windings corrected     = %d",gflip);
      fflush(stdout);
    }
    for (n=0; n < nn; n++)
    {
      delete tet_hash[n];
      delete pyr_hash[n];
      delete pri_hash[n];
      delete hex_hash[n];
      delete poly_hash[n];
    }
    free(tet_hash);
    free(pyr_hash);
    free(pri_hash);
    free(hex_hash);
    free(poly_hash);
  }

  // convert polyhedra to standard by creating centroid
  if (nply > 0 && pflag)
  {
    int f;
    double dot, dsmin;

    // ensure pointers are initialized
    for (b=0; b < nb; b++)
    {
      if (nt[b] == 0) t_n[b] = 0;
      if (nq[b] == 0) q_n[b] = 0;
      if (ngon[b] == 0) ngon_n[b] = 0;
    }

    k=0;
    for (c=0; c < nply && !k; c++)
      for (i=1; i <= poly_n[c][0][0] && !k; i++)
        if (poly_n[c][i][0] > 4)
          k=1;

    if (k)
    {
      List **p_hash;
      List tnode;
      int *nbs, ***bs;
      double *x, *y;
      Vector dr, v1, v2, norm;

      // create node to polyhedral hash table
      p_hash = (List**)malloc(nn*sizeof(List*));
      for (n=0; n < nn; n++)
        p_hash[n] = new List();
      for (c=0; c < nply; c++)
      {
        for (i=1; i <= poly_n[c][0][0]; i++)
        {
          for (j=1; j <= poly_n[c][i][0]; j++)
          {
            n = poly_n[c][i][j];
            p_hash[n]->Check_List(c);
          }
        }
      }

      // triangulate boundary polygons
      for (b=0; b < nb; b++)
      {
        for (f=0; f < ngon[b]; f++)
        {
          if (ngon_n[b][f][0] > 4)
          {
            // convert to triangles
            tnode.Redimension(0);
            Point cg = Point(0.0,0.0,0.0);
            for (i=1; i <= ngon_n[b][f][0]; i++)
            {
              n0 = ngon_n[b][f][i];
              tnode.Check_List(n0);
              cg += node[n0];
            }
            cg /= ngon_n[b][f][0];
            norm = Vector(0.0,0.0,0.0);
            nbs = new int[1];
            bs = new int**[1];
            bs[0] = new int*[ngon_n[b][f][0]];
            nbs[0] = 0;
            for (i=1; i <= ngon_n[b][f][0]; i++)
            {
              if (i==1)
                n0 = ngon_n[b][f][ngon_n[b][f][0]];
              else
                n0 = ngon_n[b][f][i-1];
              n1 = ngon_n[b][f][i];
              bs[0][nbs[0]] = new int[2];
              bs[0][nbs[0]][0] = tnode.Index(n0);
              bs[0][nbs[0]][1] = tnode.Index(n1);
              nbs[0]++;
              v1 = Vector(cg,node[n0]);
              v2 = Vector(cg,node[n1]);
              norm += v1 % v2;
            }
            norm.normalize();
            x = new double[tnode.max];
            y = new double[tnode.max];
            n0 = ngon_n[b][f][1];
            v1 = Vector(cg,node[n0]);
            dot = v1 * norm;
            v1 -= norm*dot;
            v1.normalize();
            v2 = norm % v1;
            v2.normalize();
            dsmin = 1.0e20;
            for (i=1; i <= ngon_n[b][f][0]; i++)
            {
              if (i==1)
                n0 = ngon_n[b][f][ngon_n[b][f][0]];
              else
                n0 = ngon_n[b][f][i-1];
              n1 = ngon_n[b][f][i];
              dr = Vector(cg,node[n0]);
              int i0 = tnode.Index(n0);
              x[i0] = dr*v1;
              y[i0] = dr*v2;
              dr = Vector(cg,node[n1]);
              int i1 = tnode.Index(n1);
              x[i1] = dr*v1;
              y[i1] = dr*v2;
              double dx = x[i1]-x[i0];
              double dy = y[i1]-y[i0];
              double mag = sqrt(dx*dx+dy*dy);
              dsmin = MIN(dsmin,mag);
            }
            if (dsmin < 1.0e-15)
            {
              printf("\nCONVERT: coincident points detected in triangulation process!");
              fflush(stdout);
              exit(0);
            }
            int tdim = tnode.max*5;
            int (*tri)[3] = new int[tdim][3];
            int lt, t;

            if ((lt = trimesh(tnode.max,tdim,1,0,nbs,bs,x,y,tri)) > 0)
            {
              // store in boundary triangle list
              t_n[b] = (int**)realloc((void*)t_n[b],(nt[b]+lt)*sizeof(int*));
              for (i=nt[b]; i < nt[b]+lt; i++)
                t_n[b][i] = (int*)malloc(3*sizeof(int));
              for (i=0; i < lt; i++)
              {
                t_n[b][nt[b]][0] = tnode.list[tri[i][0]];
                t_n[b][nt[b]][1] = tnode.list[tri[i][1]];
                t_n[b][nt[b]][2] = tnode.list[tri[i][2]];
                nt[b]++;
              }

              // find polyhedron
              c = -1;
              n0 = ngon_n[b][f][1];
              for (i=0; i < p_hash[n0]->max && c < 0; i++)
              {
                n = p_hash[n0]->list[i];
                int match = 1;
                for (j=2; j <= ngon_n[b][f][0] && match; j++)
                {
                  n1 = ngon_n[b][f][j];
                  match = p_hash[n1]->Is_In_List(n);
                }
                if (match)
                  c = n;
              }
              if (c < 0)
              {
                printf("\nCONVERT: Polyhedron for boundary polygon not found!");
                fflush(stdout);
                exit(0);
              }
              // find polygon face in polyhedron
              int bf = -1;
              for (i=1; i <= poly_n[c][0][0] && bf < 0; i++)
              {
                if (poly_n[c][i][0] != ngon_n[b][f][0])
                  continue;
                int *tag = new int[poly_n[c][i][0]];
                for (j=0; j < poly_n[c][i][0]; j++)
                  tag[j] = 0;
                for (j=1; j <= ngon_n[b][f][0]; j++)
                  for (k=1; k <= poly_n[c][i][0]; k++)
                    if (ngon_n[b][f][j] == poly_n[c][i][k])
                      tag[k-1] = 1;
                k=0;
                for (j=0; j < poly_n[c][i][0]; j++)
                  k+=tag[j];
                if (k == poly_n[c][i][0])
                  bf = i;
                delete[] tag;
              }
              if (bf < 0)
              {
                printf("\nCONVERT: Polyhedron face for boundary polygon not found!");
                fflush(stdout);
                exit(0);
              }
              // replace with last face
              int last = poly_n[c][0][0];
              if (bf < last)
              {
                poly_n[c][bf] = (int*)realloc((void*)poly_n[c][bf],(poly_n[c][last][0]+1)*sizeof(int));
                poly_n[c][bf][0] = poly_n[c][last][0];
                for (j=1; j <= poly_n[c][last][0]; j++)
                  poly_n[c][bf][j] = poly_n[c][last][j];
              }
              free(poly_n[c][last]);
              poly_n[c][0][0]--;
              poly_n[c] = (int**)realloc((void*)poly_n[c],(poly_n[c][0][0]+1)*sizeof(int*));
              
              int next = poly_n[c][0][0]+1;
              // now add triangles to polyhedron
              poly_n[c] = (int**)realloc((void*)poly_n[c],(poly_n[c][0][0]+lt+1)*sizeof(int*));
              for (i=0; i < lt; i++)
              {
                poly_n[c][next] = (int*)malloc(4*sizeof(int));
                poly_n[c][next][0] = 3;
                poly_n[c][next][1] = tnode.list[tri[i][2]];
                poly_n[c][next][2] = tnode.list[tri[i][1]];
                poly_n[c][next][3] = tnode.list[tri[i][0]];
                next++;
              }
              poly_n[c][0][0] += lt;
            }

            tnode.Redimension(0);
            delete[] tri;
            delete[] x;
            delete[] y;
            for (i=0; i < nbs[0]; i++)
              delete[] bs[0][i];
            delete[] nbs;
            delete[] bs[0];
            delete[] bs;
          }
        }

        // delete polygons from boundaries
        int *map = new int[ngon[b]];
        i=0;
        for (f=0; f < ngon[b]; f++)
        {
          map[f] = -1;
          if (ngon_n[b][f][0] > 0 && ngon_n[b][f][0] <= 4)
            map[f] = i++;
          if (ngon_n[b][f][0] > 4)
          {
            free(ngon_n[b][f]);
            ngon_n[b][f] = 0;
          }
        }
        // shift existing triangles and quads in ngon array down
        if (i > 0)
        {
          for (f=0; f < ngon[b]; f++)
          {
            if (map[f] >= 0 && map[f] != f)
            {
              ngon_n[b][map[f]] = (int*)realloc((void*)ngon_n[b][map[f]],(ngon_n[b][f][0]+1)*sizeof(int));
              ngon_n[b][map[f]][0] = ngon_n[b][f][0];
              for (j=1; j <= ngon_n[b][f][0]; j++)
                ngon_n[b][map[f]][j] = ngon_n[b][f][j];
              free(ngon_n[b][f]);
              ngon_n[b][f] = 0;
            }
          }
        } else
        {
          for (f=0; f < ngon[b]; f++)
          {
            if (ngon_n[b][f] != 0)
            {
              free(ngon_n[b][f]);
              ngon_n[b][f] = 0;
            }
          }
          free(ngon_n[b]);
          ngon_n[b] = 0;
          ngon[b] = 0;
        }
        delete[] map;
      }


      // convert internal face polygons
      for (c=0; c < nply; c++)
      {
        for (f=1; f <= poly_n[c][0][0]; f++)
        {
          if (poly_n[c][f][0] > 4)
          {
            // convert to triangles
            tnode.Redimension(0);
            Point cg = Point(0.0,0.0,0.0);
            for (i=1; i <= poly_n[c][f][0]; i++)
            {
              n0 = poly_n[c][f][i];
              tnode.Check_List(n0);
              cg += node[n0];
            }
            cg /= poly_n[c][f][0];
            norm = Vector(0.0,0.0,0.0);
            nbs = new int[1];
            bs = new int**[1];
            bs[0] = new int*[poly_n[c][f][0]];
            nbs[0] = 0;
            for (i=1; i <= poly_n[c][f][0]; i++)
            {
              if (i==1)
                n0 = poly_n[c][f][poly_n[c][f][0]];
              else
                n0 = poly_n[c][f][i-1];
              n1 = poly_n[c][f][i];
              bs[0][nbs[0]] = new int[2];
              bs[0][nbs[0]][0] = tnode.Index(n0);
              bs[0][nbs[0]][1] = tnode.Index(n1);
              nbs[0]++;
              v1 = Vector(cg,node[n0]);
              v2 = Vector(cg,node[n1]);
              norm += v1 % v2;
            }
            norm.normalize();
            x = new double[tnode.max];
            y = new double[tnode.max];
            n0 = poly_n[c][f][1];
            v1 = Vector(cg,node[n0]);
            dot = v1 * norm;
            v1 -= norm*dot;
            v1.normalize();
            v2 = norm % v1;
            v2.normalize();
            dsmin = 1.0e20;
            for (i=1; i <= poly_n[c][f][0]; i++)
            {
              if (i==1)
                n0 = poly_n[c][f][poly_n[c][f][0]];
              else
                n0 = poly_n[c][f][i-1];
              n1 = poly_n[c][f][i];
              dr = Vector(cg,node[n0]);
              int i0 = tnode.Index(n0);
              x[i0] = dr*v1;
              y[i0] = dr*v2;
              dr = Vector(cg,node[n1]);
              int i1 = tnode.Index(n1);
              x[i1] = dr*v1;
              y[i1] = dr*v2;
              double dx = x[i1]-x[i0];
              double dy = y[i1]-y[i0];
              double mag = sqrt(dx*dx+dy*dy);
              dsmin = MIN(dsmin,mag);
            }
            if (dsmin < 1.0e-15)
            {
              printf("\nCONVERT: coincident points detected in triangulation process!");
              fflush(stdout);
              exit(0);
            }
            int tdim = tnode.max*5;
            int (*tri)[3] = new int[tdim][3];
            int lt, t;

            if ((lt = trimesh(tnode.max,tdim,1,0,nbs,bs,x,y,tri)) > 0)
            {
              // find neighboring polyhedron
              int nc = -1;
              n0 = poly_n[c][f][1];
              for (i=0; i < p_hash[n0]->max && nc < 0; i++)
              {
                n = p_hash[n0]->list[i];
                if (n == c)
                  continue;
                int match = 1;
                for (j=2; j <= poly_n[c][f][0] && match; j++)
                {
                  n1 = poly_n[c][f][j];
                  match = p_hash[n1]->Is_In_List(n);
                }
                if (match)
                  nc = n;
              }
              if (nc < 0)
              {
                printf("\nCONVERT: Neighboring polyhedron for polygon face not found!");
                fflush(stdout);
                exit(0);
              }

              // find polygon face in neighbor
              int nf = -1;
              for (i=1; i <= poly_n[nc][0][0] && nf < 0; i++)
              {
                if (poly_n[c][f][0] != poly_n[nc][i][0])
                  continue;
                int *tag = new int[poly_n[nc][i][0]];
                for (j=0; j < poly_n[nc][i][0]; j++)
                  tag[j] = 0;
                for (j=1; j <= poly_n[c][f][0]; j++)
                  for (k=1; k <= poly_n[nc][i][0]; k++)
                    if (poly_n[c][f][j] == poly_n[nc][i][k])
                      tag[k-1] = 1;
                k=0;
                for (j=0; j < poly_n[nc][i][0]; j++)
                  k+=tag[j];
                if (k == poly_n[nc][i][0])
                  nf = i;
                delete[] tag;
              }
              if (nf < 0)
              {
                printf("\nCONVERT: Neighbor polyhedron face not found!");
                fflush(stdout);
                exit(0);
              }
              // replace existing face in cell c with last
              int last = poly_n[c][0][0];
              if (f < last)
              {
                poly_n[c][f] = (int*)realloc((void*)poly_n[c][f],(poly_n[c][last][0]+1)*sizeof(int));
                poly_n[c][f][0] = poly_n[c][last][0];
                for (j=1; j <= poly_n[c][last][0]; j++)
                  poly_n[c][f][j] = poly_n[c][last][j];
              }
              free(poly_n[c][last]);
              poly_n[c][0][0]--;
              poly_n[c] = (int**)realloc((void*)poly_n[c],(poly_n[c][0][0]+1)*sizeof(int*));

              // replace existing face in cell nc with last
              last = poly_n[nc][0][0];
              if (nf < last)
              {
                poly_n[nc][nf] = (int*)realloc((void*)poly_n[nc][nf],(poly_n[nc][last][0]+1)*sizeof(int));
                poly_n[nc][nf][0] = poly_n[nc][last][0];
                for (j=1; j <= poly_n[nc][last][0]; j++)
                  poly_n[nc][nf][j] = poly_n[nc][last][j];
              }
              free(poly_n[nc][last]);
              poly_n[nc][0][0]--;
              poly_n[nc] = (int**)realloc((void*)poly_n[nc],(poly_n[nc][0][0]+1)*sizeof(int*));

              int nextc = poly_n[c][0][0]+1;
              int nextn = poly_n[nc][0][0]+1;
              // now add triangles to both polyhedra
              poly_n[c] = (int**)realloc((void*)poly_n[c],(poly_n[c][0][0]+lt+1)*sizeof(int*));
              poly_n[nc] = (int**)realloc((void*)poly_n[nc],(poly_n[nc][0][0]+lt+1)*sizeof(int*));
              for (i=0; i < lt; i++)
              {
                poly_n[c][nextc] = (int*)malloc(4*sizeof(int));
                poly_n[c][nextc][0] = 3;
                poly_n[c][nextc][1] = tnode.list[tri[i][0]];
                poly_n[c][nextc][2] = tnode.list[tri[i][1]];
                poly_n[c][nextc][3] = tnode.list[tri[i][2]];
                nextc++;
                poly_n[nc][nextn] = (int*)malloc(4*sizeof(int));
                poly_n[nc][nextn][0] = 3;
                poly_n[nc][nextn][1] = tnode.list[tri[i][2]];
                poly_n[nc][nextn][2] = tnode.list[tri[i][1]];
                poly_n[nc][nextn][3] = tnode.list[tri[i][0]];
                nextn++;
              }
              poly_n[c][0][0] += lt;
              poly_n[nc][0][0] += lt;
            }
            tnode.Redimension(0);
            delete[] tri;
            delete[] x;
            delete[] y;
            for (i=0; i < nbs[0]; i++)
              delete[] bs[0][i];
            delete[] nbs;
            delete[] bs[0];
            delete[] bs;
          }
        }
      }

      for (n=0; n < nn; n++)
        delete p_hash[n];
      free(p_hash);
    }

    // count number of new nodes (centroids)
    node = (Point*)realloc((void*)node,(nn+nply)*sizeof(Point));
    int newpyr, newtet;
    newpyr = newtet = 0;
    for (c=0; c < nply; c++)
    {
      node[nn+c] = Point(0.0,0.0,0.0);
      k=0;
      for (i=1; i <= poly_n[c][0][0]; i++)
      {
        for (j=1; j <= poly_n[c][i][0]; j++)
        {
          n = poly_n[c][i][j];
          node[nn+c] += node[n];
          k++;
        }
        if (poly_n[c][i][0] == 3) newtet++;
        if (poly_n[c][i][0] == 4) newpyr++;
      }
      node[nn+c] /= MAX(1,k);
    }
    if (newtet > 0)
    {
      tet_n = (int**)realloc((void*)tet_n,(ntet+newtet)*sizeof(int*));
      for (c=ntet; c < ntet+newtet; c++)
        tet_n[c] = (int*)malloc(4*sizeof(int));
      if (tet_vc > 0)
        tet_vc = (int*)realloc((void*)tet_vc,(ntet+newtet)*sizeof(int));
    }
    if (newpyr > 0)
    {
      pyr_n = (int**)realloc((void*)pyr_n,(npyr+newpyr)*sizeof(int*));
      for (c=npyr; c < npyr+newpyr; c++)
        pyr_n[c] = (int*)malloc(5*sizeof(int));
      if (pyr_vc > 0)
        pyr_vc = (int*)realloc((void*)pyr_vc,(npyr+newpyr)*sizeof(int));
    }
    newtet=ntet;
    newpyr=npyr;
    for (c=0; c < nply; c++)
    {
      for (i=1; i <= poly_n[c][0][0]; i++)
      {
        if (poly_n[c][i][0] == 3)
        {
          tet_n[newtet][0] = poly_n[c][i][3];
          tet_n[newtet][1] = poly_n[c][i][2];
          tet_n[newtet][2] = poly_n[c][i][1];
          tet_n[newtet][3] = nn+c;
          if (tet_vc > 0 && poly_vc)
            tet_vc[newtet] = poly_vc[c];
          newtet++;
        }
        if (poly_n[c][i][0] == 4)
        {
          pyr_n[newpyr][0] = poly_n[c][i][4];
          pyr_n[newpyr][1] = poly_n[c][i][3];
          pyr_n[newpyr][2] = poly_n[c][i][2];
          pyr_n[newpyr][3] = poly_n[c][i][1];
          pyr_n[newpyr][4] = nn+c;
          if (pyr_vc > 0 && poly_vc)
            pyr_vc[newpyr] = poly_vc[c];
          newpyr++;
        }
      }
    }
    nn += nply;
    ntet = newtet;
    npyr = newpyr;
    for (i=0; i < nply; i++)
    {
      for (j=poly_n[i][0][0]; j >= 0; j--)
        if (poly_n[i][j] != 0) free(poly_n[i][j]);
      if (poly_n[i] != 0) free(poly_n[i]);
    }
    if (poly_n != 0) free(poly_n);
    if (poly_vc != 0) free(poly_vc);
    poly_n = 0;
    nply = 0;
  }
  
  if (strlen(vname) > 0)
  {
    if (nvol > 0)
    {
      // reset all vc numbers to 1 and store one vc name
      for (i=1; i < nvol; i++)
        free(vol_name[i]);
    } else
    {
      vol_name = (char**)malloc(sizeof(char*));
      vol_name[0] = (char*)malloc((strlen(vname)+1)*sizeof(char));
    }
    strcpy(vol_name[0],vname);
    nvol=1;
    for (i=0; i < ntet; i++)
      tet_vc[i] = 0;
    for (i=0; i < npyr; i++)
      pyr_vc[i] = 0;
    for (i=0; i < npri; i++)
      pri_vc[i] = 0;
    for (i=0; i < nhex; i++)
      hex_vc[i] = 0;
    for (i=0; i < nply; i++)
      poly_vc[i] = 0;
  } else
  {
    // determine unique volume id's
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

    if (nvol == 0)
    {
      nvol=vlist->max;
      vol_name = (char**)malloc(nvol*sizeof(char*));
      for (i=0; i < nvol; i++)
      {
        vol_name[i] = (char*)malloc(33*sizeof(char));
        sprintf(vol_name[i],"Volume_%i",i+1);
      }
    }

    // check for consistency between volume names and tags
    if (vlist->max != nvol)
    {
      printf("\nCONVERT: Number of volume tags, %i, does not equal number of volume names, %i.",vlist->max,nvol);
      fflush(stdout);
      exit(0);
    }
    delete vlist;

    // check for duplicate volume names
    // if found transfer elements to lower name
    for (b=0; b < nvol-1; b++)
    {
      c=b+1;
      while (c < nvol)
      {
        if (strcmp(vol_name[b],vol_name[c]) == 0)
        {
          printf("\nDuplicate volume names found for volumes %d & %d",b+1,c+1);
          printf("\nVolume %d is %s",b+1,vol_name[b]);
          printf("\nVolume %d is %s",c+1,vol_name[c]);
          printf("\nRetagging elements from volume %d to volume %d",c+1,b+1);
          fflush(stdout);

          for (i=0; i < ntet; i++)
            if (tet_vc[i] == c)
              tet_vc[i] = b;
          for (i=0; i < npyr; i++)
            if (pyr_vc[i] == c)
              pyr_vc[i] = b;
          for (i=0; i < npri; i++)
            if (pri_vc[i] == c)
              pri_vc[i] = b;
          for (i=0; i < nhex; i++)
            if (hex_vc[i] == c)
              hex_vc[i] = b;
          for (i=0; i < nply; i++)
            if (poly_vc[i] == c)
              poly_vc[i] = b;

          // shift remaining volumes up
          for (m=c; m < nvol-1; m++)
          {
            fprintf(stdout,"\nShifting volume %d to volume %d",m+2,m+1);
            fflush(stdout);
            vol_name[m][0] = '\0';
            strcpy(vol_name[m],vol_name[m+1]);
            for (i=0; i < ntet; i++)
              if (tet_vc[i] == m+1)
                tet_vc[i] = m;
            for (i=0; i < npyr; i++)
              if (pyr_vc[i] == m+1)
                pyr_vc[i] = m;
            for (i=0; i < npri; i++)
              if (pri_vc[i] == m+1)
                pri_vc[i] = m;
            for (i=0; i < nhex; i++)
              if (hex_vc[i] == m+1)
                hex_vc[i] = m;
            for (i=0; i < nply; i++)
              if (poly_vc[i] == m+1)
                poly_vc[i] = m;
          }
          // delete last boundary
          m = nvol-1;
          free(vol_name[m]);
          nvol--;
          vol_name = (char**)realloc((void*)vol_name,nvol*sizeof(char*));
        } else
          c++;
      }
    }
  }
    
  // check for duplicate boundary names
  // if found transfer elements to lower boundary and delete higher boundary
  for (b=0; b < nb-1; b++)
  {
    c=b+1;
    while (c < nb)
    {
      if (strcmp(b_name[b],b_name[c]) == 0)
      {
        printf("\nDuplicate boundary names found for boundaries %d & %d",b+1,c+1);
        printf("\nBoundary %d is %s",b+1,b_name[b]);
        printf("\nBoundary %d is %s",c+1,b_name[c]);
        printf("\nTransferring elements from boundary %d to boundary %d",c+1,b+1);
        fflush(stdout);

        // transfer boundary c to b
        if (nt[c] > 0)
        {
          t_n[b] = (int**)realloc((void*)t_n[b],(nt[b]+nt[c])*sizeof(int*));
          for (j=0,i=nt[b]; i < nt[b]+nt[c]; i++, j++)
          {
            t_n[b][i] = (int*)malloc(3*sizeof(int));
            t_n[b][i][0] = t_n[c][j][0];
            t_n[b][i][1] = t_n[c][j][1];
            t_n[b][i][2] = t_n[c][j][2];
          }
          nt[b] += nt[c];
        }
        if (nq[c] > 0)
        {
          q_n[b] = (int**)realloc((void*)q_n[b],(nq[b]+nq[c])*sizeof(int*));
          for (j=0,i=nq[b]; i < nq[b]+nq[c]; i++, j++)
          {
            q_n[b][i] = (int*)malloc(4*sizeof(int));
            q_n[b][i][0] = q_n[c][j][0];
            q_n[b][i][1] = q_n[c][j][1];
            q_n[b][i][2] = q_n[c][j][2];
            q_n[b][i][3] = q_n[c][j][3];
          }
          nq[b] += nq[c];
        }
        if (ngon[c] > 0)
        {
          ngon_n[b] = (int**)realloc((void*)ngon_n[b],(ngon[b]+ngon[c])*sizeof(int*));
          for (j=0,i=ngon[b]; i < ngon[b]+ngon[c]; i++, j++)
          {
            ngon_n[b][i] = (int*)malloc((ngon_n[c][j][0]+1)*sizeof(int));
            ngon_n[b][i][0] = ngon_n[c][j][0];
            for (k=1; k <= ngon_n[c][j][0]; k++)
              ngon_n[b][i][k] = ngon_n[b][j][k];
          }
        }
        
        // shift remaining boundaries up
        for (m=c; m < nb-1; m++)
        {
          fprintf(stdout,"\nShifting boundary %d to boundary %d",m+2,m+1);
          fflush(stdout);
          b_name[m][0] = '\0';
          strcpy(b_name[m],b_name[m+1]);
          for (i=0; i < nt[m]; i++)
          {
            free(t_n[m][i]);
            t_n[m][i] = 0;
          }
          for (i=0; i < nq[m]; i++)
          {
            free(q_n[m][i]);
            q_n[m][i] = 0;
          }
          for (i=0; i < ngon[m]; i++)
          {
            free(ngon_n[m][i]);
            ngon_n[m][i] = 0;
          }
          if (nt[m+1] > 0)
          {
            t_n[m] = (int**)realloc((void*)t_n[m],(nt[m+1])*sizeof(int*));
            for (i=0; i < nt[m+1]; i++)
            {
              t_n[m][i] = (int*)malloc(3*sizeof(int));
              t_n[m][i][0] = t_n[m+1][i][0];
              t_n[m][i][1] = t_n[m+1][i][1];
              t_n[m][i][2] = t_n[m+1][i][2];
            }
          } else if (nt[m] > 0)
          {
            free(t_n[m]);
            t_n[m] = 0;
          }
          nt[m] = nt[m+1];
          if (nq[m+1] > 0)
          {
            q_n[m] = (int**)realloc((void*)q_n[m],(nq[m+1])*sizeof(int*));
            for (i=0; i < nq[m+1]; i++)
            {
              q_n[m][i] = (int*)malloc(4*sizeof(int));
              q_n[m][i][0] = q_n[m+1][i][0];
              q_n[m][i][1] = q_n[m+1][i][1];
              q_n[m][i][2] = q_n[m+1][i][2];
              q_n[m][i][3] = q_n[m+1][i][3];
            }
          } else if (nq[m] > 0)
          {
            free(q_n[m]);
            q_n[m] = 0;
          }
          nq[m] = nq[m+1];
          if (ngon[m+1] > 0)
          {
            ngon_n[m] = (int**)realloc((void*)ngon_n[m],(ngon[m+1])*sizeof(int*));
            for (i=0; i < ngon[m+1]; i++)
            {
              ngon_n[m][i] = (int*)malloc((ngon_n[m+1][i][0]+1)*sizeof(int));
              ngon_n[m][i][0] = ngon_n[m+1][i][0];
              for (k=1; k <= ngon_n[m+1][i][0]; k++)
                ngon_n[m][i][k] = ngon_n[m+1][i][k];
            }
          } else if (ngon[m] > 0)
          {
            free(ngon_n[m]);
            ngon_n[m] = 0;
          }
          ngon[m] = ngon[m+1];
        }
        // delete last boundary
        m = nb-1;
        for (i=0; i < nt[m]; i++)
          free(t_n[m][i]);
        for (i=0; i < nq[m]; i++)
          free(q_n[m][i]);
        for (i=0; i < ngon[m]; i++)
          free(ngon_n[m][i]);
        if (nt[m] > 0)
          free(t_n[m]);
        if (nq[m] > 0)
          free(q_n[m]);
        if (ngon[m] > 0)
          free(ngon_n[m]);
        t_n[m] = 0;
        q_n[m] = 0;
        ngon_n[m] = 0;
        free(b_name[m]);
        nb--;
        b_name = (char**)realloc((void*)b_name,nb*sizeof(char*));
        nt = (int*)realloc((void*)nt,nb*sizeof(int));
        t_n = (int***)realloc((void*)t_n,nb*sizeof(int**));
        nq = (int*)realloc((void*)nq,nb*sizeof(int));
        q_n = (int***)realloc((void*)q_n,nb*sizeof(int**));
        ngon = (int*)realloc((void*)ngon,nb*sizeof(int));
        ngon_n = (int***)realloc((void*)ngon_n,nb*sizeof(int**));
      } else
        c++;
    }
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

  Write_Mesh(sname,b_name,nn,x,y,z,nb,nt,nq,ngon,t_n,q_n,ngon_n,ntet,tet_n,npyr,pyr_n,npri,pri_n,nhex,hex_n,nply,poly_n,
            nvol,vol_name,tet_vc,pyr_vc,pri_vc,hex_vc,poly_vc);

  // free up memory
  if (nvol > 0)
  {
    for (i=0; i < nvol; i++)
      free(vol_name[i]);
    free(vol_name);
  }
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
  for (i=0; i < npyr; i++)
    if (pyr_n[i] != 0) free(pyr_n[i]);
  if (pyr_n != 0) free(pyr_n);
  for (i=0; i < npri; i++)
    if (pri_n[i] != 0) free(pri_n[i]);
  if (pri_n != 0) free(pri_n);
  for (i=0; i < nhex; i++)
    if (hex_n[i] != 0) free(hex_n[i]);
  if (hex_n != 0) free(hex_n);
  for (i=0; i < nply; i++)
  {
    for (j=poly_n[i][0][0]; j >= 0; j--)
      if (poly_n[i][j] != 0) free(poly_n[i][j]);
    if (poly_n[i] != 0) free(poly_n[i]);
  }
  if (poly_n != 0) free(poly_n);
  if (node != 0) free(node);
  if (x != 0) free(x);
  if (y != 0) free(y);
  if (z != 0) free(z);
  if (poly_vc != 0) free(poly_vc);
  if (tet_vc != 0) free(tet_vc);
  if (pyr_vc != 0) free(pyr_vc);
  if (pri_vc != 0) free(pri_vc);
  if (hex_vc != 0) free(hex_vc);

  return(error);
}
