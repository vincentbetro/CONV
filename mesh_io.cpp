#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "Point.h"
#include "Vector.h"
#include "Util.h"
#include "Linked_List.h"
#include "List.h"
#ifdef HAVE_CGNS
#include "CGNS.h"
#endif
#include "SGIO.h"
#include "mesh_io.h"
#include "mesh_readers.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))

// Fieldview routines
#include "FV_routines.h"

#define CHUNK 100000

int Nastran(char filename[], int mode, int &nn, Point **node, int &nb, char ***b_name, int **nt, int ****t_n, 
             int **nq, int ****q_n, int **ngon, int ****ngon_n, int &ntet, int ***tet_n, int &npyr, int ***pyr_n,
             int &npri, int ***pri_n, int &nhex, int ***hex_n, int &nply, int ****poly_n)
{
  FILE* fp;
  const int bdim = 132;
  char buff[bdim], s[bdim], sx[bdim], sy[bdim], sz[bdim];
  int b, c, i, j, k, m, n, n0, n1, n2, n3, n4, n5, q, t;
  int counter;
  double x, y, z;

  fprintf(stdout,"\nFilename = <%s>",filename);

  switch (mode)
  {
    case -1:
      // Open file for read
      if ((fp = fopen(filename,"r")) == 0)
      {
        fprintf(stdout,"\nError opening file <%s>.",filename);
        exit(0);
      }

      nn = ntet = npyr = npri = nhex = nb = nply = 0;
      int shellmn, shellmx, shellnum;
      shellmn = 999999999;
      shellmx = shellnum = 0;
      // count number of nodes, triangles & quads
      while (fgets(buff,bdim,fp) != NULL)
      {
        sprintf(s,"%8s",buff);
        if (strstr(s,"GRID") != NULL)
          nn++;
        if (strstr(s,"CTETRA") != NULL)
          ntet++;
        if (strstr(s,"CHEXA") != NULL)
          nhex++;
        if (strstr(s,"CPENTA") != NULL)
        {
          sscanf(buff,"%s %d %d %d %d %d %d %d %d",sx,&i,&j,
                 &n0,&n1,&n2,&n3,&n4,&n5);
          if (n4 == n5)
            npyr++;
          else
            npri++;
        }
        if (strstr(s,"PSHELL") != NULL)
        {
          sscanf(buff,"%s %d",s,&i);
          shellmn = MIN(shellmn,i);
          shellmx = MAX(shellmx,i);
          shellnum++;
        }
      }

      if (shellnum > 0)
      {
        if (shellmx-shellmn+1 != shellnum)
        {
          fprintf(stdout,"\nNastran: shell count invalid!");
          exit(0);
        }
        nb = shellnum;
        (*nt) = (int*)malloc(nb*sizeof(int));
        (*nq) = (int*)malloc(nb*sizeof(int));
        (*ngon) = (int*)malloc(nb*sizeof(int));
        (*t_n) = (int***)malloc(nb*sizeof(int**));
        (*q_n) = (int***)malloc(nb*sizeof(int**));
        (*ngon_n) = (int***)malloc(nb*sizeof(int**));
        (*b_name) = (char**)malloc(nb*sizeof(char*));
        for (i=0; i < nb; i++)
        {
          (*nt)[i] = (*nq)[i] = (*ngon)[i] = 0;
          (*t_n)[i] = (*q_n)[i] = (*ngon_n)[i] = 0;
          (*b_name)[i] = (char*)malloc(33*sizeof(char));
          sprintf((*b_name)[i],"Boundary_%i",i+1);
        }

        rewind(fp);
      
        while (fgets(buff,bdim,fp) != NULL)
        {
          if (strstr(buff,"CTRIA3") != NULL)
          {
            sscanf(buff,"%s %d %d",s,&i,&j);
            j -= shellmn;
            if (j < 0 || j >= nb)
            {
              fprintf(stdout,"\nNastran: PSHELL out of bounds!");
              exit(0);
            }
            (*nt)[j]++;
          }
          if (strstr(buff,"CQUAD4") != NULL)
          {
            sscanf(buff,"%s %d %d",s,&i,&j);
            j -= shellmn;
            if (j < 0 || j >= nb)
            {
              fprintf(stdout,"\nNastran: PSHELL out of bounds!");
              exit(0);
            }
            (*nq)[j]++;
          }
        }
        for (b=0; b < nb; b++)
        {
          if ((*nt)[b] > 0)
          {
            fprintf(stdout,"\nBoundary %d, number of triangles = %d",b+1,(*nt)[b]);
            (*t_n)[b] = (int**)malloc((*nt)[b]*sizeof(int*));
            for (i=0; i < (*nt)[b]; i++)
              (*t_n)[b][i] = (int*)malloc(3*sizeof(int));
          }
          if ((*nq)[b] > 0)
          {
            fprintf(stdout,"\nBoundary %d, number of quadrilaterals = %d",b+1,(*nq)[b]);
            (*q_n)[b] = (int**)malloc((*nq)[b]*sizeof(int*));
            for (i=0; i < (*nq)[b]; i++)
              (*q_n)[b][i] = (int*)malloc(4*sizeof(int));
          }
        }
      }

      (*node) = (Point*)malloc(nn*sizeof(Point));
      fprintf(stdout,"\nNumber of nodes = %d",nn);
      
      if (ntet > 0)
      {
        fprintf(stdout,"\nNumber of tetrahedra = %d",ntet);
        (*tet_n) = (int**)malloc(ntet*sizeof(int*));
        for (n=0; n < ntet; n++)
          (*tet_n)[n] = (int*)malloc(4*sizeof(int));
      }
      if (npyr > 0)
      {
        fprintf(stdout,"\nNumber of pyramid = %d",npyr);
        (*pyr_n) = (int**)malloc(npyr*sizeof(int*));
        for (n=0; n < npyr; n++)
          (*pyr_n)[n] = (int*)malloc(5*sizeof(int));
      }
      if (npri > 0)
      {
        fprintf(stdout,"\nNumber of prism = %d",npri);
        (*pri_n) = (int**)malloc(npri*sizeof(int*));
        for (n=0; n < npri; n++)
          (*pri_n)[n] = (int*)malloc(6*sizeof(int));
      }
      if (nhex > 0)
      {
        fprintf(stdout,"\nNumber of hexahedra = %d",nhex);
        (*hex_n) = (int**)malloc(nhex*sizeof(int*));
        for (n=0; n < nhex; n++)
          (*hex_n)[n] = (int*)malloc(8*sizeof(int));
      }

      rewind(fp);

      nn = ntet = npyr = npri = nhex = nply = 0;
      for (i=0; i < nb; i++)
        (*nt)[i]=(*nq)[i]=(*ngon)[i]=0;

      // read nodes, volume elements and boundary elements
      while (fgets(buff,bdim,fp) != NULL)
      {
        if (strstr(buff,"GRID") != NULL)
        {
          //if ((i=sscanf(buff,"%4s%20i%8s%8s%8s",s,&n,sx,sy,sz)) != 5)
          sscanf(buff,"%s%i",s,&n);
          for (i=0; i < 8; i++)
          {
            sx[i] = buff[24+i];
            sy[i] = buff[32+i];
            sz[i] = buff[40+i];
          }
          sx[9] = '\0';
          sy[9] = '\0';
          sz[9] = '\0';
          x = atof(sx);
          y = atof(sy);
          z = atof(sz);
          nn = MAX(nn,n);
          n--;
          (*node)[n][0] = x;
          (*node)[n][1] = y;
          (*node)[n][2] = z;
        }
        if (strstr(buff,"CTETRA") != NULL)
        {
          sscanf(buff,"CTETRA %d %d %d %d %d",&j,&((*tet_n)[ntet][0]),
                                                    &((*tet_n)[ntet][1]),
                                                    &((*tet_n)[ntet][2]),
                                                    &((*tet_n)[ntet][3]));
          (*tet_n)[ntet][0]--;
          (*tet_n)[ntet][1]--;
          (*tet_n)[ntet][2]--;
          (*tet_n)[ntet][3]--;
          ntet++;
        }
        if (strstr(buff,"CPENTA") != NULL)
        {
          sscanf(buff,"%s %d %d %d %d %d %d %d %d",s,&j,&k,
                 &n0,&n1,&n2,&n3,&n4,&n5);
          if (n4 == n5)
          {
            (*pyr_n)[npyr][0] = --n0;
            (*pyr_n)[npyr][1] = --n3;
            (*pyr_n)[npyr][2] = --n4;
            (*pyr_n)[npyr][3] = --n1;
            (*pyr_n)[npyr][4] = --n2;
            npyr++;
          } else {
            (*pri_n)[npri][0] = --n0;
            (*pri_n)[npri][1] = --n1;
            (*pri_n)[npri][2] = --n2;
            (*pri_n)[npri][3] = --n3;
            (*pri_n)[npri][4] = --n4;
            (*pri_n)[npri][5] = --n5;
            npri++;
          }
        }
        if (strstr(buff,"CHEXA") != NULL)
        {
          if ((i=sscanf(buff,"CHEXA %d %d %d %d %d %d %d %d %d %d",&j,&k,
                        &((*hex_n)[nhex][0]), &((*hex_n)[nhex][1]),
                        &((*hex_n)[nhex][2]), &((*hex_n)[nhex][3]),
                        &((*hex_n)[nhex][4]), &((*hex_n)[nhex][5]),
                        &((*hex_n)[nhex][6]), &((*hex_n)[nhex][7]))) != 10)
          {
            if (i == 8)
            {
              fgets(buff,bdim,fp);
              sscanf(buff,"%d %d", &((*hex_n)[nhex][6]), &((*hex_n)[nhex][7]));
            } else
            {
              fprintf(stdout,"\nNastran: ERROR READING HEX ELEMENTS!");
              exit(0);
            }
          }
          (*hex_n)[nhex][0]--;
          (*hex_n)[nhex][1]--;
          (*hex_n)[nhex][2]--;
          (*hex_n)[nhex][3]--;
          (*hex_n)[nhex][4]--;
          (*hex_n)[nhex][5]--;
          (*hex_n)[nhex][6]--;
          (*hex_n)[nhex][7]--;
          nhex++;
        }
        if (strstr(buff,"CTRIA3") != NULL)
        {
          sscanf(buff,"%s %d %d %d %d %d",s,&j,&i,&n0,&n1,&n2);
          j--;
          i -= shellmn;
          (*t_n)[i][(*nt)[i]][0]= --n0;
          (*t_n)[i][(*nt)[i]][1]= --n1;
          (*t_n)[i][(*nt)[i]][2]= --n2;
          (*nt)[i]++;
        }
        if (strstr(buff,"CQUAD4") != NULL)
        {
          sscanf(buff,"%s %d %d %d %d %d %d",s,&j,&i,&n0,&n1,&n2,&n3);
          j--;
          i -= shellmn;
          (*q_n)[i][(*nq)[i]][0]= --n0;
          (*q_n)[i][(*nq)[i]][1]= --n1;
          (*q_n)[i][(*nq)[i]][2]= --n2;
          (*q_n)[i][(*nq)[i]][3]= --n3;
          (*nq)[i]++;
        }
      }

      fclose(fp);

      // determine boundary faces by storing un-duplicated tris and quads
      if (nb == 0)
      {
        Linked_List **tri_hash, **quad_hash;
        Linked_Node *hd0, *hd1, *hd2, *hd3;

        tri_hash = (Linked_List**)malloc(nn*sizeof(Linked_List*));
        quad_hash = (Linked_List**)malloc(nn*sizeof(Linked_List*));

        for (n=0; n < nn; n++)
        {
          tri_hash[n] = new Linked_List();
          quad_hash[n] = new Linked_List();
        }

        nb = 1;
        (*nt) = (int*)malloc(nb*sizeof(int));
        (*nq) = (int*)malloc(nb*sizeof(int));
        (*ngon) = (int*)malloc(nb*sizeof(int));
        (*t_n) = (int***)malloc(nb*sizeof(int**));
        (*q_n) = (int***)malloc(nb*sizeof(int**));
        (*ngon_n) = (int***)malloc(nb*sizeof(int**));
        (*b_name) = (char**)malloc(nb*sizeof(char*));
        for (i=0; i < nb; i++)
        {
          (*nt)[i] = (*nq)[i] = (*ngon)[i] = 0;
          (*t_n)[i] = (*q_n)[i] = (*ngon_n)[i] = 0;
          (*b_name)[i] = (char*)malloc(33*sizeof(char));
          sprintf((*b_name)[i],"Boundary_%i",i+1);
        }

        // dimension for maximum number of faces
        (*nt)[0] = ntet*4+npyr*4+npri*2;
        (*nq)[0] = npyr+npri*3+nhex*6;
        if ((*nt)[0] > 0)
          (*t_n)[0] = (int**)malloc((*nt)[0]*sizeof(int*));
        if ((*nq)[0] > 0)
          (*q_n)[0] = (int**)malloc((*nq)[0]*sizeof(int*));
        for (i=0; i < (*nt)[0]; i++)
          (*t_n)[0][i] = (int*)malloc(3*sizeof(int));
        for (i=0; i < (*nq)[0]; i++)
          (*q_n)[0][i] = (int*)malloc(4*sizeof(int));

        int tdim = (*nt)[0];
        int qdim = (*nq)[0];

        // add/delete faces
        (*nt)[0]=(*nq)[0]=0;
        for (c=0; c < ntet; c++)
        {
          for (i=0; i < 4; i++)
          {
            switch(i)
            {
              case 0:
                n0 = (*tet_n)[c][0];
                n1 = (*tet_n)[c][1];
                n2 = (*tet_n)[c][2];
                break;
              case 1:
                n0 = (*tet_n)[c][0];
                n1 = (*tet_n)[c][3];
                n2 = (*tet_n)[c][1];
                break;
              case 2:
                n0 = (*tet_n)[c][1];
                n1 = (*tet_n)[c][3];
                n2 = (*tet_n)[c][2];
                break;
              case 3:
                n0 = (*tet_n)[c][2];
                n1 = (*tet_n)[c][3];
                n2 = (*tet_n)[c][0];
                break;
            }
            t=-1;
            hd0 = tri_hash[n0]->head;
            while (hd0 && t < 0)
            {
              n = hd0->data;
              hd1 = tri_hash[n1]->head;
              while (hd1 && t < 0)
              {
                if (hd1->data == n)
                {
                  hd2 = tri_hash[n2]->head;
                  while (hd2 && t < 0)
                  {
                    if (hd2->data == n)
                      t = n;
                    hd2 = hd2->next;
                  }
                }
                hd1 = hd1->next;
              }
              hd0 = hd0->next;
            }
            if (t >= 0)
            {
              // delete face from list
              tri_hash[n0]->Remove(t);
              tri_hash[n1]->Remove(t);
              tri_hash[n2]->Remove(t);
              // move last triangle to current location
              m = (*nt)[0]-1;
              n0=(*t_n)[0][m][0];
              n1=(*t_n)[0][m][1];
              n2=(*t_n)[0][m][2];
              tri_hash[n0]->Replace(m,t);
              tri_hash[n1]->Replace(m,t);
              tri_hash[n2]->Replace(m,t);
              (*t_n)[0][t][0]=(*t_n)[0][m][0];
              (*t_n)[0][t][1]=(*t_n)[0][m][1];
              (*t_n)[0][t][2]=(*t_n)[0][m][2];
              (*nt)[0]--;
            } else
            {
              // add face to list
              (*t_n)[0][(*nt)[0]][0] = n0;
              (*t_n)[0][(*nt)[0]][1] = n1;
              (*t_n)[0][(*nt)[0]][2] = n2;
              tri_hash[n0]->Insert((*nt)[0]);
              tri_hash[n1]->Insert((*nt)[0]);
              tri_hash[n2]->Insert((*nt)[0]);
              (*nt)[0]++;
            }
          }
        }

        for (c=0; c < npyr; c++)
        {
          for (i=1; i < 5; i++)
          {
            switch(i)
            {
              case 1:
                n0 = (*pyr_n)[c][0];
                n1 = (*pyr_n)[c][4];
                n2 = (*pyr_n)[c][1];
                break;
              case 2:
                n0 = (*pyr_n)[c][1];
                n1 = (*pyr_n)[c][4];
                n2 = (*pyr_n)[c][2];
                break;
              case 3:
                n0 = (*pyr_n)[c][2];
                n1 = (*pyr_n)[c][4];
                n2 = (*pyr_n)[c][3];
                break;
              case 4:
                n0 = (*pyr_n)[c][3];
                n1 = (*pyr_n)[c][4];
                n2 = (*pyr_n)[c][0];
                break;
            }
            t=-1;
            hd0 = tri_hash[n0]->head;
            while (hd0 && t < 0)
            {
              n = hd0->data;
              hd1 = tri_hash[n1]->head;
              while (hd1 && t < 0)
              {
                if (hd1->data == n)
                {
                  hd2 = tri_hash[n2]->head;
                  while (hd2 && t < 0)
                  {
                    if (hd2->data == n)
                      t = n;
                    hd2 = hd2->next;
                  }
                }
                hd1 = hd1->next;
              }
              hd0 = hd0->next;
            }
            if (t >= 0)
            {
              // delete face from list
              tri_hash[n0]->Remove(t);
              tri_hash[n1]->Remove(t);
              tri_hash[n2]->Remove(t);
              // move last triangle to current location
              m = (*nt)[0]-1;
              n0=(*t_n)[0][m][0];
              n1=(*t_n)[0][m][1];
              n2=(*t_n)[0][m][2];
              tri_hash[n0]->Replace(m,t);
              tri_hash[n1]->Replace(m,t);
              tri_hash[n2]->Replace(m,t);
              (*t_n)[0][t][0]=(*t_n)[0][m][0];
              (*t_n)[0][t][1]=(*t_n)[0][m][1];
              (*t_n)[0][t][2]=(*t_n)[0][m][2];
              (*nt)[0]--;
            } else
            {
              // add face to list
              (*t_n)[0][(*nt)[0]][0] = n0;
              (*t_n)[0][(*nt)[0]][1] = n1;
              (*t_n)[0][(*nt)[0]][2] = n2;
              tri_hash[n0]->Insert((*nt)[0]);
              tri_hash[n1]->Insert((*nt)[0]);
              tri_hash[n2]->Insert((*nt)[0]);
              (*nt)[0]++;
            }
          }
        }

        for (c=0; c < npri; c++)
        {
          for (i=3; i < 5; i++)
          {
            switch(i)
            {
              case 3:
                n0 = (*pri_n)[c][0];
                n1 = (*pri_n)[c][1];
                n2 = (*pri_n)[c][2];
                break;
              case 4:
                n0 = (*pri_n)[c][3];
                n1 = (*pri_n)[c][5];
                n2 = (*pri_n)[c][4];
                break;
            }
            t=-1;
            hd0 = tri_hash[n0]->head;
            while (hd0 && t < 0)
            {
              n = hd0->data;
              hd1 = tri_hash[n1]->head;
              while (hd1 && t < 0)
              {
                if (hd1->data == n)
                {
                  hd2 = tri_hash[n2]->head;
                  while (hd2 && t < 0)
                  {
                    if (hd2->data == n)
                      t = n;
                    hd2 = hd2->next;
                  }
                }
                hd1 = hd1->next;
              }
              hd0 = hd0->next;
            }
            if (t >= 0)
            {
              // delete face from list
              tri_hash[n0]->Remove(t);
              tri_hash[n1]->Remove(t);
              tri_hash[n2]->Remove(t);
              // move last triangle to current location
              m = (*nt)[0]-1;
              n0=(*t_n)[0][m][0];
              n1=(*t_n)[0][m][1];
              n2=(*t_n)[0][m][2];
              tri_hash[n0]->Replace(m,t);
              tri_hash[n1]->Replace(m,t);
              tri_hash[n2]->Replace(m,t);
              (*t_n)[0][t][0]=(*t_n)[0][m][0];
              (*t_n)[0][t][1]=(*t_n)[0][m][1];
              (*t_n)[0][t][2]=(*t_n)[0][m][2];
              (*nt)[0]--;
            } else
            {
              // add face to list
              (*t_n)[0][(*nt)[0]][0] = n0;
              (*t_n)[0][(*nt)[0]][1] = n1;
              (*t_n)[0][(*nt)[0]][2] = n2;
              tri_hash[n0]->Insert((*nt)[0]);
              tri_hash[n1]->Insert((*nt)[0]);
              tri_hash[n2]->Insert((*nt)[0]);
              (*nt)[0]++;
            }
          }
        }

        for (c=0; c < npyr; c++)
        {
          n0 = (*pyr_n)[c][0];
          n1 = (*pyr_n)[c][1];
          n2 = (*pyr_n)[c][2];
          n3 = (*pyr_n)[c][3];
          q=-1;
          hd0 = quad_hash[n0]->head;
          while (hd0 && q < 0)
          {
            n = hd0->data;
            hd1 = quad_hash[n1]->head;
            while (hd1 && q < 0)
            {
              if (hd1->data == n)
              {
                hd2 = quad_hash[n2]->head;
                while (hd2 && q < 0)
                {
                  if (hd2->data == n)
                  {
                    hd3 = quad_hash[n3]->head;
                    while (hd3 && q < 0)
                    {
                      if (hd3->data == n)
                        q = n;
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
          if (q >= 0)
          {
            // delete face from list
            quad_hash[n0]->Remove(q);
            quad_hash[n1]->Remove(q);
            quad_hash[n2]->Remove(q);
            quad_hash[n3]->Remove(q);
            // move last quadrilateral to current location
            m = (*nq)[0]-1;
            n0=(*q_n)[0][m][0];
            n1=(*q_n)[0][m][1];
            n2=(*q_n)[0][m][2];
            n3=(*q_n)[0][m][3];
            quad_hash[n0]->Replace(m,q);
            quad_hash[n1]->Replace(m,q);
            quad_hash[n2]->Replace(m,q);
            quad_hash[n3]->Replace(m,q);
            (*q_n)[0][q][0]=(*q_n)[0][m][0];
            (*q_n)[0][q][1]=(*q_n)[0][m][1];
            (*q_n)[0][q][2]=(*q_n)[0][m][2];
            (*q_n)[0][q][3]=(*q_n)[0][m][3];
            (*nq)[0]--;
          } else
          {
            // add face to list
            (*q_n)[0][(*nq)[0]][0] = n0;
            (*q_n)[0][(*nq)[0]][1] = n1;
            (*q_n)[0][(*nq)[0]][2] = n2;
            (*q_n)[0][(*nq)[0]][3] = n3;
            quad_hash[n0]->Insert((*nq)[0]);
            quad_hash[n1]->Insert((*nq)[0]);
            quad_hash[n2]->Insert((*nq)[0]);
            quad_hash[n3]->Insert((*nq)[0]);
            (*nq)[0]++;
          }
        }

        for (c=0; c < npri; c++)
        {
          for (i=0; i < 3; i++)
          {
            switch(i)
            {
              case 0:
                n0 = (*pri_n)[c][0];
                n1 = (*pri_n)[c][3];
                n2 = (*pri_n)[c][4];
                n3 = (*pri_n)[c][1];
                break;
              case 1:
                n0 = (*pri_n)[c][1];
                n1 = (*pri_n)[c][4];
                n2 = (*pri_n)[c][5];
                n3 = (*pri_n)[c][2];
                break;
              case 2:
                n0 = (*pri_n)[c][2];
                n1 = (*pri_n)[c][5];
                n2 = (*pri_n)[c][3];
                n3 = (*pri_n)[c][0];
                break;
            }
            q=-1;
            hd0 = quad_hash[n0]->head;
            while (hd0 && q < 0)
            {
              n = hd0->data;
              hd1 = quad_hash[n1]->head;
              while (hd1 && q < 0)
              {
                if (hd1->data == n)
                {
                  hd2 = quad_hash[n2]->head;
                  while (hd2 && q < 0)
                  {
                    if (hd2->data == n)
                    {
                      hd3 = quad_hash[n3]->head;
                      while (hd3 && q < 0)
                      {
                        if (hd3->data == n)
                          q = n;
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
            if (q >= 0)
            {
              // delete face from list
              quad_hash[n0]->Remove(q);
              quad_hash[n1]->Remove(q);
              quad_hash[n2]->Remove(q);
              quad_hash[n3]->Remove(q);
              // move last triangle to current location
              m = (*nq)[0]-1;
              n0=(*q_n)[0][m][0];
              n1=(*q_n)[0][m][1];
              n2=(*q_n)[0][m][2];
              n3=(*q_n)[0][m][3];
              quad_hash[n0]->Replace(m,q);
              quad_hash[n1]->Replace(m,q);
              quad_hash[n2]->Replace(m,q);
              quad_hash[n3]->Replace(m,q);
              (*q_n)[0][q][0]=(*q_n)[0][m][0];
              (*q_n)[0][q][1]=(*q_n)[0][m][1];
              (*q_n)[0][q][2]=(*q_n)[0][m][2];
              (*q_n)[0][q][3]=(*q_n)[0][m][3];
              (*nq)[0]--;
            } else
            {
              // add face to list
              (*q_n)[0][(*nq)[0]][0] = n0;
              (*q_n)[0][(*nq)[0]][1] = n1;
              (*q_n)[0][(*nq)[0]][2] = n2;
              (*q_n)[0][(*nq)[0]][3] = n3;
              quad_hash[n0]->Insert((*nq)[0]);
              quad_hash[n1]->Insert((*nq)[0]);
              quad_hash[n2]->Insert((*nq)[0]);
              quad_hash[n3]->Insert((*nq)[0]);
              (*nq)[0]++;
            }
          }
        }

        for (c=0; c < nhex; c++)
        {
          for (i=0; i < 6; i++)
          {
            switch(i)
            {
              case 0:
                n0 = (*hex_n)[c][0];
                n1 = (*hex_n)[c][1];
                n2 = (*hex_n)[c][2];
                n3 = (*hex_n)[c][3];
                break;
              case 1:
                n0 = (*hex_n)[c][0];
                n1 = (*hex_n)[c][4];
                n2 = (*hex_n)[c][5];
                n3 = (*hex_n)[c][1];
                break;
              case 2:
                n0 = (*hex_n)[c][1];
                n1 = (*hex_n)[c][5];
                n2 = (*hex_n)[c][6];
                n3 = (*hex_n)[c][2];
                break;
              case 3:
                n0 = (*hex_n)[c][2];
                n1 = (*hex_n)[c][6];
                n2 = (*hex_n)[c][7];
                n3 = (*hex_n)[c][3];
                break;
              case 4:
                n0 = (*hex_n)[c][3];
                n1 = (*hex_n)[c][7];
                n2 = (*hex_n)[c][4];
                n3 = (*hex_n)[c][0];
                break;
              case 5:
                n0 = (*hex_n)[c][4];
                n1 = (*hex_n)[c][7];
                n2 = (*hex_n)[c][6];
                n3 = (*hex_n)[c][5];
                break;
            }
            q=-1;
            hd0 = quad_hash[n0]->head;
            while (hd0 && q < 0)
            {
              n = hd0->data;
              hd1 = quad_hash[n1]->head;
              while (hd1 && q < 0)
              {
                if (hd1->data == n)
                {
                  hd2 = quad_hash[n2]->head;
                  while (hd2 && q < 0)
                  {
                    if (hd2->data == n)
                    {
                      hd3 = quad_hash[n3]->head;
                      while (hd3 && q < 0)
                      {
                        if (hd3->data == n)
                          q = n;
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
            if (q >= 0)
            {
              // delete face from list
              quad_hash[n0]->Remove(q);
              quad_hash[n1]->Remove(q);
              quad_hash[n2]->Remove(q);
              quad_hash[n3]->Remove(q);
              // move last triangle to current location
              m = (*nq)[0]-1;
              n0=(*q_n)[0][m][0];
              n1=(*q_n)[0][m][1];
              n2=(*q_n)[0][m][2];
              n3=(*q_n)[0][m][3];
              quad_hash[n0]->Replace(m,q);
              quad_hash[n1]->Replace(m,q);
              quad_hash[n2]->Replace(m,q);
              quad_hash[n3]->Replace(m,q);
              (*q_n)[0][q][0]=(*q_n)[0][m][0];
              (*q_n)[0][q][1]=(*q_n)[0][m][1];
              (*q_n)[0][q][2]=(*q_n)[0][m][2];
              (*q_n)[0][q][3]=(*q_n)[0][m][3];
              (*nq)[0]--;
            } else
            {
              // add face to list
              (*q_n)[0][(*nq)[0]][0] = n0;
              (*q_n)[0][(*nq)[0]][1] = n1;
              (*q_n)[0][(*nq)[0]][2] = n2;
              (*q_n)[0][(*nq)[0]][3] = n3;
              quad_hash[n0]->Insert((*nq)[0]);
              quad_hash[n1]->Insert((*nq)[0]);
              quad_hash[n2]->Insert((*nq)[0]);
              quad_hash[n3]->Insert((*nq)[0]);
              (*nq)[0]++;
            }
          }
        }

        fprintf(stdout,"\nNumber of boundary triangles      = %d",(*nt)[0]);
        fprintf(stdout,"\nNumber of boundary quadrilaterals = %d",(*nq)[0]);

        // free up extra face memory
        for (t=(*nt)[0]; t < tdim; t++)
          free((*t_n)[0][t]);
        for (q=(*nq)[0]; q < qdim; q++)
          free((*q_n)[0][q]);
        (*t_n)[0] = (int**)realloc((void*)(*t_n)[0],(*nt)[0]*sizeof(int*));
        (*q_n)[0] = (int**)realloc((void*)(*q_n)[0],(*nq)[0]*sizeof(int*));

        for (n=0; n < nn; n++)
        {
          delete tri_hash[n];
          delete quad_hash[n];
        }
        free(tri_hash);
        free(quad_hash);
      }

      break;
    case 1:
      // Open file for write
      if ((fp = fopen(filename,"w")) == 0)
      {
        fprintf(stdout,"\nError opening file <%s>.",filename);
        exit(0);
      }

      counter=1;
      fprintf(fp,"PSOLID,%7d,1\n",counter);
      for (b=0; b < nb; b++)
      {
        counter++;
        fprintf(fp,"PSHELL,%7d,1,1,,,,,0\n",counter);
      }
      for (n=0; n < nn; n++)
      {
        fprintf(fp,"GRID  ,%7d,,%16.9e,%16.9e,%16.9e\n",n+1,(*node)[n][0],
                                                            (*node)[n][1],
                                                            (*node)[n][2]);
      }

      // volume elements
      counter=0;
      for (c=0; c < nhex; c++)
      {
        counter++;
        fprintf(fp,"CHEXA ,%7d,1,%7d,%7d,%7d,%7d,%7d,%7d\n      ,          %7d,%7d\n",
                    counter,(*hex_n)[c][0]+1,
                            (*hex_n)[c][1]+1,
                            (*hex_n)[c][2]+1,
                            (*hex_n)[c][3]+1,
                            (*hex_n)[c][4]+1,
                            (*hex_n)[c][5]+1,
                            (*hex_n)[c][6]+1,
                            (*hex_n)[c][7]+1);
      }

      for (c=0; c < npri; c++)
      {
        counter++;
        fprintf(fp,"CPENTA,%7d,1,%7d,%7d,%7d,%7d,%7d,%7d\n",
                    counter,(*pri_n)[c][0]+1,
                            (*pri_n)[c][1]+1,
                            (*pri_n)[c][2]+1,
                            (*pri_n)[c][3]+1,
                            (*pri_n)[c][4]+1,
                            (*pri_n)[c][5]+1);
      }

      for (c=0; c < ntet; c++)
      {
        counter++;
        fprintf(fp,"CTETRA,%7d,1,%7d,%7d,%7d,%7d\n",
                    counter,(*tet_n)[c][0]+1,
                            (*tet_n)[c][1]+1,
                            (*tet_n)[c][2]+1,
                            (*tet_n)[c][3]+1);
      }

      for (c=0; c < npyr; c++)
      {
        counter++;
        fprintf(fp,"CPENTA,%7d,1,%7d,%7d,%7d,%7d,%7d,%7d\n",
                    counter,(*pyr_n)[c][0]+1,
                            (*pyr_n)[c][3]+1,
                            (*pyr_n)[c][4]+1,
                            (*pyr_n)[c][1]+1,
                            (*pyr_n)[c][2]+1,
                            (*pyr_n)[c][2]+1);
      }

      // boundary elements
      for (b=0; b < nb; b++)
      {
        for (c=0; c < (*nt)[b]; c++)
        {
          counter++;
          fprintf(fp,"CTRIA3,%7d,%3d,%7d,%7d,%7d,,0.000000\n",
                      counter,b+2,(*t_n)[b][c][0]+1,
                                  (*t_n)[b][c][1]+1,
                                  (*t_n)[b][c][2]+1);
        }
      }
      for (b=0; b < nb; b++)
      {
        for (c=0; c < (*nq)[b]; c++)
        {
          counter++;
          fprintf(fp,"CQUAD4,%7d,%3d,%7d,%7d,%7d,%7d,,0.000000\n",
                      counter,b+2,(*q_n)[b][c][0]+1,
                                  (*q_n)[b][c][1]+1,
                                  (*q_n)[b][c][2]+1,
                                  (*q_n)[b][c][3]+1);

        }
      }

      fclose(fp);

      break;
  }

  return(0);
}

int Ensight(char filename[], int mode, int &nn, Point **node, int &nb, char ***b_name, int **nt, int ****t_n, 
             int **nq, int ****q_n, int **ngon, int ****ngon_n, int &ntet, int ***tet_n, int &npyr, int ***pyr_n,
             int &npri, int ***pri_n, int &nhex, int ***hex_n, int &nply, int ****poly_n)
{
  FILE *fp, *casefp;
  int const bdim = 80;
  char buff[bdim], gname[bdim];
  int b, i, j, k, n, part;
  float fx, fy, fz;
  bool ascii = false;

  //
  // read/write Ensight files
  //
  switch(mode)
  {
    case -1: // Read

      break;

    case 1: // Write
      fprintf(stdout,"\nFilename = <%s>",filename);
      // Open file for write
      if ((casefp = fopen(filename,"w")) == 0)
      {
        fprintf(stdout,"\nEnsight: Error opening file <%s>.",filename);
        exit(0);
      }

      fprintf(casefp,"# Ensight Gold Case file: %s\n\n",filename);

      fprintf(casefp,"FORMAT\n\n");
      fprintf(casefp,"type:  ensight gold\n\n");

      sprintf(buff,"%s",filename);
      char *ptr = strstr(buff,".case");
      if (ptr == NULL)
      {
        fprintf(stdout,"\nEnsight suffix <.case> not found in file name!");
        fflush(stdout);
        exit(0);
      } else
        *ptr = '\0';
      sprintf(gname,"%s.geo",buff);

      fprintf(casefp,"GEOMETRY\n\n");
      fprintf(casefp,"model:     %s\n\n",gname);

      fclose(casefp);

      fprintf(stdout,"\nGeometry filename = <%s>",gname);
      // Open file for write
      if (ascii)
      {
        if ((fp = fopen(gname,"w")) == 0)
        {
          fprintf(stdout,"\nEnsight: Error opening geometry file <%s>.",gname);
          exit(0);
        }
      } else
        if ((fp = fopen(gname,"wb")) == 0)
        {
          fprintf(stdout,"\nEnsight: Error opening geometry file <%s>.",gname);
          exit(0);
        }

      buff[sizeof(buff)-1] = '\0';
      if (ascii)
      {
        fprintf(fp,"Ensight Geometry file written from Convert.\n");
        fprintf(fp,"Case file name = %s\n",filename);
      } else
      {
        strncpy(buff,"C Binary",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
        strncpy(buff,"Ensight Geometry file written from Convert.",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
        sprintf(buff,"Case file name = %s",filename);
        fwrite(buff,sizeof(char),80,fp);
      }

      if (ascii)
      {
        fprintf(fp,"node id given\n");
        fprintf(fp,"element id given\n");
      } else
      {
        strncpy(buff,"node id given",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
        strncpy(buff,"element id given",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
      }

      part = 0;
      if (ascii)
      {
        fprintf(fp,"part\n");
        fprintf(fp,"%10d\n",++part);
        fprintf(fp,"3D uns-elements\n");
      } else
      {
        strncpy(buff,"part",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
        ++part;
        fwrite(&part,sizeof(int),1,fp);
        strncpy(buff,"3D uns-elements",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
      }

      if (ascii)
      {
        fprintf(fp,"coordinates\n");
        fprintf(fp,"%10d\n",nn);
        for (n=0; n < nn; n++)
          fprintf(fp,"%10d\n",n+1);
        for (n=0; n < nn; n++)
          fprintf(fp,"%12.5e\n",(*node)[n][0]);
        for (n=0; n < nn; n++)
          fprintf(fp,"%12.5e\n",(*node)[n][1]);
        for (n=0; n < nn; n++)
          fprintf(fp,"%12.5e\n",(*node)[n][2]);
      } else
      {
        strncpy(buff,"coordinates",sizeof(buff)-1);
        fwrite(buff,sizeof(char),80,fp);
        fwrite(&nn,sizeof(int),1,fp);
        for (n=1; n <= nn; n++)
          fwrite(&n,sizeof(int),1,fp);
        for (n=0; n < nn; n++)
        {
          fx = (float)(*node)[n][0];
          fwrite(&fx,sizeof(float),1,fp);
        }
        for (n=0; n < nn; n++)
        {
          fy = (float)(*node)[n][1];
          fwrite(&fy,sizeof(float),1,fp);
        }
        for (n=0; n < nn; n++)
        {
          fz = (float)(*node)[n][2];
          fwrite(&fz,sizeof(float),1,fp);
        }
      }

      if (ntet > 0)
      {
        if (ascii)
        {
          fprintf(fp,"tetra4\n");
          fprintf(fp,"%10d\n",ntet);
          for (n=0; n < ntet; n++)
            fprintf(fp,"%10d\n",n+1);
          for (n=0; n < ntet; n++)
          {
            for (j=0; j < 4; j++)
              fprintf(fp,"%10d",(*tet_n)[n][j]+1);
            fprintf(fp,"\n");
          }
        } else
        {
          strncpy(buff,"tetra4",sizeof(buff)-1);
          fwrite(buff,sizeof(char),80,fp);
          fwrite(&ntet,sizeof(int),1,fp);
          for (n=1; n <= ntet; n++)
            fwrite(&n,sizeof(int),1,fp);
          for (n=0; n < ntet; n++)
          {
            for (j=0; j < 4; j++)
            {
              i=(*tet_n)[n][j]+1;
              fwrite(&i,sizeof(int),1,fp);
            }
          }
        }
      }
      if (npyr > 0)
      {
        if (ascii)
        {
          fprintf(fp,"pyramid5\n");
          fprintf(fp,"%10d\n",npyr);
          for (n=0; n < npyr; n++)
            fprintf(fp,"%10d\n",n+1);
          for (n=0; n < npyr; n++)
          {
            for (j=0; j < 5; j++)
              fprintf(fp,"%10d",(*pyr_n)[n][j]+1);
            fprintf(fp,"\n");
          }
        } else
        {
          strncpy(buff,"pyramid5",sizeof(buff)-1);
          fwrite(buff,sizeof(char),80,fp);
          fwrite(&npyr,sizeof(int),1,fp);
          for (n=1; n <= npyr; n++)
            fwrite(&n,sizeof(int),1,fp);
          for (n=0; n < npyr; n++)
          {
            for (j=0; j < 5; j++)
            {
              i=(*pyr_n)[n][j]+1;
              fwrite(&i,sizeof(int),1,fp);
            }
          }
        }
      }
      if (npri > 0)
      {
        if (ascii)
        {
          fprintf(fp,"penta6\n");
          fprintf(fp,"%10d\n",npri);
          for (n=0; n < npri; n++)
            fprintf(fp,"%10d\n",n+1);
          for (n=0; n < npri; n++)
          {
            for (j=0; j < 6; j++)
              fprintf(fp,"%10d",(*pri_n)[n][j]+1);
            fprintf(fp,"\n");
          }
        } else
        {
          strncpy(buff,"penta6",sizeof(buff)-1);
          fwrite(buff,sizeof(char),80,fp);
          fwrite(&npri,sizeof(int),1,fp);
          for (n=1; n <= npri; n++)
            fwrite(&n,sizeof(int),1,fp);
          for (n=0; n < npri; n++)
          {
            for (j=0; j < 6; j++)
            {
              i=(*pri_n)[n][j]+1;
              fwrite(&i,sizeof(int),1,fp);
            }
          }
        }
      }
      if (nhex > 0)
      {
        if (ascii)
        {
          fprintf(fp,"hexa8\n");
          fprintf(fp,"%10d\n",nhex);
          for (n=0; n < nhex; n++)
            fprintf(fp,"%10d\n",n+1);
          for (n=0; n < nhex; n++)
          {
            for (j=0; j < 8; j++)
              fprintf(fp,"%10d",(*hex_n)[n][j]+1);
            fprintf(fp,"\n");
          }
        } else
        {
          strncpy(buff,"hexa8",sizeof(buff)-1);
          fwrite(buff,sizeof(char),80,fp);
          fwrite(&nhex,sizeof(int),1,fp);
          for (n=1; n <= nhex; n++)
            fwrite(&n,sizeof(int),1,fp);
          for (n=0; n < nhex; n++)
          {
            for (j=0; j < 8; j++)
            {
              i=(*hex_n)[n][j]+1;
              fwrite(&i,sizeof(int),1,fp);
            }
          }
        }
      }
      if (nply > 0)
      {
        if (ascii)
        {
          fprintf(fp,"nfaced\n");
          fprintf(fp,"%10d\n",nply);
          for (n=0; n < nply; n++)
            fprintf(fp,"%10d\n",n+1);
          for (n=0; n < nply; n++)
            fprintf(fp,"%10d\n",(*poly_n)[n][0][0]);
          for (n=0; n < nply; n++)
            for (j=1; j <= (*poly_n)[n][0][0]; j++)
              fprintf(fp,"%10d\n",(*poly_n)[n][j][0]);
          for (n=0; n < nply; n++)
          {
            for (j=1; j <= (*poly_n)[n][0][0]; j++)
            {
              for (k=1; k <= (*poly_n)[n][j][0]; k++)
                fprintf(fp,"%10d",(*poly_n)[n][j][k]+1);
              fprintf(fp,"\n");
            }
          }
        } else
        {
          strncpy(buff,"nfaced",sizeof(buff)-1);
          fwrite(buff,sizeof(char),80,fp);
          fwrite(&nply,sizeof(int),1,fp);
          for (n=1; n <= nply; n++)
            fwrite(&n,sizeof(int),1,fp);
          for (n=0; n < nply; n++)
            fwrite(&((*poly_n)[n][0][0]),sizeof(int),1,fp);
          for (n=0; n < nply; n++)
            for (j=1; j <= (*poly_n)[n][0][0]; j++)
              fwrite(&((*poly_n)[n][j][0]),sizeof(int),1,fp);
          for (n=0; n < nply; n++)
          {
            for (j=1; j <= (*poly_n)[n][0][0]; j++)
            {
              for (k=1; k <= (*poly_n)[n][j][0]; k++)
              {
                i=(*poly_n)[n][j][k]+1;
                fwrite(&i,sizeof(int),1,fp);
              }
            }
          }
        }
      }
      
      int *map = new int[nn];
      for (b=0; b < nb; b++)
      {
        if (ascii)
        {
          fprintf(fp,"part\n%10d\n",++part);
          fprintf(fp,"Boundary %d %s\n",b+1,(*b_name)[b]);
        } else
        {
          strncpy(buff,"part",sizeof(buff)-1);
          fwrite(buff,sizeof(char),80,fp);
          ++part;
          fwrite(&part,sizeof(int),1,fp);
          strncpy(buff," ",sizeof(buff)-1);
          sprintf(buff,"Boundary %d %s\n",b+1,(*b_name)[b]);
          fwrite(buff,sizeof(char),80,fp);
        }

        for (n=0; n < nn; n++)
          map[n] = 0;
        for (i=0; i < (*nt)[b]; i++)
          for (j=0; j < 3; j++)
            map[(*t_n)[b][i][j]] = 1;
        for (i=0; i < (*nq)[b]; i++)
          for (j=0; j < 4; j++)
            map[(*q_n)[b][i][j]] = 1;
        for (i=0; i < (*ngon)[b]; i++)
          for (j=1; j <= (*ngon_n)[b][i][0]; j++)
            map[(*ngon_n)[b][i][j]] = 1;
        i=0;
        for (n=0; n < nn; n++)
          if (map[n] > 0) map[n] = ++i;

        if (ascii)
        {
          fprintf(fp,"coordinates\n");
          fprintf(fp,"%10d\n",i);
          for (n=0; n < nn; n++)
            if (map[n] > 0)
              fprintf(fp,"%10d\n",n+1);
          for (n=0; n < nn; n++)
            if (map[n] > 0)
              fprintf(fp,"%12.5e\n",(*node)[n][0]);
          for (n=0; n < nn; n++)
            if (map[n] > 0)
              fprintf(fp,"%12.5e\n",(*node)[n][1]);
          for (n=0; n < nn; n++)
            if (map[n] > 0)
              fprintf(fp,"%12.5e\n",(*node)[n][2]);
        } else
        {
          strncpy(buff,"coordinates",sizeof(buff)-1);
          fwrite(buff,sizeof(char),80,fp);
          fwrite(&i,sizeof(int),1,fp);
          for (n=0; n < nn; n++)
            if (map[n] > 0)
              fwrite(&map[n],sizeof(int),1,fp);
          for (n=0; n < nn; n++)
          {
            if (map[n] > 0)
            {
              fx = (float)(*node)[n][0];
              fwrite(&fx,sizeof(float),1,fp);
            }
          }
          for (n=0; n < nn; n++)
          {
            if (map[n] > 0)
            {
              fy = (float)(*node)[n][1];
              fwrite(&fy,sizeof(float),1,fp);
            }
          }
          for (n=0; n < nn; n++)
          {
            if (map[n] > 0)
            {
              fz = (float)(*node)[n][2];
              fwrite(&fz,sizeof(float),1,fp);
            }
          }
        }

        if ((*nt)[b] > 0)
        {
          if (ascii)
          {
            fprintf(fp,"tria3\n");
            fprintf(fp,"%10d\n",(*nt)[b]);
            for (n=0; n < (*nt)[b]; n++)
              fprintf(fp,"%10d\n",n+1);
            for (n=0; n < (*nt)[b]; n++)
            {
              for (j=0; j < 3; j++)
                fprintf(fp,"%10d",map[(*t_n)[b][n][j]]);
              fprintf(fp,"\n");
            }
          } else
          {
            strncpy(buff,"tria3",sizeof(buff)-1);
            fwrite(buff,sizeof(char),80,fp);
            fwrite(&((*nt)[b]),sizeof(int),1,fp);
            for (n=1; n <= (*nt)[b]; n++)
              fwrite(&n,sizeof(int),1,fp);
            for (n=0; n < (*nt)[b]; n++)
              for (j=0; j < 3; j++)
                fwrite(&map[(*t_n)[b][n][j]],sizeof(int),1,fp);
          }
        }
        if ((*nq)[b] > 0)
        {
          if (ascii)
          {
            fprintf(fp,"quad4\n");
            fprintf(fp,"%10d\n",(*nq)[b]);
            for (n=0; n < (*nq)[b]; n++)
              fprintf(fp,"%10d\n",n+1);
            for (n=0; n < (*nq)[b]; n++)
            {
              for (j=0; j < 4; j++)
                fprintf(fp,"%10d",map[(*q_n)[b][n][j]]);
              fprintf(fp,"\n");
            }
          } else
          {
            strncpy(buff,"quad4",sizeof(buff)-1);
            fwrite(buff,sizeof(char),80,fp);
            fwrite(&((*nq)[b]),sizeof(int),1,fp);
            for (n=1; n <= (*nq)[b]; n++)
              fwrite(&n,sizeof(int),1,fp);
            for (n=0; n < (*nq)[b]; n++)
              for (j=0; j < 4; j++)
                fwrite(&map[(*q_n)[b][n][j]],sizeof(int),1,fp);
          }
        }
        if ((*ngon)[b] > 0)
        {
          if (ascii)
          {
            fprintf(fp,"nsided\n");
            fprintf(fp,"%10d\n",(*ngon)[b]);
            for (n=0; n < (*ngon)[b]; n++)
              fprintf(fp,"%10d\n",n+1);
            for (n=0; n < (*ngon)[b]; n++)
              fprintf(fp,"%10d\n",(*ngon_n)[b][n][0]);
            for (n=0; n < (*ngon)[b]; n++)
            {
              for (j=1; j <= (*ngon_n)[b][n][0]; j++)
                fprintf(fp,"%10d",map[(*ngon_n)[b][n][j]]);
              fprintf(fp,"\n");
            }
          } else
          {
            strncpy(buff,"nsided",sizeof(buff)-1);
            fwrite(buff,sizeof(char),80,fp);
            fwrite(&((*ngon)[b]),sizeof(int),1,fp);
            for (n=1; n <= (*ngon)[b]; n++)
              fwrite(&n,sizeof(int),1,fp);
            for (n=0; n < (*ngon)[b]; n++)
              fwrite(&((*ngon_n)[b][n][0]),sizeof(int),1,fp);
            for (n=0; n < (*ngon)[b]; n++)
              for (j=1; j <= (*ngon_n)[b][n][0]; j++)
                fwrite(&map[(*ngon_n)[b][n][j]],sizeof(int),1,fp);
          }
        }
      }
      delete[] map;
      fclose(fp);

      break;

  }

  return(0);

}

int Generic_mesh(char filename[], int mode, int &nn, Point **node, int &nb, char ***b_name, int **nt, int ****t_n, 
             int **nq, int ****q_n, int **ngon, int ****ngon_n, int &ntet, int ***tet_n, int &npyr, int ***pyr_n,
             int &npri, int ***pri_n, int &nhex, int ***hex_n, int &nply, int ****poly_n)
{
  FILE* fp;
  int const bdim = 80;
  char buff[bdim];
  int i, j, nc, n, eflag;
  int ibuf[10];
  double dx, dy, dz;

  eflag = 0;
  //
  // read/write Generic mesh file
  //
  switch(mode)
  {
    case -1: // Read ASCII
      fprintf(stdout,"\nFilename = <%s>",filename);
      // Open file for write
      if ((fp = fopen(filename,"r")) == 0)
      {
        fprintf(stdout,"\nGeneric_mesh: Error opening file <%s>.",filename);
        exit(0);
      }
      fgets(buff,bdim,fp);
      fgets(buff,bdim,fp);
      sscanf(buff,"%d",&nn);
      fprintf(stdout,"\nNumber of points in file = %d",nn);
      (*node) = (Point*)malloc(nn*sizeof(Point));
      for (n=0; n < nn; n++)
      {
        fgets(buff,bdim,fp);
        sscanf(buff,"%lf %lf %lf",&dx,&dy,&dz);
        (*node)[n] = Point(dx,dy,dz);
      }

      fgets(buff,bdim,fp);
      fgets(buff,bdim,fp);
      sscanf(buff,"%d %d %d %d",&ntet,&npyr,&npri,&nhex);
      fprintf(stdout,"\n # of tetrahedra = %d",ntet);
      fprintf(stdout,"\n # of pyramid    = %d",npyr);
      fprintf(stdout,"\n # of prism      = %d",npri);
      fprintf(stdout,"\n # of hexahedra  = %d",nhex);
      fflush(stdout);

      nply = 0;

      if (ntet > 0)
      {
        int nc;
        (*tet_n) = (int**)malloc(ntet*sizeof(int*));
        int numtets = ntet;
        ntet = 0;
        for (n=0; n < numtets; n++)
        {
          fgets(buff,bdim,fp);
          nc = sscanf(buff,"%d %d %d %d %d %d %d %d %d %d",
                 &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3], &ibuf[4],
                 &ibuf[5], &ibuf[6], &ibuf[7], &ibuf[8], &ibuf[9]);
          if (nc == 4)
          {
            (*tet_n)[ntet] = (int*)malloc(4*sizeof(int));
            for (i=0; i < nc; i++)
              (*tet_n)[n][i] = ibuf[i]-1;
            ntet++;
          } else if (nc == 10)
          {
            for (i=0; i < nc; i++)
              ibuf[i]--;
            // add polyhedra with 4 polygonal faces
            (*poly_n) = (int***)realloc((void*)(*poly_n),(nply+1)*sizeof(int**));
            (*poly_n)[nply] = (int**)malloc(5*sizeof(int*));
            (*poly_n)[nply][0] = (int*)malloc(sizeof(int));
            (*poly_n)[nply][0][0] = 4;
            for (i=1; i <= 4; i++)
            {
              (*poly_n)[nply][i] = (int*)malloc(7*sizeof(int));
              (*poly_n)[nply][i][0] = 6;
              switch(i)
              {
                case 1:
                  (*poly_n)[nply][i][1] = ibuf[0];
                  (*poly_n)[nply][i][2] = ibuf[6];
                  (*poly_n)[nply][i][3] = ibuf[2];
                  (*poly_n)[nply][i][4] = ibuf[5];
                  (*poly_n)[nply][i][5] = ibuf[1];
                  (*poly_n)[nply][i][6] = ibuf[4];
                  break;
                case 2:
                  (*poly_n)[nply][i][1] = ibuf[0];
                  (*poly_n)[nply][i][2] = ibuf[4];
                  (*poly_n)[nply][i][3] = ibuf[1];
                  (*poly_n)[nply][i][4] = ibuf[8];
                  (*poly_n)[nply][i][5] = ibuf[3];
                  (*poly_n)[nply][i][6] = ibuf[7];
                  break;
                case 3:
                  (*poly_n)[nply][i][1] = ibuf[1];
                  (*poly_n)[nply][i][2] = ibuf[5];
                  (*poly_n)[nply][i][3] = ibuf[2];
                  (*poly_n)[nply][i][4] = ibuf[9];
                  (*poly_n)[nply][i][5] = ibuf[3];
                  (*poly_n)[nply][i][6] = ibuf[8];
                  break;
                case 4:
                  (*poly_n)[nply][i][1] = ibuf[2];
                  (*poly_n)[nply][i][2] = ibuf[6];
                  (*poly_n)[nply][i][3] = ibuf[0];
                  (*poly_n)[nply][i][4] = ibuf[7];
                  (*poly_n)[nply][i][5] = ibuf[3];
                  (*poly_n)[nply][i][6] = ibuf[9];
                  break;
              }
            }
            nply++;
          } else
          {
            fprintf(stderr,"\nMESH3D tet format invalid.");
            fflush(stderr);
            exit(0);
          }
        }
        if (ntet < numtets)
        {
          (*tet_n) = (int**)realloc((void*)(*tet_n),ntet*sizeof(int*));
          fprintf(stdout,"\n # of tetrahedra = %d",ntet);
          fprintf(stdout,"\n # of polyhedra  = %d",nply);
        }
      }
      if (npyr > 0)
      {
        (*pyr_n) = (int**)malloc(npyr*sizeof(int*));
        for (n=0; n < npyr; n++)
        {
          (*pyr_n)[n] = (int*)malloc(5*sizeof(int));
          fgets(buff,bdim,fp);
          sscanf(buff,"%d %d %d %d %d", &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3], &ibuf[4]);
          (*pyr_n)[n][0] = ibuf[0]-1;
          (*pyr_n)[n][1] = ibuf[1]-1;
          (*pyr_n)[n][2] = ibuf[2]-1;
          (*pyr_n)[n][3] = ibuf[3]-1;
          (*pyr_n)[n][4] = ibuf[4]-1;
        }
      }
      if (npri > 0)
      {
        (*pri_n) = (int**)malloc(npri*sizeof(int*));
        for (n=0; n < npri; n++)
        {
          (*pri_n)[n] = (int*)malloc(6*sizeof(int));
          fgets(buff,bdim,fp);
          sscanf(buff,"%d %d %d %d %d %d", &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3], &ibuf[4], &ibuf[5]);
          (*pri_n)[n][0] = ibuf[0]-1;
          (*pri_n)[n][1] = ibuf[1]-1;
          (*pri_n)[n][2] = ibuf[2]-1;
          (*pri_n)[n][3] = ibuf[3]-1;
          (*pri_n)[n][4] = ibuf[4]-1;
          (*pri_n)[n][5] = ibuf[5]-1;
        }
      }
      if (nhex > 0)
      {
        (*hex_n) = (int**)malloc(nhex*sizeof(int*));
        for (n=0; n < nhex; n++)
        {
          (*hex_n)[n] = (int*)malloc(8*sizeof(int));
          fgets(buff,bdim,fp);
          sscanf(buff,"%d %d %d %d %d %d %d %d", &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3],
                                                 &ibuf[4], &ibuf[5], &ibuf[6], &ibuf[7]);
          (*hex_n)[n][0] = ibuf[0]-1;
          (*hex_n)[n][1] = ibuf[1]-1;
          (*hex_n)[n][2] = ibuf[2]-1;
          (*hex_n)[n][3] = ibuf[3]-1;
          (*hex_n)[n][4] = ibuf[4]-1;
          (*hex_n)[n][5] = ibuf[5]-1;
          (*hex_n)[n][6] = ibuf[6]-1;
          (*hex_n)[n][7] = ibuf[7]-1;
        }
      }

      fgets(buff,bdim,fp);
      fgets(buff,bdim,fp);
      sscanf(buff,"%d %d",&nb,&eflag);
      (*nt) = (int*)malloc(nb*sizeof(int));
      (*t_n) = (int***)malloc(nb*sizeof(int**));
      (*nq) = (int*)malloc(nb*sizeof(int));
      (*q_n) = (int***)malloc(nb*sizeof(int**));
      (*b_name) = (char**)malloc(nb*sizeof(char*));
      (*ngon) = (int*)malloc(nb*sizeof(int));
      (*ngon_n) = (int***)malloc(nb*sizeof(int**));
      fprintf(stdout,"\nNumber of boundaries = %d",nb);
      fflush(stdout);
      for (n=0; n < nb; n++)
      {
        (*nt)[n]=0;
        (*nq)[n]=0;
        (*t_n)[n]=0;
        (*q_n)[n]=0;
        (*ngon)[n]=0;
        (*ngon_n)[n]=0;
        (*b_name)[n] = (char*)malloc(33*sizeof(char));
        fgets(buff,bdim,fp);
        sscanf(buff,"%s",(*b_name)[n]);
        fprintf(stdout,"\nBoundary %d = %s",n+1,(*b_name)[n]);
        fflush(stdout);
        fgets(buff,bdim,fp); // boundary type
        fgets(buff,bdim,fp); // geometry index
        fgets(buff,bdim,fp);
        sscanf(buff,"%d %d",&((*nt)[n]),&((*nq)[n]));
        (*ngon)[n]=0;

        if ((*nt)[n] > 0)
        {
          int numtri = (*nt)[n];
          (*t_n)[n] =(int**)realloc((void*)(*t_n)[n],((*nt)[n]+1)*sizeof(int*));
          (*nt)[n]=0;
          for (i=0; i < numtri; i++)
          {
            fgets(buff,bdim,fp);
            nc = sscanf(buff,"%d %d %d %d %d %d %d",
                 &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3], &ibuf[4], &ibuf[5], &ibuf[6]);
            for (j=0; j < nc; j++)
              ibuf[j]--;
            if ((nc==3 && !eflag) || (nc==4 && eflag))  // store triangle
            {
              (*t_n)[n][(*nt)[n]] = (int*)malloc(3*sizeof(int));
              for (j=0; j < 3; j++)
                (*t_n)[n][(*nt)[n]][j] = ibuf[j];
              (*nt)[n]++;
            } else if ((nc==6 && !eflag) || (nc==7 && eflag)) // store 6 node polygon
            {
              (*ngon_n)[n] =(int**)realloc((void*)(*ngon_n)[n],((*ngon)[n]+1)*sizeof(int*));
              (*ngon_n)[n][(*ngon)[n]] = (int*)malloc(7*sizeof(int));
              (*ngon_n)[n][(*ngon)[n]][0] = 6;
              (*ngon_n)[n][(*ngon)[n]][1] = ibuf[0];
              (*ngon_n)[n][(*ngon)[n]][2] = ibuf[3];
              (*ngon_n)[n][(*ngon)[n]][3] = ibuf[1];
              (*ngon_n)[n][(*ngon)[n]][4] = ibuf[4];
              (*ngon_n)[n][(*ngon)[n]][5] = ibuf[2];
              (*ngon_n)[n][(*ngon)[n]][6] = ibuf[5];
              (*ngon)[n]++;
            }
          }
          if ((*nt)[n] < numtri)
            (*t_n)[n] =(int**)realloc((void*)(*t_n)[n],((*nt)[n]+1)*sizeof(int*));
        }

        if ((*nq)[n] > 0)
        {
          (*q_n)[n] =(int**)realloc((void*)(*q_n)[n],((*nq)[n]+1)*sizeof(int*));
          for (i=0; i < (*nq)[n]; i++)
          {
            (*q_n)[n][i] = (int*)malloc(4*sizeof(int));
            fgets(buff,bdim,fp);
            if (eflag)
              sscanf(buff,"%d %d %d %d %d", &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3], &ibuf[4]);
            else
              sscanf(buff,"%d %d %d %d", &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3]);
            (*q_n)[n][i][0] = ibuf[0]-1;
            (*q_n)[n][i][1] = ibuf[1]-1;
            (*q_n)[n][i][2] = ibuf[2]-1;
            (*q_n)[n][i][3] = ibuf[3]-1;
          }
        }
      }
      fclose(fp);
      break;
    case 1: // Write ASCII
      fprintf(stdout,"\nFilename = <%s>",filename);
      // Open file for write
      if ((fp = fopen(filename,"w")) == 0)
      {
        fprintf(stdout,"\nGeneric_mesh: Error opening file <%s>.",filename);
        exit(0);
      }

      fprintf(fp,"#Number of Verts\n");
      fprintf(fp,"%d\n",nn);
      for (n=0; n < nn; n++)
        fprintf(fp,"%22.15e %22.15e %22.15e\n",(*node)[n][0],(*node)[n][1],(*node)[n][2]);

      fprintf(fp,"#Number of Elements\n");
      fprintf(fp,"%d %d %d %d\n",ntet,npyr,npri,nhex);

      for (n=0; n < ntet; n++)
        fprintf(fp,"%d %d %d %d\n",(*tet_n)[n][0]+1,(*tet_n)[n][1]+1,(*tet_n)[n][2]+1,(*tet_n)[n][3]+1);
      for (n=0; n < npyr; n++)
        fprintf(fp,"%d %d %d %d %d\n",(*pyr_n)[n][0]+1,(*pyr_n)[n][1]+1,(*pyr_n)[n][2]+1,(*pyr_n)[n][3]+1,
                                      (*pyr_n)[n][4]+1);
      for (n=0; n < npri; n++)
        fprintf(fp,"%d %d %d %d %d %d\n",(*pri_n)[n][0]+1,(*pri_n)[n][1]+1,(*pri_n)[n][2]+1,(*pri_n)[n][3]+1,
                                         (*pri_n)[n][4]+1,(*pri_n)[n][5]+1);
      for (n=0; n < nhex; n++)
        fprintf(fp,"%d %d %d %d %d %d %d %d\n",(*hex_n)[n][0]+1,(*hex_n)[n][1]+1,
                                               (*hex_n)[n][2]+1,(*hex_n)[n][3]+1,
                                               (*hex_n)[n][4]+1,(*hex_n)[n][5]+1,
                                               (*hex_n)[n][6]+1,(*hex_n)[n][7]+1);

      fprintf(fp,"#Number of Boundaries, element flag\n");
      fprintf(fp,"%d %d\n",nb,eflag);
      for (n=0; n < nb; n++)
      {
        fprintf(fp,"%s\n",(*b_name)[n]);
        fprintf(fp,"%s_type\n",(*b_name)[n]);
        fprintf(fp,"%s_geom_index\n",(*b_name)[n]);
        fprintf(fp,"%d %d\n",(*nt)[n],(*nq)[n]);

        for (i=0; i < (*nt)[n]; i++)
          fprintf(fp,"%d %d %d\n",(*t_n)[n][i][0]+1,(*t_n)[n][i][1]+1,(*t_n)[n][i][2]+1);
        for (i=0; i < (*nq)[n]; i++)
          fprintf(fp,"%d %d %d %d\n",(*q_n)[n][i][0]+1,(*q_n)[n][i][1]+1,(*q_n)[n][i][2]+1,(*q_n)[n][i][3]+1);
      }
      fclose(fp);

      break;
    default:
      break;
  }

  return(0);

}

int Plotfile(char filename[], int mode, int &nn, Point **node, int &nb, char ***b_name, int **nt, int ****t_n, 
             int **nq, int ****q_n, int **ngon, int ****ngon_n, int &ntet, int ***tet_n, int &npyr, int ***pyr_n,
             int &npri, int ***pri_n, int &nhex, int ***hex_n, int &nply, int ****poly_n)
{
  FILE* fp;
  int const bdim = 80;
  char buff[bdim],buff1[bdim],buff2[bdim],buff3[bdim],buff4[bdim];
  int b, c, i, j, k, m, n, nc, nf, t, q, major, minor, ftype;
  bool grid_only;
  int nv, n0, n1, n2, n3;
  int csize, isize, fsize;
  int ibuf[257], fvtag;
  float *ftmp;
  double dx, dy, dz;

  csize = sizeof(char);
  isize = sizeof(int);
  fsize = sizeof(float);

  //
  // read/write Fieldview file
  //
  switch(mode)
  {
    case -2: // Read ASCII
      fprintf(stdout,"\nFilename = <%s>",filename);
      // Open file for write
      if ((fp = fopen(filename,"r")) == 0)
      {
        fprintf(stdout,"\nError opening file <%s>.",filename);
        exit(0);
      }
      //fgets(buff,bdim,fp);
      //sscanf(buff,"FIELDVIEW %d %d",&major,&minor);
      fscanf(fp,"%s %d %d",buff,&major,&minor);
      for (i=0; i < (int)strlen(buff); i++)
        buff[i] = toupper(buff[i]);

      if (strcmp(buff,"FIELDVIEW_GRIDS") == 0)
        grid_only = true;
      else
      {
        grid_only = false;

        if (strcmp(buff,"FIELDVIEW"))
        {
          fprintf(stdout,"\nPlotfile(): Wrong Moniker -> %s",buff);
          fclose(fp);
          exit(0);
        }
      }
      fprintf(stdout,"\nFieldview version number %d.%d",major,minor);

      if (!grid_only)
      {
        ftmp = (float*)malloc(4*fsize);
        fscanf(fp,"%s %f %f %f %f",buff,&(ftmp[0]),&(ftmp[1]),&(ftmp[2]),&(ftmp[3]));
        free(ftmp);
      }

      do
      {
        if (fgets(buff,bdim,fp) == NULL)
        {
          fprintf(stdout,"\nPlotfile(): GRIDS marker not found!\n");
          fclose(fp);
          exit(0);
        }
        for (i=0; i < (int)strlen(buff); i++)
          buff[i] = toupper(buff[i]);
      } while (strstr(buff,"GRIDS") == NULL);

      // Number of grids
      j=sscanf(buff,"%s %d",buff1,&n);
      if (j != 2)
      {
        fgets(buff,bdim,fp);
        sscanf(buff,"%d",&n);
      }
      if (n != 1)
      {
        fprintf(stdout,"\nFieldview file contains multiple grids!");
        exit(0);
      }

      // Boundary table
      fgets(buff,bdim,fp);
      j=sscanf(buff,"%s %s %d",buff1,buff2,&k);
      for (i=0; i < (int)strlen(buff1); i++)
        buff1[i] = toupper(buff1[i]);
      for (i=0; i < (int)strlen(buff2); i++)
        buff2[i] = toupper(buff2[i]);
      if (strstr(buff1,"BOUNDARY") == NULL || strstr(buff2,"TABLE") == NULL)
      {
        fprintf(stdout,"\nPlotfile(): BOUNDARY TABLE marker not found -> %s",buff);
        fclose(fp);
        exit(0);
      }
      if (j != 3)
      {
        fgets(buff,bdim,fp);
        sscanf(buff,"%d",&k);
      }
      nb = k;
      (*nt) = (int*)malloc(nb*sizeof(int));
      (*t_n) = (int***)malloc(nb*sizeof(int**));
      (*nq) = (int*)malloc(nb*sizeof(int));
      (*q_n) = (int***)malloc(nb*sizeof(int**));
      (*ngon) = (int*)malloc(nb*sizeof(int));
      (*ngon_n) = (int***)malloc(nb*sizeof(int**));
      (*b_name) = (char**)malloc(nb*sizeof(char*));
      //wall = new int[nb];
      //bdata = new int[nb];
      //clockness = new int[nb];
      fprintf(stdout,"\nNumber of boundaries = %d",nb);
      fflush(stdout);
      for (n=0; n < nb; n++)
      {
        (*nt)[n]=0;
        (*nq)[n]=0;
        (*ngon)[n]=0;
        (*t_n)[n]=0;
        (*q_n)[n]=0;
        (*ngon_n)[n]=0;
        fgets(buff,bdim,fp);
        //printf("\nbuffer = %s",buff);
        //fflush(stdout);
        j = sscanf(buff,"%s %s %s %s",buff1, buff2, buff3, buff4);
        if (j == 1)
        {
          fprintf(stdout,"\nBody %d is named %s",n+1,buff1);
          strcpy(buff,buff1);
        } else if (j == 2)
        {
          fprintf(stdout,"\nBody %d is named %s",n+1,buff2);
          strcpy(buff,buff2);
        } else if (j == 4)
        {
          fprintf(stdout,"\nBody %d is named %s",n+1,buff4);
          strcpy(buff,buff4);
        } else
        {
          fprintf(stdout,"\nError reading boundary table! j=%d\n",j);
          exit(0);
        }
        (*b_name)[n] = (char*)malloc(33*sizeof(char));
        buff[33] = '\0';
        strcpy((*b_name)[n],buff);
        fprintf(stdout,"\nBody %d saved named = %s",n+1,(*b_name)[n]);
        fflush(stdout);
      }

      if (!grid_only)
      {
        // Variable Names
        do
        {
          fgets(buff,bdim,fp);
          for (i=0; i < (int)strlen(buff); i++)
            buff[i] = toupper(buff[i]);
          j=sscanf(buff,"VARIABLE NAMES %d",&n);
        } while(strstr(buff,"VARIABLE NAMES") == NULL);
        if (j != 1)
        {
          fgets(buff,bdim,fp);
          sscanf(buff,"%d",&n);
        }

        for (i=0; i < n; i++)
          fgets(buff,bdim,fp);
          //fscanf(fp,"%s",buff);

        // Boundary Variable Names
        if (major > 2 || minor >= 5)
        {
          //fgets(buff,bdim,fp);
          //fgets(buff,bdim,fp);
          //sscanf(buff,"%d",&n);
          fscanf(fp,"%s %s %s %d",buff,buff2,buff3,&n);
          for (i=0; i < (int)strlen(buff); i++)
            buff[i] = toupper(buff[i]);
          for (i=0; i < (int)strlen(buff2); i++)
            buff2[i] = toupper(buff2[i]);
          for (i=0; i < (int)strlen(buff3); i++)
            buff3[i] = toupper(buff3[i]);
          if (strstr(buff,"BOUNDARY") == NULL ||
              strstr(buff2,"VARIABLE") == NULL ||
              strstr(buff3,"NAMES") == NULL)
          {
            fprintf(stdout,"\nPlotfile(): BOUNDARY VARIABLE NAMES marker not found -> %s",buff);
            fclose(fp);
            exit(0);
          }
          for (i=0; i < n; i++)
            fgets(buff,bdim,fp);
        }
      }

      // Nodes
      //fgets(buff,bdim,fp);
      //fgets(buff,bdim,fp);
      //sscanf(buff,"%d",&nn);
      fscanf(fp,"%s %d",buff,&k);
      nn = k;
      fprintf(stdout,"\n # of nodes = %d",nn);
      for (i=0; i < (int)strlen(buff); i++)
        buff[i] = toupper(buff[i]);
      if (strstr(buff,"NODES") == NULL)
      {
        fprintf(stdout,"\nPlotfile(): NODES marker not found -> %s",buff);
        fclose(fp);
        exit(0);
      }
      (*node) = (Point*)malloc(nn*sizeof(Point));
      for (n=0; n < nn; n++)
      {
        fscanf(fp,"%lf %lf %lf",&dx,&dy,&dz);
        (*node)[n] = Point(dx,dy,dz);
      }

      // Boundary Faces
      //fgets(buff,bdim,fp);
      //fgets(buff,bdim,fp);
      //sscanf(buff,"%d",&n);
      fscanf(fp,"%s %s %d",buff,buff2,&nf);
      for (i=0; i < (int)strlen(buff); i++)
        buff[i] = toupper(buff[i]);
      for (i=0; i < (int)strlen(buff2); i++)
        buff2[i] = toupper(buff2[i]);
      if (strstr(buff,"BOUNDARY") == NULL || strstr(buff2,"FACES") == NULL)
      {
        fprintf(stdout,"\nPlotfile(): BOUNDARY FACES marker not found -> %s",buff);
        fclose(fp);
        exit(0);
      }
      for (i=0; i < nf; i++)
      {
        fscanf(fp,"%d %d",&b,&j);
        if (b <= 0 || b > nb)
        {
          fprintf(stdout,"\nInvalid boundary number -> %d",b);
          exit(0);
        }
        b--;
        if (j==3)
        {
          fscanf(fp,"%d %d %d",&(ibuf[0]),
                               &(ibuf[1]),
                               &(ibuf[2]));
          (*t_n)[b] =(int**)realloc((void*)(*t_n)[b],((*nt)[b]+1)*sizeof(int*));
          (*t_n)[b][(*nt)[b]] = (int*)malloc(3*sizeof(int));
          (*t_n)[b][(*nt)[b]][0] = ibuf[0]-1;
          (*t_n)[b][(*nt)[b]][1] = ibuf[1]-1;
          (*t_n)[b][(*nt)[b]][2] = ibuf[2]-1;
          (*nt)[b]++;
        }
        if (j==4)
        {
          fscanf(fp,"%d %d %d %d",&(ibuf[0]),
                                  &(ibuf[1]),
                                  &(ibuf[2]),
                                  &(ibuf[3]));
          (*q_n)[b] =(int**)realloc((void*)(*q_n)[b],((*nq)[b]+1)*sizeof(int*));
          (*q_n)[b][(*nq)[b]] = (int*)malloc(4*sizeof(int));
          (*q_n)[b][(*nq)[b]][0] = ibuf[0]-1;
          (*q_n)[b][(*nq)[b]][1] = ibuf[1]-1;
          (*q_n)[b][(*nq)[b]][2] = ibuf[2]-1;
          (*q_n)[b][(*nq)[b]][3] = ibuf[3]-1;
          (*nq)[b]++;
        }
        if (j > 4)
        {
          (*ngon_n)[b] =(int**)realloc((void*)(*ngon_n)[b],((*ngon)[b]+1)*sizeof(int*));
          (*ngon_n)[b][(*ngon)[b]] = (int*)malloc(j*sizeof(int));
          (*ngon_n)[b][(*ngon)[b]][0] = j;
          for (k=1; k <= j; k++)
          {
            fscanf(fp,"%d",ibuf);
            (*ngon_n)[b][(*ngon)[b]][k] = ibuf[0]-1;
          }
          (*ngon)[b]++;
        }
      }
      b=0;
      while (b < nb)
      {
        if ((*nt)[b] > 0) fprintf(stdout,"\nBody %d has %d triangles.",b+1,(*nt)[b]);
        if ((*nq)[b] > 0) fprintf(stdout,"\nBody %d has %d quadrilaterals.",b+1,(*nq)[b]);
        if ((*ngon)[b] > 0) fprintf(stdout,"\nBody %d has %d polygons.",b+1,(*ngon)[b]);
        fflush(stdout);
        if ((*nt)[b] <= 0 && (*nq)[b] <= 0 && (*ngon)[b] <= 0)
        {
          fprintf(stdout,"\nNo boundary elements. Eliminating boundary!");
          fflush(stdout);
          for (i=b; i < nb-1; i++)
          {
            if ((*nt)[i+1] > 0)
            {
              (*t_n)[i] =(int**)realloc((void*)(*t_n)[i],((*nt)[i+1])*sizeof(int*));
              for (j=0; j < (*nt)[i+1]; j++)
              {
                (*t_n)[i][j] = (int*)malloc(3*sizeof(int));
                (*t_n)[i][j][0] = (*t_n)[i+1][j][0];
                (*t_n)[i][j][1] = (*t_n)[i+1][j][1];
                (*t_n)[i][j][2] = (*t_n)[i+1][j][2];
                free((*t_n)[i+1][j]);
              }
              free((*t_n)[i+1]);
              (*t_n)[i+1] = 0;
              (*nt)[i] = (*nt)[i+1];
              (*nt)[i+1] = 0;
            }
            if ((*nq)[i+1] > 0)
            {
              (*q_n)[i] =(int**)realloc((void*)(*q_n)[i],((*nq)[i+1])*sizeof(int*));
              for (j=0; j < (*nq)[i+1]; j++)
              {
                (*q_n)[i][j] = (int*)malloc(4*sizeof(int));
                (*q_n)[i][j][0] = (*q_n)[i+1][j][0];
                (*q_n)[i][j][1] = (*q_n)[i+1][j][1];
                (*q_n)[i][j][2] = (*q_n)[i+1][j][2];
                (*q_n)[i][j][3] = (*q_n)[i+1][j][3];
                free((*q_n)[i+1][j]);
              }
              free((*q_n)[i+1]);
              (*q_n)[i+1] = 0;
              (*nq)[i] = (*nq)[i+1];
              (*nq)[i+1] = 0;
            }
            if ((*ngon)[i+1] > 0)
            {
              (*ngon_n)[i] =(int**)realloc((void*)(*ngon_n)[i],((*ngon)[i+1])*sizeof(int*));
              for (j=0; j < (*ngon)[i+1]; j++)
              {
                (*ngon_n)[i][j] = (int*)malloc(((*ngon_n)[i+1][j][0]+1)*sizeof(int));
                (*ngon_n)[i][j][0] = (*ngon_n)[i+1][j][0];
                for (k=1; k <= (*ngon_n)[i+1][j][0]; k++)
                  (*ngon_n)[i][j][k] = (*ngon_n)[i+1][j][k];
                free((*ngon_n)[i+1][j]);
              }
              free((*ngon_n)[i+1]);
              (*ngon_n)[i+1] = 0;
              (*ngon)[i] = (*ngon)[i+1];
              (*ngon)[i+1] = 0;
            }
            for (j=0; j < 32; j++)
              (*b_name)[i][j]=(*b_name)[i+1][j];
          }
          nb--;
          free((*b_name)[nb]);
          (*nt) = (int*)realloc((void*)(*nt),nb*sizeof(int));
          (*t_n) = (int***)realloc((void*)(*t_n),nb*sizeof(int**));
          (*nq) = (int*)realloc((void*)(*nq),nb*sizeof(int));
          (*q_n) = (int***)realloc((void*)(*q_n),nb*sizeof(int**));
          (*ngon) = (int*)realloc((void*)(*ngon),nb*sizeof(int));
          (*ngon_n) = (int***)realloc((void*)(*ngon_n),nb*sizeof(int**));
          (*b_name) = (char**)realloc((void*)(*b_name),nb*sizeof(char*));
        } else
          b++;
      }

      int tetdim, pyrdim, pridim, hexdim, polydim;

      tetdim = pyrdim = pridim = hexdim = polydim = 0;
      ntet=npri=npyr=nhex=nply=0;
      // Elements
      //fgets(buff,bdim,fp);
      //fgets(buff,bdim,fp);
      //sscanf(buff,"%d %d",&i,&j);
      fscanf(fp,"%s",buff);
      for (i=0; i < (int)strlen(buff); i++)
        buff[i] = toupper(buff[i]);
      if (strstr(buff,"ELEMENTS") == NULL)
      {
        fprintf(stdout,"\nPlotfile(): ELEMENTS marker not found -> %s",buff);
        fclose(fp);
        exit(0);
      }
      fvtag=1;
      while (fvtag != 0)
      {
        if (fscanf(fp,"%s",buff) == EOF)
          fvtag = 0;
        else
        {
          for (i=0; i < (int)strlen(buff); i++)
            buff[i] = toupper(buff[i]);
          fvtag=-1;
          if (strstr(buff,"VARIABLES") != NULL)
            fvtag=0;
          else
            sscanf(buff,"%d",&fvtag);
        }
        switch (fvtag)
        {
          case 1:
            if (ntet+1 > tetdim)
            {
              (*tet_n) =(int**)realloc((void*)(*tet_n),(tetdim+CHUNK)*sizeof(int*));
              for (n=tetdim; n < tetdim+CHUNK; n++)
                (*tet_n)[n] = (int*)malloc(4*sizeof(int));
              tetdim += CHUNK;
            }
            fscanf(fp,"%d",&j);
            fscanf(fp,"%d",&ibuf[0]);
            fscanf(fp,"%d",&ibuf[1]);
            fscanf(fp,"%d",&ibuf[2]);
            fscanf(fp,"%d",&ibuf[3]);
            n0 = ibuf[0]-1;
            n1 = ibuf[2]-1;
            n2 = ibuf[1]-1;
            n3 = ibuf[3]-1;
            (*tet_n)[ntet][0] = ibuf[0]-1;
            (*tet_n)[ntet][1] = ibuf[2]-1;
            (*tet_n)[ntet][2] = ibuf[1]-1;
            (*tet_n)[ntet][3] = ibuf[3]-1;
            ntet++;
            break;
          case 2:
            if (nhex+1 > hexdim)
            {
              (*hex_n) =(int**)realloc((void*)(*hex_n),(hexdim+CHUNK)*sizeof(int*));
              for (n=hexdim; n < hexdim+CHUNK; n++)
                (*hex_n)[n] = (int*)malloc(8*sizeof(int));
              hexdim += CHUNK;
            }
            fscanf(fp,"%d",&j);
            fscanf(fp,"%d",&ibuf[0]);
            fscanf(fp,"%d",&ibuf[1]);
            fscanf(fp,"%d",&ibuf[2]);
            fscanf(fp,"%d",&ibuf[3]);
            fscanf(fp,"%d",&ibuf[4]);
            fscanf(fp,"%d",&ibuf[5]);
            fscanf(fp,"%d",&ibuf[6]);
            fscanf(fp,"%d",&ibuf[7]);
            (*hex_n)[nhex][0]= ibuf[0]-1;
            (*hex_n)[nhex][1]= ibuf[1]-1;
            (*hex_n)[nhex][2]= ibuf[3]-1;
            (*hex_n)[nhex][3]= ibuf[2]-1;
            (*hex_n)[nhex][4]= ibuf[4]-1;
            (*hex_n)[nhex][5]= ibuf[5]-1;
            (*hex_n)[nhex][6]= ibuf[7]-1;
            (*hex_n)[nhex][7]= ibuf[6]-1;
            nhex++;
            break;
          case 3:
            if (npri+1 > pridim)
            {
              (*pri_n) =(int**)realloc((void*)(*pri_n),(pridim+CHUNK)*sizeof(int*));
              for (n=pridim; n < pridim+CHUNK; n++)
                (*pri_n)[n] = (int*)malloc(6*sizeof(int));
              pridim += CHUNK;
            }
            fscanf(fp,"%d",&j);
            fscanf(fp,"%d",&ibuf[0]);
            fscanf(fp,"%d",&ibuf[1]);
            fscanf(fp,"%d",&ibuf[2]);
            fscanf(fp,"%d",&ibuf[3]);
            fscanf(fp,"%d",&ibuf[4]);
            fscanf(fp,"%d",&ibuf[5]);
            (*pri_n)[npri][0]= ibuf[0]-1;
            (*pri_n)[npri][1]= ibuf[3]-1;
            (*pri_n)[npri][2]= ibuf[5]-1;
            (*pri_n)[npri][3]= ibuf[1]-1;
            (*pri_n)[npri][4]= ibuf[2]-1;
            (*pri_n)[npri][5]= ibuf[4]-1;
            npri++;
            break;
          case 4:
            if (npyr+1 > pyrdim)
            {
              (*pyr_n) =(int**)realloc((void*)(*pyr_n),(pyrdim+CHUNK)*sizeof(int*));
              for (n=pyrdim; n < pyrdim+CHUNK; n++)
                (*pyr_n)[n] = (int*)malloc(5*sizeof(int));
              pyrdim += CHUNK;
            }
            fscanf(fp,"%d",&j);
            fscanf(fp,"%d",&ibuf[0]);
            fscanf(fp,"%d",&ibuf[1]);
            fscanf(fp,"%d",&ibuf[2]);
            fscanf(fp,"%d",&ibuf[3]);
            fscanf(fp,"%d",&ibuf[4]);
            (*pyr_n)[npyr][0]= ibuf[0]-1;
            (*pyr_n)[npyr][1]= ibuf[1]-1;
            (*pyr_n)[npyr][2]= ibuf[2]-1;
            (*pyr_n)[npyr][3]= ibuf[3]-1;
            (*pyr_n)[npyr][4]= ibuf[4]-1;
            npyr++;
            break;
          case 5:
            if (nply+1 > polydim)
            {
              (*poly_n) =(int***)realloc((void*)(*poly_n),(polydim+CHUNK)*sizeof(int**));
              polydim += CHUNK;
            }
            fscanf(fp,"%d",&j);
            fscanf(fp,"%d",&j);
            (*poly_n)[nply] = (int**)malloc((j+1)*sizeof(int*));
            (*poly_n)[nply][0] = (int*)malloc(sizeof(int));
            (*poly_n)[nply][0][0] = j;
            fscanf(fp,"%d",&j);
            fscanf(fp,"%d",&j);
            for (i=1; i <= (*poly_n)[nply][0][0]; i++)
            {
              fscanf(fp,"%d",&j);
              (*poly_n)[nply][i] = (int*)malloc((j+1)*sizeof(int));
              (*poly_n)[nply][i][0] = j;
              for (k=1; k <= (*poly_n)[nply][i][0]; k++)
              {
                fscanf(fp,"%d",&j);
                (*poly_n)[nply][i][k] = j-1;
              }
              fscanf(fp,"%d",&j);
            }
            nply++;
            break;
          default:
            break;
        }
      }

      fclose(fp);

      if (ntet < tetdim)
      {
        for (n=ntet; n < tetdim; n++)
          free((*tet_n)[n]);
        (*tet_n) =(int**)realloc((void*)(*tet_n),ntet*sizeof(int*));
      }
      if (npyr < pyrdim)
      {
        for (n=npyr; n < pyrdim; n++)
          free((*pyr_n)[n]);
        (*pyr_n) =(int**)realloc((void*)(*pyr_n),npyr*sizeof(int*));
      }
      if (npri < pridim)
      {
        for (n=npri; n < pridim; n++)
          free((*pri_n)[n]);
        (*pri_n) =(int**)realloc((void*)(*pri_n),npri*sizeof(int*));
      }
      if (nhex < hexdim)
      {
        for (n=nhex; n < hexdim; n++)
          free((*hex_n)[n]);
        (*hex_n) =(int**)realloc((void*)(*hex_n),nhex*sizeof(int*));
      }
      if (nply < polydim)
      {
        (*poly_n) =(int***)realloc((void*)(*poly_n),nply*sizeof(int**));
      }

      fprintf(stdout,"\n # of tetrahedra = %d",ntet);
      fprintf(stdout,"\n # of pyramid    = %d",npyr);
      fprintf(stdout,"\n # of prism      = %d",npri);
      fprintf(stdout,"\n # of hexahedra  = %d",nhex);
      fprintf(stdout,"\n # of polyhedra  = %d",nply);

      //delete wall;
      //delete bdata;
      //delete clockness;

      break;
    case -1: // Read Binary
      fprintf(stdout,"\nFilename = <%s>",filename);
      // Open file for write
      if ((fp = fopen(filename,"rb")) == 0)
      {
        fprintf(stdout,"\nError opening file <%s>.",filename);
        exit(0);
      }
      // Fieldview bit pattern
      fread(ibuf,isize,1,fp);
      if (ibuf[0] != FV_MAGIC)
      {
        fprintf(stdout,"\nPlotfile(): Wrong MAGIC number -> %d",ibuf[0]);
        fclose(fp);
        exit(0);
      }

      // name and version number
      fread(buff,csize,80,fp);
      if (strstr(buff,"FIELDVIEW") == NULL && strstr(buff,"FieldView") == NULL)
      {
        fprintf(stdout,"\nPlotfile(): Wrong Moniker -> %s",buff);
        fclose(fp);
        exit(0);
      }
      fread(ibuf,isize,2,fp);
      major= ibuf[0];
      minor= ibuf[1];
      fprintf(stdout,"\nFieldview version number %d.%d",major,minor);

      ftmp = (float*)malloc(4*fsize);

      // File type
      if (major > 2)
      {
        fread(&ftype,isize,1,fp);
        fread(ibuf,isize,1,fp);
        if (ftype == 2)
        {
          fprintf(stdout,"\nFieldview file contains only solution!");
          exit(0);
        }
      } else
        ftype = 3;

      if (ftype == 3)
      {
        fprintf(stdout,"\nFieldview file contains grids & solution!");
        //nc = 4;
        //constant = (double*)malloc(nc*sizeof(double));
        fread(ftmp,fsize,4,fp);
        //constant[0] = ftmp[0];
        //constant[1] = ftmp[1];
        //constant[2] = ftmp[2];
        //constant[3] = ftmp[3];
        //printf("\nFloating constants = %g, %g, %g, %g",constant[0],constant[1],constant[2],constant[3]);
        fflush(stdout);
      }

      // Number of grids
      fread(ibuf,isize,1,fp);
      if (ibuf[0] != 1)
      {
        fprintf(stdout,"\nFieldview file contains multiple grids!");
        exit(0);
      }

      // Number of boundary types
      fread(&k,isize,1,fp);
      nb = k;

      fprintf(stdout,"\nNumber of boundaries = %d",nb);
      fflush(stdout);

      (*b_name) = (char**)malloc(nb*sizeof(char*));
      // Boundary types
      for (n=0; n < nb; n++)
      {
        (*b_name)[n] = (char*)malloc(33*sizeof(char));
        fread(ibuf,isize,2,fp);
        fread(buff,csize,80,fp);
        sprintf((*b_name)[n],"%-32s",buff);
        fprintf(stdout,"\n Input boundary %d name = %s",n+1,buff);
        fprintf(stdout,"\nStored boundary %d name = %s",n+1,(*b_name)[n]);
      }
      fflush(stdout);

      if (ftype == 3)
      {
        // Number of variables
        fread(&nv,isize,1,fp);
        fprintf(stdout,"\nNumber of variables = %d",nv);
        //var_name = (char**)malloc(nv*sizeof(char*));
        for (i=0; i < nv; i++)
        {
          //var_name[i] = (char*)malloc(32*sizeof(char));
          fread(buff,csize,80,fp);
          //sprintf(var_name[i],"%-32s",buff);
          //printf("\nVariable %d, name = %s",i+1,var_name[i]);
        }

        // Number of boundary variables
        if (major > 2 || minor >= 5)
        {
          fread(ibuf,isize,1,fp);
          fprintf(stdout,"\nNumber of boundary variables = %d",ibuf[0]);
          for (i=0; i < ibuf[0]; i++)
          {
            fread(buff,csize,80,fp);
            fprintf(stdout,"\nBoundary variable %d, name = %s",i+1,buff);
          }
        }
      }

      // nodes
      fread(ibuf,isize,2,fp);
      if (ibuf[0] != FV_NODES)
      {
        fprintf(stdout,"\nPlotfile(): Wrong FV_NODES number -> %d",ibuf[0]);
        fclose(fp);
        exit(0);
      }
      nn = ibuf[1];
      fprintf(stdout,"\nNumber of nodes = %d",nn);
      fflush(stdout);
      (*node) = (Point*)malloc(nn*sizeof(Point));
      ftmp = (float*)realloc((void*)ftmp,nn*sizeof(float));
      // X coordinate
      fread(ftmp,fsize,nn,fp);
      for (n=0; n < nn; n++)
        (*node)[n][0] = ftmp[n];
      // Y coordinate
      fread(ftmp,fsize,nn,fp);
      for (n=0; n < nn; n++)
        (*node)[n][1] = ftmp[n];
      // Z coordinate
      fread(ftmp,fsize,nn,fp);
      for (n=0; n < nn; n++)
      {
        (*node)[n][2] = ftmp[n];
      }

      (*nt) = (int*)malloc(nb*sizeof(int));
      (*t_n) = (int***)malloc(nb*sizeof(int**));
      (*nq) = (int*)malloc(nb*sizeof(int));
      (*q_n) = (int***)malloc(nb*sizeof(int**));
      (*ngon) = (int*)malloc(nb*sizeof(int));
      (*ngon_n) = (int***)malloc(nb*sizeof(int**));
      for (b=0; b < nb; b++)
      {
        (*nt)[b] = (*nq)[b] = (*ngon)[b] = 0;
        (*t_n)[b] = 0;
        (*q_n)[b] = 0;
        (*ngon_n)[b] = 0;
      }

      //if (nv > 0)
      //{
      //  variable = (double**)malloc(nn*sizeof(double*));
      //  for (n=0; n < nn; n++)
      //    variable[n] = (double*)malloc(nv*sizeof(double));
      //}

      // Boundary section
      while ((fread(&fvtag,isize,1,fp)==1) && fvtag == FV_FACES)
      {
        fread(ibuf,isize,2,fp);
        b = ibuf[0];
        i = ibuf[1];
        fprintf(stdout,"\nBody %d, # of faces= %d",b,i);
        b--;
        for (j=0; j < i; j++)
        {
          fread(ibuf,isize,4,fp);
          if (ibuf[3] == 0)
          {
            (*t_n)[b] =(int**)realloc((void*)(*t_n)[b],((*nt)[b]+1)*sizeof(int*));
            (*t_n)[b][(*nt)[b]] = (int*)malloc(3*sizeof(int));
            (*t_n)[b][(*nt)[b]][0] = ibuf[0]-1;
            (*t_n)[b][(*nt)[b]][1] = ibuf[1]-1;
            (*t_n)[b][(*nt)[b]][2] = ibuf[2]-1;
            (*nt)[b]++;
          } else
          {
            (*q_n)[b] =(int**)realloc((void*)(*q_n)[b],((*nq)[b]+1)*sizeof(int*));
            (*q_n)[b][(*nq)[b]] = (int*)malloc(4*sizeof(int));
            (*q_n)[b][(*nq)[b]][0] = ibuf[0]-1;
            (*q_n)[b][(*nq)[b]][1] = ibuf[1]-1;
            (*q_n)[b][(*nq)[b]][2] = ibuf[2]-1;
            (*q_n)[b][(*nq)[b]][3] = ibuf[3]-1;
            (*nq)[b]++;
          }
        }
      }
      while (major > 2 && fvtag == FV_ARB_POLY_FACES)
      {
        fread(ibuf,isize,2,fp);
        b = ibuf[0];
        i = ibuf[1];
        fprintf(stdout,"\nBody %d, # of faces= %d",b,i);
        b--;
        for (j=0; j < i; j++)
        {
          fread(ibuf,isize,1,fp);
          k=ibuf[0];
          (*ngon_n)[b] = (int**)realloc((void*)(*ngon_n)[b],((*ngon)[b]+1)*sizeof(int*));
          (*ngon_n)[b][(*ngon)[b]] = (int*)malloc((k+1)*sizeof(int));
          (*ngon_n)[b][(*ngon)[b]][0] = k;
          fread(ibuf,isize,k,fp);
          for (m=1; m <= k; m++)
            (*ngon_n)[b][(*ngon)[b]][m] = ibuf[m-1]-1;
          (*ngon)[b]++;
        }
      }
      fflush(stdout);

      b=0;
      while (b < nb)
      {
        if ((*nt)[b] > 0) fprintf(stdout,"\nBody %d has %d triangles.",b+1,(*nt)[b]);
        if ((*nq)[b] > 0) fprintf(stdout,"\nBody %d has %d quadrilaterals.",b+1,(*nq)[b]);
        if ((*ngon)[b] > 0) fprintf(stdout,"\nBody %d has %d polygons.",b+1,(*ngon)[b]);
        fflush(stdout);
        if ((*nt)[b] <= 0 && (*nq)[b] <= 0 && (*ngon)[b] <= 0)
        {
          fprintf(stdout,"\nNo boundary elements. Eliminating boundary!");
          fflush(stdout);
          for (i=b; i < nb-1; i++)
          {
            if ((*nt)[i+1] > 0)
            {
              (*t_n)[i] =(int**)realloc((void*)(*t_n)[i],((*nt)[i+1])*sizeof(int*));
              for (j=0; j < (*nt)[i+1]; j++)
              {
                (*t_n)[i][j] = (int*)malloc(3*sizeof(int));
                (*t_n)[i][j][0] = (*t_n)[i+1][j][0];
                (*t_n)[i][j][1] = (*t_n)[i+1][j][1];
                (*t_n)[i][j][2] = (*t_n)[i+1][j][2];
                free((*t_n)[i+1][j]);
              }
              free((*t_n)[i+1]);
              (*t_n)[i+1] = 0;
              (*nt)[i] = (*nt)[i+1];
              (*nt)[i+1] = 0;
            }
            if ((*nq)[i+1] > 0)
            {
              (*q_n)[i] =(int**)realloc((void*)(*q_n)[i],((*nq)[i+1])*sizeof(int*));
              for (j=0; j < (*nq)[i+1]; j++)
              {
                (*q_n)[i][j] = (int*)malloc(4*sizeof(int));
                (*q_n)[i][j][0] = (*q_n)[i+1][j][0];
                (*q_n)[i][j][1] = (*q_n)[i+1][j][1];
                (*q_n)[i][j][2] = (*q_n)[i+1][j][2];
                (*q_n)[i][j][3] = (*q_n)[i+1][j][3];
                free((*q_n)[i+1][j]);
              }
              free((*q_n)[i+1]);
              (*q_n)[i+1] = 0;
              (*nq)[i] = (*nq)[i+1];
              (*nq)[i+1] = 0;
            }
            if ((*ngon)[i+1] > 0)
            {
              (*ngon_n)[i] =(int**)realloc((void*)(*ngon_n)[i],((*ngon)[i+1])*sizeof(int*));
              for (j=0; j < (*ngon)[i+1]; j++)
              {
                (*ngon_n)[i][j] = (int*)malloc(((*ngon_n)[i+1][j][0]+1)*sizeof(int));
                (*ngon_n)[i][j][0] = (*ngon_n)[i+1][j][0];
                for (k=1; k <= (*ngon_n)[i+1][j][0]; k++)
                  (*ngon_n)[i][j][k] = (*ngon_n)[i+1][j][k];
                free((*ngon_n)[i+1][j]);
              }
              free((*ngon_n)[i+1]);
              (*ngon_n)[i+1] = 0;
              (*ngon)[i] = (*ngon)[i+1];
              (*ngon)[i+1] = 0;
            }
            for (j=0; j < 32; j++)
              (*b_name)[i][j]=(*b_name)[i+1][j];
          }
          nb--;
          free((*b_name)[nb]);
          (*nt) = (int*)realloc((void*)(*nt),nb*sizeof(int));
          (*t_n) = (int***)realloc((void*)(*t_n),nb*sizeof(int**));
          (*nq) = (int*)realloc((void*)(*nq),nb*sizeof(int));
          (*q_n) = (int***)realloc((void*)(*q_n),nb*sizeof(int**));
          (*ngon) = (int*)realloc((void*)(*ngon),nb*sizeof(int));
          (*ngon_n) = (int***)realloc((void*)(*ngon_n),nb*sizeof(int**));
          (*b_name) = (char**)realloc((void*)(*b_name),nb*sizeof(char*));
        } else
          b++;
      }

      // Element section
      if (fvtag != FV_ELEMENTS)
      {
        fprintf(stdout,"\nRead_Plotfile(): Wrong FV_ELEMENTS number -> %d",ibuf[0]);
        fclose(fp);
        exit(0);
      }
      fread(ibuf,isize,4,fp);
      ntet = ibuf[0]; // # of tetrahedron
      nhex = ibuf[1]; // # of hexahedron
      npri = ibuf[2]; // # of prism
      npyr = ibuf[3]; // # of pyramid
      fprintf(stdout,"\n # of tetrahedra = %d",ntet);
      fprintf(stdout,"\n # of pyramid    = %d",npyr);
      fprintf(stdout,"\n # of prism      = %d",npri);
      fprintf(stdout,"\n # of hexahedra  = %d",nhex);

      if (ntet > 0)
      {
        (*tet_n) = (int**)malloc(ntet*sizeof(int*));
        for (n=0; n < ntet; n++)
          (*tet_n)[n] = (int*)malloc(4*sizeof(int));
      }
      if (npyr > 0)
      {
        (*pyr_n) = (int**)malloc(npyr*sizeof(int*));
        for (n=0; n < npyr; n++)
          (*pyr_n)[n] = (int*)malloc(5*sizeof(int));
      }
      if (npri > 0)
      {
        (*pri_n) = (int**)malloc(npri*sizeof(int*));
        for (n=0; n < npri; n++)
          (*pri_n)[n] = (int*)malloc(6*sizeof(int));
      }
      if (nhex > 0)
      {
        (*hex_n) = (int**)malloc(nhex*sizeof(int*));
        for (n=0; n < nhex; n++)
          (*hex_n)[n] = (int*)malloc(8*sizeof(int));
      }

      nc = ntet+nhex+npri+npyr;
      ntet = npyr = npri = nhex = 0;
      for (c=0; c < nc; c++)
      {
        fread(ibuf,isize,1,fp);
        switch(ibuf[0] >> 18)
        {
          case 1:
            fread(ibuf,isize,4,fp);
            (*tet_n)[ntet][0]= ibuf[0]-1;
            (*tet_n)[ntet][1]= ibuf[2]-1;
            (*tet_n)[ntet][2]= ibuf[1]-1;
            (*tet_n)[ntet][3]= ibuf[3]-1;
            ntet++;
            break;
          case 2:
            fread(ibuf,isize,5,fp);
            (*pyr_n)[npyr][0]= ibuf[0]-1;
            (*pyr_n)[npyr][1]= ibuf[1]-1;
            (*pyr_n)[npyr][2]= ibuf[2]-1;
            (*pyr_n)[npyr][3]= ibuf[3]-1;
            (*pyr_n)[npyr][4]= ibuf[4]-1;
            npyr++;
            break;
          case 3:
            fread(ibuf,isize,6,fp);
            (*pri_n)[npri][0]= ibuf[0]-1;
            (*pri_n)[npri][1]= ibuf[3]-1;
            (*pri_n)[npri][2]= ibuf[5]-1;
            (*pri_n)[npri][3]= ibuf[1]-1;
            (*pri_n)[npri][4]= ibuf[2]-1;
            (*pri_n)[npri][5]= ibuf[4]-1;
            npri++;
            break;
          case 4:
            fread(ibuf,isize,8,fp);
            (*hex_n)[nhex][0]= ibuf[0]-1;
            (*hex_n)[nhex][1]= ibuf[1]-1;
            (*hex_n)[nhex][2]= ibuf[3]-1;
            (*hex_n)[nhex][3]= ibuf[2]-1;
            (*hex_n)[nhex][4]= ibuf[4]-1;
            (*hex_n)[nhex][5]= ibuf[5]-1;
            (*hex_n)[nhex][6]= ibuf[7]-1;
            (*hex_n)[nhex][7]= ibuf[6]-1;
            nhex++;
            break;
          default:
            break;
        }
      }

      // Arbitrary element section
      if (major > 2 && (fread(&fvtag,isize,1,fp)==1) && fvtag == FV_ARB_POLY_ELEMENTS)
      {
        fread(ibuf,isize,1,fp);
        nply = ibuf[0]; // # of polyhedra
        fprintf(stdout,"\n # of polyhedra = %d",nply);
        (*poly_n) = (int***)malloc(nply*sizeof(int**));
        for (n=0; n < nply; n++)
        {
          fread(ibuf,isize,3,fp); // # faces, # nodes in element(ignored), center node(ignored)
          int nf = ibuf[0];
          (*poly_n)[n] = (int**)malloc((nf+1)*sizeof(int*));
          (*poly_n)[n][0] = (int*)malloc(sizeof(int));
          (*poly_n)[n][0][0] = nf;
          for (i=1; i <= nf; i++)
          {
            fread(ibuf,isize,1,fp); // wall value
            fread(ibuf,isize,1,fp); // number of face nodes
            (*poly_n)[n][i] = (int*)malloc((ibuf[0]+1)*sizeof(int));
            (*poly_n)[n][i][0] = ibuf[0];
            for (j=1; j <= ibuf[0]; j++)
            {
              fread(ibuf,isize,1,fp); // node index
              (*poly_n)[n][i][j] = ibuf[0]-1;
            }
            fread(ibuf,isize,1,fp); // number of hanging nodes
            k=ibuf[0];
            for (j=0; j < k; j++)
              fread(ibuf,isize,1,fp); // hanging nodes(ignored)
          }
        }
      }

      // Variable section
      //if (nv > 0)
      //{
      //  fread(ibuf,isize,1,fp);
      //  if (ibuf[0] != FV_VARIABLES)
      //  {
      //    fprintf(stdout,"\nRead_Plotfile(): Wrong FV_VARIABLES number -> %d",ibuf[0]);
      //    fclose(fp);
      //    exit(0);
      //  }

      //  for (j=0; j < nv; j++)
      //  {
      //    fread(ftmp,fsize,nn,fp);
      //    for (n=0; n < nn; n++)
      //      variable[n][j] = ftmp[n];
      //  }
      //}

      free(ftmp);

      fclose(fp);

      break;
    case 1: // Write Binary
      fprintf(stdout,"\nFilename = <%s>",filename);
      // Open file for write
      if ((fp = fopen(filename,"wb")) == 0)
      {
        fprintf(stdout,"\nError opening file <%s>.",filename);
        exit(0);
      }

      // Fieldview bit pattern
      ibuf[0] = FV_MAGIC;
      fwrite(ibuf,isize,1,fp);

      // name and version number
      sprintf(buff,"FIELDVIEW");
      fwrite(buff,csize,80,fp);
      ibuf[0] = 3;
      ibuf[1] = 0;
      fwrite(ibuf,isize,2,fp);
      ibuf[0] = FV_GRIDS_FILE;
      fwrite(ibuf,isize,1,fp);

      //if (nvars > 0 || nv > 0)
      //{
      //  ibuf[0] = FV_COMBINED_FILE;
      //  fwrite(ibuf,isize,1,fp);
      //}

      ibuf[0] = 0;
      fwrite(ibuf,isize,1,fp);

      //if (nvars > 0)
      //{
      //  ftmp = (float*)malloc(4*fsize);
      //  ftmp[0] = 1.0;
      //  ftmp[1] = 0.0;
      //  ftmp[2] = 0.0;
      //  ftmp[3] = 0.0;
      //  fwrite(ftmp,fsize,4,fp);
      //  free(ftmp);
      //} else if (nc > 0)
      //{
      //  ftmp = (float*)malloc(nc*fsize);
      //  for (n=0; n < nc; n++)
      //    ftmp[n] = constant[n];
      //  fwrite(ftmp,fsize,4,fp);
      //  free(ftmp);
      //}

      // Number of grids
      ibuf[0] = 1;
      fwrite(ibuf,isize,1,fp);

      // Number of boundary types
      ibuf[0] = nb;
      fwrite(ibuf,isize,1,fp);

      // Boundary types
      for (n=0; n < nb; n++)
      {
        for (i=0; i < 80; i++)
          buff[i] = ' ';
        sprintf(buff,"%s",(*b_name)[n]);
        ibuf[0] = 0;
        ibuf[1] = 1;
        fwrite(ibuf,isize,2,fp);
        fwrite(buff,csize,80,fp);
      }

      //if (nvars > 0)
      //{
      //  ibuf[0] = nvars;
      //  fwrite(ibuf,isize,1,fp);
      //  for (n=0; n < nvars; n++)
      //  {
      //    for (i=0; i < 80; i++)
      //      buff[i] = ' ';
      //    sprintf(buff,"%s",vname[n]);
      //    fwrite(buff,csize,80,fp);
      //  }
      //  ibuf[0] = 0;
      //  fwrite(ibuf,isize,1,fp);
      //} else if (nv > 0)
      //{
      //  ibuf[0] = nv;
      //  fwrite(ibuf,isize,1,fp);
      //  for (n=0; n < nv; n++)
      //  {
      //    for (i=0; i < 80; i++)
      //      buff[i] = ' ';
      //    sprintf(buff,"%s",var_name[n]);
      //    fwrite(buff,csize,80,fp);
      //  }
      //  ibuf[0] = 0;
      //  fwrite(ibuf,isize,1,fp);
      //}

      ibuf[0] = FV_NODES;
      ibuf[1] = nn;
      fwrite(ibuf,isize,2,fp);
      ftmp = (float*)malloc(nn*fsize);
      // X coordinate
      for (n=0; n < nn; n++)
        ftmp[n] = (float)((*node)[n][0]);
      fwrite(ftmp,fsize,nn,fp);
      // Y coordinate
      for (n=0; n < nn; n++)
        ftmp[n] = (float)((*node)[n][1]);
      fwrite(ftmp,fsize,nn,fp);
      // Z coordinate
      for (n=0; n < nn; n++)
        ftmp[n] = (float)((*node)[n][2]);
      fwrite(ftmp,fsize,nn,fp);

      free(ftmp);

      for (b=0; b < nb; b++)
      {
        ibuf[0] = FV_FACES;
        ibuf[1] = b+1;
        ibuf[2] = (*nt)[b]+(*nq)[b];
        fwrite(ibuf,isize,3,fp);
        for (t=0; t < (*nt)[b]; t++)
        {
          ibuf[0] = (*t_n)[b][t][0]+1;
          ibuf[1] = (*t_n)[b][t][1]+1;
          ibuf[2] = (*t_n)[b][t][2]+1;
          ibuf[3] = 0;
          fwrite(ibuf,isize,4,fp);
        }
        for (q=0; q < (*nq)[b]; q++)
        {
          ibuf[0] = (*q_n)[b][q][0]+1;
          ibuf[1] = (*q_n)[b][q][1]+1;
          ibuf[2] = (*q_n)[b][q][2]+1;
          ibuf[3] = (*q_n)[b][q][3]+1;
          fwrite(ibuf,isize,4,fp);
        }
      }
      for (b=0; b < nb; b++)
      {
        if ((*ngon)[b] == 0)
          continue;
        ibuf[0] = FV_ARB_POLY_FACES;
        ibuf[1] = b+1;
        ibuf[2] = (*ngon)[b];
        fwrite(ibuf,isize,3,fp);
        for (i=0; i < (*ngon)[b]; i++)
        {
          ibuf[0] = (*ngon_n)[b][i][0];
          for (j=1; j <= (*ngon_n)[b][i][0]; j++)
            ibuf[j] = (*ngon_n)[b][i][j]+1;
          fwrite(ibuf,isize,((*ngon_n)[b][i][0]+1),fp);
        }
      }

      // Element section
      int walls[6];
      walls[0] = NOT_A_WALL;
      walls[1] = NOT_A_WALL;
      walls[2] = NOT_A_WALL;
      walls[3] = NOT_A_WALL;
      walls[4] = NOT_A_WALL;
      walls[5] = NOT_A_WALL;
      ibuf[0] = FV_ELEMENTS;
      ibuf[1] = ntet; // # of tetrahedron
      ibuf[2] = nhex; // # of hexahedron
      ibuf[3] = npri; // # of prism
      ibuf[4] = npyr; // # of pyramid
      fwrite(ibuf,isize,5,fp);

      for (c=0; c < ntet; c++)
      {
        ibuf[0]=fv_encode_elem_header(FV_TET_ELEM_ID,walls);
        fwrite(ibuf,isize,1,fp);
        ibuf[0]=(*tet_n)[c][0]+1;
        ibuf[1]=(*tet_n)[c][2]+1;
        ibuf[2]=(*tet_n)[c][1]+1;
        ibuf[3]=(*tet_n)[c][3]+1;
        fwrite(ibuf,isize,4,fp);
      }

      for (c=0; c < npyr; c++)
      {
        ibuf[0]=fv_encode_elem_header(FV_PYRA_ELEM_ID,walls);
        fwrite(ibuf,isize,1,fp);
        ibuf[0]=(*pyr_n)[c][0]+1;
        ibuf[1]=(*pyr_n)[c][1]+1;
        ibuf[2]=(*pyr_n)[c][2]+1;
        ibuf[3]=(*pyr_n)[c][3]+1;
        ibuf[4]=(*pyr_n)[c][4]+1;
        fwrite(ibuf,isize,5,fp);
      }

      for (c=0; c < npri; c++)
      {
        ibuf[0]=fv_encode_elem_header(FV_PRISM_ELEM_ID,walls);
        fwrite(ibuf,isize,1,fp);
        ibuf[0]=(*pri_n)[c][0]+1;
        ibuf[1]=(*pri_n)[c][3]+1;
        ibuf[2]=(*pri_n)[c][4]+1;
        ibuf[3]=(*pri_n)[c][1]+1;
        ibuf[4]=(*pri_n)[c][5]+1;
        ibuf[5]=(*pri_n)[c][2]+1;
        fwrite(ibuf,isize,6,fp);
      }

      for (c=0; c < nhex; c++)
      {
        ibuf[0]=fv_encode_elem_header(FV_HEX_ELEM_ID,walls);
        fwrite(ibuf,isize,1,fp);
        ibuf[0]=(*hex_n)[c][0]+1;
        ibuf[1]=(*hex_n)[c][1]+1;
        ibuf[2]=(*hex_n)[c][3]+1;
        ibuf[3]=(*hex_n)[c][2]+1;
        ibuf[4]=(*hex_n)[c][4]+1;
        ibuf[5]=(*hex_n)[c][5]+1;
        ibuf[6]=(*hex_n)[c][7]+1;
        ibuf[7]=(*hex_n)[c][6]+1;
        fwrite(ibuf,isize,8,fp);
      }

      if (nply > 0)
      {
        List nlist;
        ibuf[0] = FV_ARB_POLY_ELEMENTS;
        ibuf[1] = nply; // # of tetrahedron
        fwrite(ibuf,isize,2,fp);
        for (c=0; c < nply; c++)
        {
          ibuf[0] = (*poly_n)[c][0][0];  // number of faces
          nlist.Redimension(0);
          for (i=1; i <= (*poly_n)[c][0][0]; i++)
            for (j=1; j <= (*poly_n)[c][i][0]; j++)
              nlist.Check_List((*poly_n)[c][i][j]);
          ibuf[1] = nlist.max; // number of nodes in element
          ibuf[2] = -1;  // center node
          fwrite(ibuf,isize,3,fp);
          for (i=1; i <= (*poly_n)[c][0][0]; i++)
          {
            ibuf[0] = NOT_A_WALL;
            fwrite(ibuf,isize,1,fp);  // wall value
            ibuf[0] = (*poly_n)[c][i][0];  // number of face nodes
            for (j=1; j <= (*poly_n)[c][i][0]; j++)
              ibuf[j] = (*poly_n)[c][i][j]+1;  // face nodes
            fwrite(ibuf,isize,((*poly_n)[c][i][0]+1),fp);
            ibuf[0] = 0;
            fwrite(ibuf,isize,1,fp); // number of hanging nodes
          }
        }
      }

      //if (nvars > 0)
      //{
      //  ibuf[0] = FV_VARIABLES;
      //  fwrite(ibuf,isize,1,fp);
      //  for (i=0; i < nvars; i++)
      //    fwrite(var[i],fsize,nn,fp);
      //} else if (nv > 0)
      //{
      //  ftmp = (float*)malloc(nn*fsize);
      //  ibuf[0] = FV_VARIABLES;
      //  fwrite(ibuf,isize,1,fp);
      //  for (i=0; i < nv; i++)
      //  {
      //    for (n=0; n < nn; n++)
      //      ftmp[n] = variable[n][i];
      //    fwrite(ftmp,fsize,nn,fp);
      //  }
      //  free(ftmp);
      //}

      fclose(fp);

      break;
    case 2: // Write ASCII
      fprintf(stdout,"\nFilename = <%s>",filename);
      // Open file for write
      if ((fp = fopen(filename,"w")) == 0)
      {
        fprintf(stdout,"\nError opening file <%s>.",filename);
        exit(0);
      }
      grid_only = true;
      // name and version number
      if (nply > 0)
        fprintf(fp,"FieldView_Grids 3 0\n");
      else
      {
        fprintf(fp,"FieldView 2 5\n");
        grid_only = false;
      }

      if (!grid_only)
      {
        fprintf(fp,"Constants\n");
        ftmp = (float*)malloc(4*fsize);
        //if (nc == 4)
        //{
        //  ftmp[0] = constant[0]; // Time
        //  ftmp[1] = constant[1]; // Mach
        //  ftmp[2] = constant[2]; // Alpha
        //  ftmp[3] = constant[3]; // Reynold's number
        //} else
        {
          ftmp[0] = 0.0; // Time
          ftmp[1] = 0.0; // Mach
          ftmp[2] = 0.0; // Alpha
          ftmp[3] = 0.0; // Reynold's number
        }
        fprintf(fp,"%g %g %g %g\n",ftmp[0],ftmp[1],ftmp[2],ftmp[3]);
        free(ftmp);
      }

      // Number of grids
      fprintf(fp,"Grids\n");
      fprintf(fp,"1\n");

      // Number of boundary types
      fprintf(fp,"Boundary Table\n");
      fprintf(fp,"%d\n",nb);
      fprintf(stdout,"\nPlotfile(): Number of boundaries = %d",nb);

      // Boundary types
      for (n=0; n < nb; n++)
      {
        for (i=0; i < 80; i++)
          buff[i] = ' ';
        sprintf(buff,"%s",(*b_name)[n]);
        fprintf(fp,"1 0 1 %s\n",buff);
      }

      // Number of variables
      if (!grid_only)
      {
        fprintf(fp,"Variable Names\n");
        //if (nvars > 0)
        //{
        //  fprintf(fp,"%d\n",nvars);
        //  for (n=0; n < nvars; n++)
        //    fprintf(fp,"%s\n",vname[n]);
        //} else
        //{
        //  fprintf(fp,"%d\n",nv);
        //  for (n=0; n < nv; n++)
        //    fprintf(fp,"%s\n",var_name[n]);
        //}
        nv = 0;
        fprintf(fp,"%d\n",nv);
        fprintf(fp,"Boundary Variable Names\n");
        fprintf(fp,"%d\n",nv);
      }

      // Node information  
      fprintf(fp,"Nodes\n");
      fprintf(fp,"%d\n",nn);
      fprintf(stdout,"\nPlotfile(): Number of nodes = %d",nn);
      for (n=0; n < nn; n++)
        fprintf(fp,"%22.15e %22.15e %22.15e\n",(*node)[n][0],(*node)[n][1],(*node)[n][2]);

      // Boundary section
      n=0;
      for (b=0; b < nb; b++)
        n += (*nt)[b]+(*nq)[b]+(*ngon)[b];
      fprintf(fp,"Boundary Faces\n");
      fprintf(fp,"%d\n",n);
      fprintf(stdout,"\nPlotfile(): Number of boundary faces = %d",n);
      for (b=0; b < nb; b++)
      {
        for (t=0; t < (*nt)[b]; t++)
        {
          fprintf(fp,"%d 3 %d %d %d\n",b+1,(*t_n)[b][t][0]+1,
                                           (*t_n)[b][t][1]+1,
                                           (*t_n)[b][t][2]+1);
        }
        for (q=0; q < (*nq)[b]; q++)
        {
          fprintf(fp,"%d 4 %d %d %d %d\n",b+1,(*q_n)[b][q][0]+1,
                                              (*q_n)[b][q][1]+1,
                                              (*q_n)[b][q][2]+1,
                                              (*q_n)[b][q][3]+1);
        }
        for (i=0; i < (*ngon)[b]; i++)
        {
          fprintf(fp,"%d %d",b+1,(*ngon_n)[b][i][0]);
          for (j=1; j <= (*ngon_n)[b][i][0]; j++)
            fprintf(fp," %d",(*ngon_n)[b][i][j]+1);
          fprintf(fp,"\n");
        }
      }

      // Element section
      fprintf(fp,"Elements\n");
      for (c=0; c < ntet; c++)
      {
        fprintf(fp,"1 1\n%d %d %d %d\n",(*tet_n)[c][0]+1,
                                        (*tet_n)[c][2]+1,
                                        (*tet_n)[c][1]+1,
                                        (*tet_n)[c][3]+1);
      }

      for (c=0; c < npyr; c++)
      {
        fprintf(fp,"4 1\n%d %d %d %d %d\n",(*pyr_n)[c][0]+1,
                                           (*pyr_n)[c][1]+1,
                                           (*pyr_n)[c][2]+1,
                                           (*pyr_n)[c][3]+1,
                                           (*pyr_n)[c][4]+1);
      }

      for (c=0; c < npri; c++)
      {
        fprintf(fp,"3 1\n%d %d %d %d %d %d\n",(*pri_n)[c][0]+1,
                                              (*pri_n)[c][3]+1,
                                              (*pri_n)[c][4]+1,
                                              (*pri_n)[c][1]+1,
                                              (*pri_n)[c][5]+1,
                                              (*pri_n)[c][2]+1);
      }

      for (c=0; c < nhex; c++)
      {
        fprintf(fp,"2 1\n%d %d %d %d %d %d %d %d\n",(*hex_n)[c][0]+1,
                                                    (*hex_n)[c][1]+1,
                                                    (*hex_n)[c][3]+1,
                                                    (*hex_n)[c][2]+1,
                                                    (*hex_n)[c][4]+1,
                                                    (*hex_n)[c][5]+1,
                                                    (*hex_n)[c][7]+1,
                                                    (*hex_n)[c][6]+1);
      }

      for (c=0; c < nply; c++)
      {
        List nlist;
        nlist.Redimension(0);
        for (i=1; i <= (*poly_n)[c][0][0]; i++)
          for (j=1; j <= (*poly_n)[c][i][0]; j++)
            nlist.Check_List((*poly_n)[c][i][j]);
        fprintf(fp,"5 1\n");
        fprintf(fp,"%d %d -1\n",(*poly_n)[c][0][0],nlist.max);
        for (i=1; i <= (*poly_n)[c][0][0]; i++)
        {
          fprintf(fp,"%d ",(*poly_n)[c][i][0]);
          for (j=1; j <= (*poly_n)[c][i][0]; j++)
            fprintf(fp,"%d ",(*poly_n)[c][i][j]+1);
          fprintf(fp,"0\n");
        }
      }

      // Variable section
      if (!grid_only)
        fprintf(fp,"Variables\n");

      //if (nvars > 0)
      //{
      //  for (i=0; i < nvars; i++)
      //    for (n=0; n < nn; n++)
      //      fprintf(fp,"%22.15e\n",var[i][n]);
      //} else
      //{
      //  for (i=0; i < nv; i++)
      //    for (n=0; n < nn; n++)
      //      fprintf(fp,"%22.15e\n",variable[n][i]);
      //}

      // Boundary variable section
      if (!grid_only)
        fprintf(fp,"Boundary Variables\n");

      fclose(fp);

      break;
    default:
      break;
  }

  return(0);

}

int Read_Mesh(char sname[], char **&b_name, int &nn, double *&x, double *&y, double *&z,
              int &nb, int *&nt, int *&nq, int *&ngon, int ***&t_n, int ***&q_n, int ***&ngon_n,
              int &ntet, int **&tet_n, int &npyr, int **&pyr_n,
              int &npri, int **&pri_n, int &nhex, int **&hex_n, int &nply, int ***&poly_n)
{
  int nvol;
  int *tet_vc, *pyr_vc, *pri_vc, *hex_vc, *poly_vc;
  char **vol_name;
  int error = 0;

  nvol = 0;
  tet_vc = 0;
  pyr_vc = 0;
  pri_vc = 0;
  hex_vc = 0;
  poly_vc = 0;
  vol_name = 0;

  error = Read_Mesh(sname,b_name,nn,x,y,z,nb,nt,nq,ngon,t_n,q_n,ngon_n,ntet,tet_n,npyr,pyr_n,npri,pri_n,nhex,hex_n,nply,poly_n,
            nvol,vol_name,tet_vc,pyr_vc,pri_vc,hex_vc,poly_vc);

  return error;
}

int Read_Mesh(char sname[], char **&b_name, int &nn, double *&x, double *&y, double *&z,
              int &nb, int *&nt, int *&nq, int *&ngon, int ***&t_n, int ***&q_n, int ***&ngon_n,
              int &ntet, int **&tet_n, int &npyr, int **&pyr_n,
              int &npri, int **&pri_n, int &nhex, int **&hex_n, int &nply, int ***&poly_n,
              int &nvol, char **&vol_name, int *&tet_vc, int *&pyr_vc, int *&pri_vc, int *&hex_vc, int *&poly_vc)
{
  int b, n;
  Point *node;
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

  if (strstr(sname,".sg") != NULL)
  {
    // read serial mesh file
    if ((error = SimGrid_read(sname, nn, &node, nb, &b_name, nvol, &vol_name, &nt, &t_n, &nq, &q_n, &ngon, &ngon_n,
         ntet, &tet_n, &tet_vc, npyr, &pyr_n, &pyr_vc, npri, &pri_n, &pri_vc, nhex, &hex_n, &hex_vc, nply, &poly_n, &poly_vc)) != 0)
    {
      fprintf(stderr,"\nError reading serial SG file <%s>.",sname);
      fflush(stderr);
      exit(0);
    }
  } else if (strstr(sname,".cgns") != NULL)
  {
    // read serial mesh file
#ifdef HAVE_CGNS
    if ((error = CGNS_read(sname, nn, &node, nb, &b_name, &nt, &t_n, &nq, &q_n, &ngon, &ngon_n,
                ntet, &tet_n, npyr, &pyr_n, npri, &pri_n, nhex, &hex_n, nply, &poly_n)) != 0)
    {
      fprintf(stderr,"\nError reading serial CGNS file <%s>.",sname);
      fflush(stderr);
      exit(0);
    }
#else
      fprintf(stderr,"\nCGNS not enabled!");
      fflush(stderr);
      exit(0);
#endif
  } else if (strstr(sname,".uns") != NULL || strstr(sname,".fv") != NULL)
  {
    if ((error = Plotfile(sname, -1, nn, &node, nb, &b_name, &nt, &t_n, &nq, &q_n, &ngon, &ngon_n,
                ntet, &tet_n, npyr, &pyr_n, npri, &pri_n, nhex, &hex_n, nply, &poly_n)) != 0)
    {
      fprintf(stderr,"\nError reading binary Fieldview file <%s>.",sname);
      fflush(stderr);
      exit(0);
    }
  } else if (strstr(sname,".crunch") != NULL || strstr(sname,".grd") != NULL)
  {
    if ((error = Plotfile(sname, -2, nn, &node, nb, &b_name, &nt, &t_n, &nq, &q_n, &ngon, &ngon_n,
                ntet, &tet_n, npyr, &pyr_n, npri, &pri_n, nhex, &hex_n, nply, &poly_n)) != 0)
    {
      fprintf(stderr,"\nError reading ASCII Fieldview file <%s>.",sname);
      fflush(stderr);
      exit(0);
    }
  } else if (strstr(sname,".inp") != NULL)
  {
    if ((error = StarCD_read(sname, nn, &node, nb, &b_name, &nt, &t_n, &nq, &q_n, &ngon, &ngon_n,
                nvol, &vol_name, ntet, &tet_n, &tet_vc, npyr, &pyr_n, &pyr_vc, npri, &pri_n, &pri_vc,
                nhex, &hex_n, &hex_vc, nply, &poly_n, &poly_vc)) != 0)
    {
      fprintf(stderr,"\nError reading STAR-CD file <%s>.",sname);
      fflush(stderr);
      exit(0);
    }
  } else if (strstr(sname,".nas") != NULL)
  {
    if ((error = Nastran(sname, -1, nn, &node, nb, &b_name, &nt, &t_n, &nq, &q_n, &ngon, &ngon_n,
                ntet, &tet_n, npyr, &pyr_n, npri, &pri_n, nhex, &hex_n, nply, &poly_n)) != 0)
    {
      fprintf(stderr,"\nError reading Nastran file <%s>.",sname);
      fflush(stderr);
      exit(0);
    }
  } else if (strstr(sname,".mesh3D") != NULL)
  {
    if ((error = Generic_mesh(sname, -1, nn, &node, nb, &b_name, &nt, &t_n, &nq, &q_n, &ngon, &ngon_n,
                ntet, &tet_n, npyr, &pyr_n, npri, &pri_n, nhex, &hex_n, nply, &poly_n)) != 0)
    {
      fprintf(stderr,"\nError reading ASCII Generic file <%s>.",sname);
      fflush(stderr);
      exit(0);
    }
  } else
  {
    fprintf(stdout,"\nInput file type not identified for file <%s>",sname);
    error = -1;
  }

  fprintf(stdout,"\n\nInput mesh statistics:");
  fprintf(stdout,"\n # of nodes = %d",nn);
  if (ntet > 0) fprintf(stdout,"\n # of tetrahedra = %d",ntet);
  if (npyr > 0) fprintf(stdout,"\n # of pyramids   = %d",npyr);
  if (npri > 0) fprintf(stdout,"\n # of prisms     = %d",npri);
  if (nhex > 0) fprintf(stdout,"\n # of hexahedra  = %d",nhex);
  if (nply > 0) fprintf(stdout,"\n # of polyhedra  = %d",nply);
  if (nvol > 0) fprintf(stdout,"\n # of volume conditions  = %d",nvol);
  for (n=0; n < nvol; n++)
    fprintf(stdout,"\nVolume %i name = %s",n+1,vol_name[n]);
  fprintf(stdout,"\n # of boundaries = %d",nb);
  int nttot = 0;
  int nqtot = 0;
  int ngtot = 0;
  for (b=0; b < nb; b++)
  {
    fprintf(stdout,"\nBoundary %d name = %s",b+1,b_name[b]);
    if (  nt[b] > 0) fprintf(stdout,"\n  # of triangles      = %d",nt[b]);
    if (  nq[b] > 0) fprintf(stdout,"\n  # of quadrilaterals = %d",nq[b]);
    if (ngon[b] > 0) fprintf(stdout,"\n  # of polygons       = %d",ngon[b]);
    nttot += nt[b];
    nqtot += nq[b];
    ngtot += ngon[b];
  }
  if (nttot > 0) fprintf(stdout,"\nTotal number of surface triangles      = %d",nttot);
  if (nqtot > 0) fprintf(stdout,"\nTotal number of surface quadrilaterals = %d",nqtot);
  if (ngtot > 0) fprintf(stdout,"\nTotal number of surface polygons       = %d",ngtot);
  fprintf(stdout,"\n\n");
  fprintf(stdout,"\n");
  fflush(stdout);

  x = (double*)malloc(nn*sizeof(double));
  y = (double*)malloc(nn*sizeof(double));
  z = (double*)malloc(nn*sizeof(double));
  for (n=0; n < nn; n++)
  {
    x[n] = node[n][0];
    y[n] = node[n][1];
    z[n] = node[n][2];
  }

  free(node);

  return error;
}

int Write_Mesh(char sname[], char **&b_name, int &nn, double *&x, double *&y, double *&z,
               int &nb, int *&nt, int *&nq, int *&ngon, int ***&t_n, int ***&q_n, int ***&ngon_n,
               int &ntet, int **&tet_n, int &npyr, int **&pyr_n,
               int &npri, int **&pri_n, int &nhex, int **&hex_n, int &nply, int ***&poly_n)
{
  int nvol;
  int *tet_vc, *pyr_vc, *pri_vc, *hex_vc, *poly_vc;
  char **vol_name;
  int error = 0;

  tet_vc = 0;
  pyr_vc = 0;
  pri_vc = 0;
  hex_vc = 0;
  poly_vc = 0;
  vol_name = (char**)malloc(sizeof(char*));
  vol_name[0] = (char*)malloc(7*sizeof(char));
  sprintf(vol_name[0],"Volume");
  nvol = 1;

  error = Write_Mesh(sname,b_name,nn,x,y,z,nb,nt,nq,ngon,t_n,q_n,ngon_n,ntet,tet_n,npyr,pyr_n,npri,pri_n,nhex,hex_n,nply,poly_n,
             nvol,vol_name,tet_vc,pyr_vc,pri_vc,hex_vc,poly_vc);

  free(vol_name[0]);
  free(vol_name);

  return error;
}

int Write_Mesh(char sname[], char **&b_name, int &nn, double *&x, double *&y, double *&z,
                int &nb, int *&nt, int *&nq, int *&ngon, int ***&t_n, int ***&q_n, int ***&ngon_n,
                int &ntet, int **&tet_n, int &npyr, int **&pyr_n,
                int &npri, int **&pri_n, int &nhex, int **&hex_n, int &nply, int ***&poly_n,
                int &nvol, char **&vol_name, int *&tet_vc, int *&pyr_vc, int *&pri_vc, int *&hex_vc, int *&poly_vc)
{
  int b, i, j, n, n0, n1, n2, n3;
  FILE *fp;
  char boundary_name[132];
  int error = 0;

  Point *node = (Point*)malloc(nn*sizeof(Point));

  for (n=0; n < nn; n++)
    node[n] = Point(x[n],y[n],z[n]);

  if (strstr(sname,".sg") != NULL)
  {
    // write serial mesh file
    error = SimGrid_write(sname, nn, node, nb, b_name, nvol, vol_name, nt, t_n, nq, q_n, ngon, ngon_n,
                ntet, tet_n, tet_vc, npyr, pyr_n, pyr_vc, npri, pri_n, pri_vc, nhex, hex_n, hex_vc, nply, poly_n, poly_vc);
  } else if (strstr(sname,".cgns") != NULL)
  {
    // write serial mesh file
#ifdef HAVE_CGNS
    error = CGNS_write(sname, nn, node, nb, b_name, nt, t_n, nq, q_n, ngon, ngon_n,
                ntet, tet_n, npyr, pyr_n, npri, pri_n, nhex, hex_n, nply, poly_n);
#else
    fprintf(stderr,"\nCGNS not enabled!");
    fflush(stderr);
    exit(0);
#endif
  } else if (strstr(sname,".uns") != NULL || strstr(sname,".fv") != NULL)
  {
    error = Plotfile(sname, 1, nn, &node, nb, &b_name, &nt, &t_n, &nq, &q_n, &ngon, &ngon_n,
              ntet, &tet_n, npyr, &pyr_n, npri, &pri_n, nhex, &hex_n, nply, &poly_n);
  } else if (strstr(sname,".crunch") != NULL || strstr(sname,".grd") != NULL)
  {
    error = Plotfile(sname, 2, nn, &node, nb, &b_name, &nt, &t_n, &nq, &q_n, &ngon, &ngon_n,
              ntet, &tet_n, npyr, &pyr_n, npri, &pri_n, nhex, &hex_n, nply, &poly_n);
  } else if (strstr(sname,".case") != NULL)
  {
    error = Ensight(sname, 1, nn, &node, nb, &b_name, &nt, &t_n, &nq, &q_n, &ngon, &ngon_n,
              ntet, &tet_n, npyr, &pyr_n, npri, &pri_n, nhex, &hex_n, nply, &poly_n);
  } else if (strstr(sname,".inp") != NULL)
  {
    error = StarCD_write(sname, nn, node, nb, b_name, nt, t_n, nq, q_n, ngon, ngon_n,
                         nvol, vol_name, ntet, tet_n, tet_vc, npyr, pyr_n, pyr_vc,
                         npri, pri_n, pri_vc, nhex, hex_n, hex_vc, nply, poly_n, poly_vc);
  } else if (strstr(sname,".nas") != NULL)
  {
    error = Nastran(sname, 1, nn, &node, nb, &b_name, &nt, &t_n, &nq, &q_n, &ngon, &ngon_n,
              ntet, &tet_n, npyr, &pyr_n, npri, &pri_n, nhex, &hex_n, nply, &poly_n);
  } else if (strstr(sname,".mesh3D") != NULL)
  {
    error = Generic_mesh(sname, 1, nn, &node, nb, &b_name, &nt, &t_n, &nq, &q_n, &ngon, &ngon_n,
              ntet, &tet_n, npyr, &pyr_n, npri, &pri_n, nhex, &hex_n, nply, &poly_n);
  } else if (strstr(sname,".tri") != NULL)
  {
    // generate geometry file from boundary data
    if ((fp=fopen(sname,"w")) == NULL)
    {
      printf("\nCould not open geometry file <%s>",sname);
      exit(0);
    }

    int nf, np, *map, *component;
    map = new int[nn];
    for (n=0; n < nn; n++)
      map[n] = 0;
    nf = 0;
    for (b=0; b < nb; b++)
    {
      for (i=0; i < nt[b]; i++)
      {
        for (j=0; j < 3; j++)
        {
          n = t_n[b][i][j];
          map[n] = 1;
        }
        nf++;
      }
      for (i=0; i < nq[b]; i++)
      {
        for (j=0; j < 4; j++)
        {
          n = q_n[b][i][j];
          map[n] = 1;
        }
        n0 = q_n[b][i][0];
        n1 = q_n[b][i][1];
        n2 = q_n[b][i][2];
        n3 = q_n[b][i][3];
        if (n0 == n1)
        {
          nf++;
        } else if (n1 == n2)
        {
          nf++;
        } else if (n2 == n3)
        {
          nf++;
        } else if (n3 == n0)
        {
          nf++;
        } else
        {
          nf+=2;
        }
      }
    }
    np = 0;
    for (n=0; n < nn; n++)
      if (map[n] == 1)
        map[n] = ++np;

    fprintf(fp,"%d %d\n", np, nf);
    for (n=0; n < nn; n++)
    {
      if (map[n] > 0)
        fprintf(fp,"%.8g %.8g %.8g %i\n",node[n][0],node[n][1],node[n][2],n+1);
    }

    component = new int[nf];
    nf = 0;
    for (b=0; b < nb; b++)
    {
      printf("\nWriting geometry for boundary #%d, name = %s",b+1,b_name[b]);
      for (i=0; i < nt[b]; i++)
      {
        n0 = t_n[b][i][0];
        n1 = t_n[b][i][1];
        n2 = t_n[b][i][2];
        fprintf(fp,"%d %d %d\n",map[n0],map[n1],map[n2]);
        component[nf++] = b+1;
      }
      for (i=0; i < nq[b]; i++)
      {
        n0 = q_n[b][i][0];
        n1 = q_n[b][i][1];
        n2 = q_n[b][i][2];
        n3 = q_n[b][i][3];
        if (n0 == n1)
        {
          fprintf(fp,"%d %d %d\n",map[n1],map[n2],map[n3]);
          component[nf++] = b+1;
        } else if (n1 == n2)
        {
          fprintf(fp,"%d %d %d\n",map[n2],map[n3],map[n0]);
          component[nf++] = b+1;
        } else if (n2 == n3)
        {
          fprintf(fp,"%d %d %d\n",map[n3],map[n0],map[n1]);
          component[nf++] = b+1;
        } else if (n3 == n0)
        {
          fprintf(fp,"%d %d %d\n",map[n0],map[n1],map[n2]);
          component[nf++] = b+1;
        } else
        {
          fprintf(fp,"%d %d %d\n",map[n0],map[n1],map[n2]);
          component[nf++] = b+1;
          fprintf(fp,"%d %d %d\n",map[n0],map[n2],map[n3]);
          component[nf++] = b+1;
        }
      }
    }
    for (n=0; n < nf; n++)
      fprintf(fp,"%d\n",component[n]);

    fclose(fp);

    printf("\n");

    if ((fp=fopen("geometry.params","w")) == NULL)
    {
      printf("\nCould not open parameter file");
      exit(0);
    }
    fprintf(fp,"Number of boundaries\n");
    fprintf(fp,"%d\n",nb);
    fprintf(fp,"Boundary #   Layers      g_space       n_space              boundary name\n");
    fprintf(fp,"---------- ---------- -------------- -------------- --------------------------------\n");
    for (i=0; i < nb; i++)
    {
      sprintf(boundary_name,"%s",b_name[i]);
      fprintf(fp,"%10d %10d %14.7g %14.7g %s\n",
                   i+1,0,1.0,0.000001,boundary_name);
    }
    fclose(fp);

    delete map;
    delete component;
  } else if (strstr(sname,".facet") != NULL)
  {
    // generate geometry file from boundary data
    if ((fp=fopen(sname,"w")) == NULL)
    {
      printf("\nCould not open geometry file <%s>",sname);
      exit(0);
    }

    int nf, np, *map;
    map = new int[nn];
    for (n=0; n < nn; n++)
      map[n] = 0;
    nf = 0;
    for (b=0; b < nb; b++)
    {
      for (i=0; i < nt[b]; i++)
      {
        for (j=0; j < 3; j++)
        {
          n = t_n[b][i][j];
          map[n] = 1;
        }
        nf++;
      }
      for (i=0; i < nq[b]; i++)
      {
        for (j=0; j < 4; j++)
        {
          n = q_n[b][i][j];
          map[n] = 1;
        }
        n0 = q_n[b][i][0];
        n1 = q_n[b][i][1];
        n2 = q_n[b][i][2];
        n3 = q_n[b][i][3];
        if (n0 == n1)
        {
          nf++;
        } else if (n1 == n2)
        {
          nf++;
        } else if (n2 == n3)
        {
          nf++;
        } else if (n3 == n0)
        {
          nf++;
        } else
        {
          nf+=2;
        }
      }
    }
    np = 0;
    for (n=0; n < nn; n++)
      if (map[n] == 1)
        map[n] = ++np;

    fprintf(fp,"FACET FILE V3.0   Boundary data from Voxel_Mesh in Xpatch format\n");
    fprintf(fp,"1\n");
    fprintf(fp,"Grid\n");
    fprintf(fp,"0, 0.00 0.00 0.00 0.00\n");
    fprintf(fp,"%i\n",np);

    for (n=0; n < nn; n++)
    {
      if (map[n] > 0)
        fprintf(fp,"%22.15e %22.15e %22.15e\n",node[n][0],node[n][1],node[n][2]);
    }

    fprintf(fp,"1\n");
    fprintf(fp,"Triangles\n");
    fprintf(fp,"%i 3\n",nf);

    nf = 0;
    for (b=0; b < nb; b++)
    {
      printf("\nWriting geometry for boundary #%d, name = %s",b+1,b_name[b]);
      for (i=0; i < nt[b]; i++)
      {
        n0 = t_n[b][i][0];
        n1 = t_n[b][i][1];
        n2 = t_n[b][i][2];
        fprintf(fp,"%7i %7i %7i 0  %04i  %i\n",map[n0],map[n1],map[n2],b+1,++nf);
      }
      for (i=0; i < nq[b]; i++)
      {
        n0 = q_n[b][i][0];
        n1 = q_n[b][i][1];
        n2 = q_n[b][i][2];
        n3 = q_n[b][i][3];
        if (n0 == n1)
        {
          fprintf(fp,"%7i %7i %7i 0  %04i  %i\n",map[n1],map[n2],map[n3],b+1,++nf);
        } else if (n1 == n2)
        {
          fprintf(fp,"%7i %7i %7i 0  %04i  %i\n",map[n2],map[n3],map[n0],b+1,++nf);
        } else if (n2 == n3)
        {
          fprintf(fp,"%7i %7i %7i 0  %04i  %i\n",map[n3],map[n0],map[n1],b+1,++nf);
        } else if (n3 == n0)
        {
          fprintf(fp,"%7i %7i %7i 0  %04i  %i\n",map[n0],map[n1],map[n2],b+1,++nf);
        } else
        {
          fprintf(fp,"%7i %7i %7i 0  %04i  %i\n",map[n0],map[n1],map[n2],b+1,++nf);
          fprintf(fp,"%7i %7i %7i 0  %04i  %i\n",map[n0],map[n2],map[n3],b+1,++nf);
        }
      }
    }
    fclose(fp);

    printf("\n");

    if ((fp=fopen("geometry.params","w")) == NULL)
    {
      printf("\nCould not open parameter file");
      exit(0);
    }
    fprintf(fp,"Number of boundaries\n");
    fprintf(fp,"%d\n",nb);
    fprintf(fp,"Boundary #   Layers      g_space       n_space              boundary name\n");
    fprintf(fp,"---------- ---------- -------------- -------------- --------------------------------\n");
    for (i=0; i < nb; i++)
    {
      sprintf(boundary_name,"%s",b_name[i]);
      fprintf(fp,"%10d %10d %14.7g %14.7g %s\n",
                   i+1,0,1.0,0.000001,boundary_name);
    }
    fclose(fp);

    delete map;
  } else
  {
    fprintf(stdout,"\nOutput file type unidentified for file <%s>",sname);
    error = -1;
  }

  free(node);

  fprintf(stdout,"\n");

  return error;
}
