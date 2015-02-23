#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <sys/param.h>
#include "List.h"
#include "Point.h"
#include "Vector.h"
#include "journal.h"
#include "metis.h"
#ifdef HAVE_CGNS
#include "CGNS.h"
#endif
#include "SGIO.h"
#include "convert.h"
#include "merge.h"
#include "split_tree_partition.h"

int decomp(char *fname, int nparts, int model, int pmode)
{
  int b, c, i, j, k, n, n0, n1, n2, n3, n4, n5, n6, n7, p;
  int nn, nb, ntet, npyr, npri, nhex, nply, nvol;
  int *nt, *nq, *ngon, ***tri_conn, ***quad_conn, ***ngon_conn;
  int **tet_conn, **pyr_conn, **pri_conn, **hex_conn, ***poly_conn;
  int **node_map, **tri_map, **quad_map, **ngon_map, *tet_map, *pyr_map, *pri_map, *hex_map, *poly_map;
  int *tet_vc, *pyr_vc, *pri_vc, *hex_vc, *poly_vc;
  int *ntmp, *index;
  Point *node;
  char **b_name, **vol_name;
  const int bdim = 132;
  char buff[bdim], pname[bdim];
  int error = 0;
  FILE *part_f;

  // set all counts to 0
  nvol = nn = nb = ntet = npyr = npri = nhex = nply = 0;
  // set all pointers to null
  b_name = 0;
  vol_name = 0;
  nt = 0;
  nq = 0;
  ngon = 0;
  tri_conn = 0;
  quad_conn = 0;
  ngon_conn = 0;
  tet_conn = 0;
  pyr_conn = 0;
  pri_conn = 0;
  hex_conn = 0;
  poly_conn = 0;
  node_map = 0;
  tri_map = 0;
  quad_map = 0;
  ngon_map = 0;
  tet_map = 0;
  pyr_map = 0;
  pri_map = 0;
  hex_map = 0;
  poly_map = 0;
  node = 0;
  tet_vc = 0;
  pyr_vc = 0;
  pri_vc = 0;
  hex_vc = 0;
  poly_vc = 0;

  // read serial mesh file
  if (strstr(fname,".sg") != NULL)
  {
    if ((error = SimGrid_read(fname,nn,&node,nb,&b_name,nvol,&vol_name,&nt,&tri_conn,
         &nq,&quad_conn,&ngon,&ngon_conn,ntet,&tet_conn,&tet_vc,npyr,&pyr_conn,&pyr_vc,
         npri,&pri_conn,&pri_vc,nhex,&hex_conn,&hex_vc,nply,&poly_conn,&poly_vc)) != 0)
    {
      fprintf(stderr,"\nError reading serial SG file <%s>.",fname);
      fflush(stderr);
      return(error);
      //exit(0);
    }
  } else if (strstr(fname,".cgns") != NULL)
  {
#ifdef HAVE_CGNS
    if ((error = CGNS_read(fname, nn, &node, nb, &b_name, &nt, &tri_conn, &nq, &quad_conn, &ngon, &ngon_conn,
                ntet, &tet_conn, npyr, &pyr_conn, npri, &pri_conn, nhex, &hex_conn, nply, &poly_conn)) != 0)
    {
      fprintf(stderr,"\nError reading serial CGNS file <%s>.",fname);
      fflush(stderr);
      return(error);
      //exit(0);
    }
#else
      fprintf(stderr,"\nCGNS not enabled!");
      fflush(stderr);
      exit(0);
#endif
    if (ntet > 0)
      tet_vc = (int*)calloc(ntet,sizeof(int));
    if (npyr > 0)
      pyr_vc = (int*)calloc(npyr,sizeof(int));
    if (npri > 0)
      pri_vc = (int*)calloc(npri,sizeof(int));
    if (nhex > 0)
      hex_vc = (int*)calloc(nhex,sizeof(int));
    if (nply > 0)
      poly_vc = (int*)calloc(nply,sizeof(int));
  } else
  {
    fprintf(stderr,"\nFile with correct suffix not found. File = %s",fname);
    fflush(stderr);
    return(-1);
  }
  
  int *part = 0; // partition array
  part = new int[nn];

  if (pmode == 1)
  {
    if (model == 0)
    {
      // partition mesh using nodes with Metis

      // create node-to-node hash table for partitioning
      List **hash;
      hash = new List*[nn];
      for (n=0; n < nn; n++)
      {
        hash[n] = new List();
        hash[n]->Add_To_List(n);
      }
      for (c=0; c < ntet; c++)
      {
        n0 = tet_conn[c][0];
        n1 = tet_conn[c][1];
        n2 = tet_conn[c][2];
        n3 = tet_conn[c][3];
        hash[n0]->Check_List(n1);
        hash[n0]->Check_List(n2);
        hash[n0]->Check_List(n3);
        hash[n1]->Check_List(n0);
        hash[n1]->Check_List(n2);
        hash[n1]->Check_List(n3);
        hash[n2]->Check_List(n0);
        hash[n2]->Check_List(n1);
        hash[n2]->Check_List(n3);
        hash[n3]->Check_List(n0);
        hash[n3]->Check_List(n1);
        hash[n3]->Check_List(n2);
      }
      for (c=0; c < npyr; c++)
      {
        n0 = pyr_conn[c][0];
        n1 = pyr_conn[c][1];
        n2 = pyr_conn[c][2];
        n3 = pyr_conn[c][3];
        n4 = pyr_conn[c][4];
        hash[n0]->Check_List(n1);
        hash[n0]->Check_List(n3);
        hash[n0]->Check_List(n4);
        hash[n1]->Check_List(n0);
        hash[n1]->Check_List(n2);
        hash[n1]->Check_List(n4);
        hash[n2]->Check_List(n1);
        hash[n2]->Check_List(n3);
        hash[n2]->Check_List(n4);
        hash[n3]->Check_List(n0);
        hash[n3]->Check_List(n2);
        hash[n3]->Check_List(n4);
        hash[n4]->Check_List(n0);
        hash[n4]->Check_List(n1);
        hash[n4]->Check_List(n2);
        hash[n4]->Check_List(n3);
      }
      for (c=0; c < npri; c++)
      {
        n0 = pri_conn[c][0];
        n1 = pri_conn[c][1];
        n2 = pri_conn[c][2];
        n3 = pri_conn[c][3];
        n4 = pri_conn[c][4];
        n5 = pri_conn[c][5];
        hash[n0]->Check_List(n1);
        hash[n0]->Check_List(n2);
        hash[n0]->Check_List(n3);
        hash[n1]->Check_List(n0);
        hash[n1]->Check_List(n2);
        hash[n1]->Check_List(n4);
        hash[n2]->Check_List(n0);
        hash[n2]->Check_List(n1);
        hash[n2]->Check_List(n5);
        hash[n3]->Check_List(n0);
        hash[n3]->Check_List(n4);
        hash[n3]->Check_List(n5);
        hash[n4]->Check_List(n1);
        hash[n4]->Check_List(n3);
        hash[n4]->Check_List(n5);
        hash[n5]->Check_List(n2);
        hash[n5]->Check_List(n3);
        hash[n5]->Check_List(n4);
      }
      for (c=0; c < nhex; c++)
      {
        n0 = hex_conn[c][0];
        n1 = hex_conn[c][1];
        n2 = hex_conn[c][2];
        n3 = hex_conn[c][3];
        n4 = hex_conn[c][4];
        n5 = hex_conn[c][5];
        n6 = hex_conn[c][6];
        n7 = hex_conn[c][7];
        hash[n0]->Check_List(n1);
        hash[n0]->Check_List(n3);
        hash[n0]->Check_List(n4);
        hash[n1]->Check_List(n0);
        hash[n1]->Check_List(n2);
        hash[n1]->Check_List(n5);
        hash[n2]->Check_List(n1);
        hash[n2]->Check_List(n3);
        hash[n2]->Check_List(n6);
        hash[n3]->Check_List(n0);
        hash[n3]->Check_List(n2);
        hash[n3]->Check_List(n7);
        hash[n4]->Check_List(n0);
        hash[n4]->Check_List(n5);
        hash[n4]->Check_List(n7);
        hash[n5]->Check_List(n1);
        hash[n5]->Check_List(n4);
        hash[n5]->Check_List(n6);
        hash[n6]->Check_List(n2);
        hash[n6]->Check_List(n5);
        hash[n6]->Check_List(n7);
        hash[n7]->Check_List(n3);
        hash[n7]->Check_List(n4);
        hash[n7]->Check_List(n6);
      }
      for (c=0; c < nply; c++)
      {
        for (i=1; i <= poly_conn[c][0][0]; i++)
        {
          for (j=1; j <= poly_conn[c][i][0]; j++)
          {
            n0 = poly_conn[c][i][j];
            if (j < poly_conn[c][i][0])
              n1 = poly_conn[c][i][j+1];
            else
              n1 = poly_conn[c][i][1];
            hash[n0]->Check_List(n1);
            hash[n1]->Check_List(n0);
          }
        }
      }

      // transfer hash table to compressed row format
      int *ia, *ja, *adjwgt = NULL;
      ia = new int[nn+1];
      int mdim;
      for (mdim=n=0; n < nn; n++)
      {
        k = hash[n]->max;
        ia[n] = mdim;
        mdim += k;
      }
      ia[nn] = mdim;
      
      ja = new int[mdim];
      for (mdim=n=0; n < nn; n++)
        for (j=0; j < hash[n]->max; j++)
          ja[mdim++] = hash[n]->list[j];

      // free hash table memory
      for (n=0; n < nn; n++)
        delete hash[n];
      delete[] hash;

      int wgtflag = 2; // weights on vertices only
      int numflag = 0; // c++ style
      int *options = 0; // set up other options
      options = new int[5];
      for (i=0; i < 5; i++)
        options[i] = 0;
      int edgecut = 0; // number of edges cut by partition
      int *vwgt = new int[nn]; // weight array
      for (n=0; n < nn; n++)
      {
        part[n] = 0;
        vwgt[n] = 1;
      }
      // weight boundary nodes more
      //for (b=0; b < nb; b++)
      //{
      //  for (c=0; c < nt[b]; c++)
      //    for (i=0; i < 3; i++)
      //      vwgt[tri_conn[b][c][i]] = 2;
      //  for (c=0; c < nq[b]; c++)
      //    for (i=0; i < 4; i++)
      //      vwgt[quad_conn[b][c][i]] = 2;
      //}

      //if (nparts <= 8)
      //  METIS_PartGraphRecursive(&nn,ia,ja,vwgt,adjwgt,&wgtflag,&numflag,&nparts,options,&edgecut,part);
      //else
        METIS_PartGraphKway(&nn,ia,ja,vwgt,adjwgt,&wgtflag,&numflag,&nparts,options,&edgecut,part);

      // free memory
      delete[] vwgt;
      delete[] ia;
      delete[] ja;
      delete[] options;

    } else if (model == 1)
    {
      // use split tree partitioning
      printf("\nUsing balanced split-tree partitioning for %d nodes.",nn);
      fflush(stdout);
      n = 0;
      if (split_tree_partition(nn,node,nparts,part,n,stdout) == 0)
      {
        printf("\nPartitioning failed.!");
        fflush(stdout);
        exit(0);
      }
    } else
    {
      // use split tree partitioning
      printf("\nUsing spacial split-tree partitioning for %d nodes.",nn);
      fflush(stdout);
      n = 1;
      if (split_tree_partition(nn,node,nparts,part,n,stdout) == 0)
      {
        printf("\nPartitioning failed.!");
        fflush(stdout);
        exit(0);
      }
    }

    printf("\nPartitioning complete.\n");
    fflush(stdout);

    // create the parallel filename
    sprintf(buff,"%s",fname);
    if (strstr(fname,".sg") != NULL)
    {
      char *ptr = strstr(buff,".sg");
      if (ptr == NULL)
      {
        printf("\nSG suffix <.sg> not found in file name!");
        fflush(stdout);
        exit(0);
      } else
        *ptr = '\0';
    } else if (strstr(fname,".cgns") != NULL)
    {
      char *ptr = strstr(buff,".cgns");
      if (ptr == NULL)
      {
        printf("\nCGNS suffix <.cgns> not found in file name!");
        fflush(stdout);
        exit(0);
      } else
        *ptr = '\0';
    }
    sprintf(pname,"%s.partitions",buff);

    if ((part_f=fopen(pname,"w")) == NULL)
    {
      printf("\nCouldn't open file partition file < %s >",pname);
      exit(0);
    }
    for (n=0; n < nn; n++)
      fprintf(part_f,"%i\n",part[n]);

    fclose(part_f);

  } else if (pmode == 2)
  {
  
    // create the parallel filename
    sprintf(buff,"%s",fname);
    if (strstr(fname,".sg") != NULL)
    {
      char *ptr = strstr(buff,".sg");
      if (ptr == NULL)
      {
        printf("\nSG suffix <.sg> not found in file name!");
        fflush(stdout);
        exit(0);
      } else
        *ptr = '\0';
    } else if (strstr(fname,".cgns") != NULL)
    {
      char *ptr = strstr(buff,".cgns");
      if (ptr == NULL)
      {
        printf("\nCGNS suffix <.cgns> not found in file name!");
        fflush(stdout);
        exit(0);
      } else
        *ptr = '\0';
    }
    sprintf(pname,"%s.partitions",buff);

    if ((part_f=fopen(pname,"r")) == NULL)
    {
      printf("\nCouldn't open file partition file < %s >",pname);
      exit(0);
    }
    nparts = 0;
    for (n=0; n < nn; n++)
    {
      fscanf(part_f,"%i",&part[n]);
      nparts = MAX(nparts,part[n]+1);
    }

    fclose(part_f);

    printf("\nPartitions read from file, number of partitions = %i\n",nparts);
    fflush(stdout);

    // create list of nodes owned by each processor
    List **plist;
    plist = new List*[nparts];
    for (p=0; p < nparts; p++)
      plist[p] = new List();
    for (n=0; n < nn; n++)
      plist[part[n]]->Add_To_List(n);

    // allocate space for partitioned node and element arrays
    // maximum per processor is the total number in the mesh
    Point *pnode;
    int *pnt, *pnq, *png, ***ptri, ***pquad, ***pngon, **ptet, **ppyr, **ppri, **phex, ***ppoly;
    int *ptet_vc, *ppyr_vc, *ppri_vc, *phex_vc, *ppoly_vc;
    pnode = (Point*)malloc(nn*sizeof(Point));
    node_map = (int**)malloc(nn*sizeof(int*));
    for (i=0; i < nn; i++)
      node_map[i] = (int*)malloc(3*sizeof(int));
    ntmp = new int[nn];
    index = new int[nn];

    pnt = (int*)malloc(nb*sizeof(int));
    pnq = (int*)malloc(nb*sizeof(int));
    png = (int*)malloc(nb*sizeof(int));
    ptri = (int***)malloc(nb*sizeof(int**));
    pquad = (int***)malloc(nb*sizeof(int**));
    pngon = (int***)malloc(nb*sizeof(int**));
    tri_map = (int**)malloc(nb*sizeof(int*));
    quad_map = (int**)malloc(nb*sizeof(int*));
    ngon_map = (int**)malloc(nb*sizeof(int*));
    for (b=0; b < nb; b++)
    {
      if (nt[b] > 0)
      {
        ptri[b] = (int**)malloc(nt[b]*sizeof(int*));
        for (i=0; i < nt[b]; i++)
          ptri[b][i] = (int*)malloc(3*sizeof(int));
        tri_map[b] = (int*)malloc(nt[b]*sizeof(int));
      } else
      {
        ptri[b] = 0;
        tri_map[b] = 0;
      }
      if (nq[b] > 0)
      {
        pquad[b] = (int**)malloc(nq[b]*sizeof(int*));
        for (i=0; i < nq[b]; i++)
          pquad[b][i] = (int*)malloc(4*sizeof(int));
        quad_map[b] = (int*)malloc(nq[b]*sizeof(int));
      } else
      {
        pquad[b] = 0;
        quad_map[b] = 0;
      }
      if (ngon[b] > 0)
      {
        pngon[b] = (int**)malloc(ngon[b]*sizeof(int*));
        for (i=0; i < ngon[b]; i++)
          pngon[b][i] = 0;  // will be allocated later for each partition list
        ngon_map[b] = (int*)malloc(ngon[b]*sizeof(int));
      } else
      {
        pngon[b] = 0;
        ngon_map[b] = 0;
      }
    }
      
    if (ntet > 0)
    {
      ptet_vc = (int*)malloc(ntet*sizeof(int));
      ptet = (int**)malloc(ntet*sizeof(int*));
      for (i=0; i < ntet; i++)
        ptet[i] = (int*)malloc(4*sizeof(int));
      tet_map = (int*)malloc(ntet*sizeof(int));
    } else
    {
      ptet_vc = 0;
      ptet = 0;
      tet_map = 0;
    }

    if (npyr > 0)
    {
      ppyr_vc = (int*)malloc(npyr*sizeof(int));
      ppyr = (int**)malloc(npyr*sizeof(int*));
      for (i=0; i < npyr; i++)
        ppyr[i] = (int*)malloc(5*sizeof(int));
      pyr_map = (int*)malloc(npyr*sizeof(int));
    } else
    {
      ppyr_vc = 0;
      ppyr = 0;
      pyr_map = 0;
    }

    if (npri > 0)
    {
      ppri_vc = (int*)malloc(npri*sizeof(int));
      ppri = (int**)malloc(npri*sizeof(int*));
      for (i=0; i < npri; i++)
        ppri[i] = (int*)malloc(6*sizeof(int));
      pri_map = (int*)malloc(npri*sizeof(int));
    } else
    {
      ppri_vc = 0;
      ppri = 0;
      pri_map = 0;
    }

    if (nhex > 0)
    {
      phex_vc = (int*)malloc(nhex*sizeof(int));
      phex = (int**)malloc(nhex*sizeof(int*));
      for (i=0; i < nhex; i++)
        phex[i] = (int*)malloc(8*sizeof(int));
      hex_map = (int*)malloc(nhex*sizeof(int));
    } else
    {
      phex_vc = 0;
      phex = 0;
      hex_map = 0;
    }

    if (nply > 0)
    {
      ppoly_vc = (int*)malloc(nply*sizeof(int));
      ppoly = (int***)malloc(nply*sizeof(int**));
      for (i=0; i < nply; i++)
        ppoly[i] = 0; // will be allocated later for each partition list
      poly_map = (int*)malloc(nply*sizeof(int));
    } else
    {
      ppoly_vc = 0;
      ppoly = 0;
      poly_map = 0;
    }

    int pnn, pntet, pnpyr, pnpri, pnhex, pnply;

    for (p = 0; p < nparts; p++)
    {

      for (b=0; b < nb; b++)
      {
        pnt[b] = pnq[b] = png[b] = 0;
        for (c=0; c < nt[b]; c++)
        {
          for (j=i=0; i < 3 && !j; i++)
            if (part[tri_conn[b][c][i]] == p)
              j=1;
          if (j)
          {
            ptri[b][pnt[b]][0] = tri_conn[b][c][0];
            ptri[b][pnt[b]][1] = tri_conn[b][c][1];
            ptri[b][pnt[b]][2] = tri_conn[b][c][2];
            tri_map[b][pnt[b]] = c;
            pnt[b]++;
          }
        }
        for (c=0; c < nq[b]; c++)
        {
          for (j=i=0; i < 4 && !j; i++)
            if (part[quad_conn[b][c][i]] == p)
              j=1;
          if (j)
          {
            pquad[b][pnq[b]][0] = quad_conn[b][c][0];
            pquad[b][pnq[b]][1] = quad_conn[b][c][1];
            pquad[b][pnq[b]][2] = quad_conn[b][c][2];
            pquad[b][pnq[b]][3] = quad_conn[b][c][3];
            quad_map[b][pnq[b]] = c;
            pnq[b]++;
          }
        }
        for (c=0; c < ngon[b]; c++)
        {
          for (j=0,i=1; i <= ngon_conn[b][c][0] && !j; i++)
            if (part[ngon_conn[b][c][i]] == p)
              j=1;
          if (j)
          {
            pngon[b][png[b]] = (int*)realloc((void*)pngon[b][png[b]],(ngon_conn[b][c][0]+1)*sizeof(int));
            pngon[b][png[b]][0] = ngon_conn[b][c][0];
            for (k=1; k <= ngon_conn[b][c][0]; k++)
              pngon[b][png[b]][k] = ngon_conn[b][c][k];
            ngon_map[b][png[b]] = c;
            png[b]++;
          }
        }
      }

      pnn = pntet = pnpyr = pnpri = pnhex = pnply = 0;

      for (c=0; c < ntet; c++)
      {
        for (j=i=0; i < 4 && !j; i++)
          if (part[tet_conn[c][i]] == p)
            j=1;
        if (j)
        {
          ptet[pntet][0] = tet_conn[c][0];
          ptet[pntet][1] = tet_conn[c][1];
          ptet[pntet][2] = tet_conn[c][2];
          ptet[pntet][3] = tet_conn[c][3];
          ptet_vc[pntet] = tet_vc[c];
          tet_map[pntet] = c;
          pntet++;
        }
      }
          
      for (c=0; c < npyr; c++)
      {
        for (j=i=0; i < 5 && !j; i++)
          if (part[pyr_conn[c][i]] == p)
            j=1;
        if (j)
        {
          ppyr[pnpyr][0] = pyr_conn[c][0];
          ppyr[pnpyr][1] = pyr_conn[c][1];
          ppyr[pnpyr][2] = pyr_conn[c][2];
          ppyr[pnpyr][3] = pyr_conn[c][3];
          ppyr[pnpyr][4] = pyr_conn[c][4];
          ppyr_vc[pnpyr] = pyr_vc[c];
          pyr_map[pnpyr] = c;
          pnpyr++;
        }
      }
          
      for (c=0; c < npri; c++)
      {
        for (j=i=0; i < 6 && !j; i++)
          if (part[pri_conn[c][i]] == p)
            j=1;
        if (j)
        {
          ppri[pnpri][0] = pri_conn[c][0];
          ppri[pnpri][1] = pri_conn[c][1];
          ppri[pnpri][2] = pri_conn[c][2];
          ppri[pnpri][3] = pri_conn[c][3];
          ppri[pnpri][4] = pri_conn[c][4];
          ppri[pnpri][5] = pri_conn[c][5];
          ppri_vc[pnpri] = pri_vc[c];
          pri_map[pnpri] = c;
          pnpri++;
        }
      }
          
      for (c=0; c < nhex; c++)
      {
        for (j=i=0; i < 8 && !j; i++)
          if (part[hex_conn[c][i]] == p)
            j=1;
        if (j)
        {
          phex[pnhex][0] = hex_conn[c][0];
          phex[pnhex][1] = hex_conn[c][1];
          phex[pnhex][2] = hex_conn[c][2];
          phex[pnhex][3] = hex_conn[c][3];
          phex[pnhex][4] = hex_conn[c][4];
          phex[pnhex][5] = hex_conn[c][5];
          phex[pnhex][6] = hex_conn[c][6];
          phex[pnhex][7] = hex_conn[c][7];
          phex_vc[pnhex] = hex_vc[c];
          hex_map[pnhex] = c;
          pnhex++;
        }
      }
          
      for (c=0; c < nply; c++)
      {
        for (k=0, i=1; i <= poly_conn[c][0][0] && !k; i++)
          for (j=1; j <= poly_conn[c][i][0] && !k; j++)
            if (part[poly_conn[c][i][j]] == p)
              k=1;
        if (k)
        {
          ppoly[pnply] = (int**)realloc((void*)ppoly[pnply],(poly_conn[c][0][0]+1)*sizeof(int*));
          ppoly[pnply][0] = (int*)malloc(sizeof(int));
          ppoly[pnply][0][0] = poly_conn[c][0][0];
          for (i=1; i <= poly_conn[c][0][0]; i++)
          {
            ppoly[pnply][i] = (int*)malloc((poly_conn[c][i][0]+1)*sizeof(int));
            ppoly[pnply][i][0] = poly_conn[c][i][0];
            for (j=1; j <= poly_conn[c][i][0]; j++)
              ppoly[pnply][i][j] = poly_conn[c][i][j];
          }
          ppoly_vc[pnply] = poly_vc[c];
          poly_map[pnply] = c;
          pnply++;
        }
      }

      // create node list for current processor
      List nlist, overlap;

      overlap.Redimension(0);
      for (n=0; n < nn; n++)
        index[n] = -1;

      // create unique list of overlap nodes
      for (n=0; n < nn; n++)
        ntmp[n] = 0;
      for (b=0; b < nb; b++)
      {
        for (c=0; c < pnt[b]; c++)
          for (i=0; i < 3; i++)
            if (part[n = ptri[b][c][i]] != p)
              ntmp[n] = 1;
        for (c=0; c < pnq[b]; c++)
          for (i=0; i < 4; i++)
            if (part[n = pquad[b][c][i]] != p)
              ntmp[n] = 1;
        for (c=0; c < png[b]; c++)
          for (i=1; i <= pngon[b][c][0]; i++)
            if (part[n = pngon[b][c][i]] != p)
              ntmp[n] = 1;
      }
      for (c=0; c < pntet; c++)
        for (i=0; i < 4; i++)
          if (part[n = ptet[c][i]] != p)
            ntmp[n] = 1;
      for (c=0; c < pnpyr; c++)
        for (i=0; i < 5; i++)
          if (part[n = ppyr[c][i]] != p)
            ntmp[n] = 1;
      for (c=0; c < pnpri; c++)
        for (i=0; i < 6; i++)
          if (part[n = ppri[c][i]] != p)
            ntmp[n] = 1;
      for (c=0; c < pnhex; c++)
        for (i=0; i < 8; i++)
          if (part[n = phex[c][i]] != p)
            ntmp[n] = 1;
      for (c=0; c < pnply; c++)
        for (i=1; i <= ppoly[c][0][0]; i++)
          for (j=1; j <= ppoly[c][i][0]; j++)
            if (part[n = ppoly[c][i][j]] != p)
              ntmp[n] = 1;
      for (n=0; n < nn; n++)
        if (ntmp[n] == 1)
          overlap.Add_To_List(n);

      // create list of partitioned nodes for this process
      nlist.Redimension(0);
      for (i=0; i < plist[p]->max; i++)
      {
        n=plist[p]->list[i];
        index[n] = nlist.max;
        nlist.Add_To_List(n);
      }

      // add overlap nodes to list
      for (i=0; i < overlap.max; i++)
      {
        n=overlap.list[i];
        index[n] = nlist.max;
        nlist.Add_To_List(n);
      }

      // modify partitioned connectivities with node map
      for (b=0; b < nb; b++)
      {
        for (c=0; c < pnt[b]; c++)
          for (i=0; i < 3; i++)
            ptri[b][c][i] = index[ptri[b][c][i]];
        for (c=0; c < pnq[b]; c++)
          for (i=0; i < 4; i++)
            pquad[b][c][i] = index[pquad[b][c][i]];
        for (c=0; c < png[b]; c++)
          for (i=1; i <= pngon[b][c][0]; i++)
            pngon[b][c][i] = index[pngon[b][c][i]];
      }
      for (c=0; c < pntet; c++)
        for (i=0; i < 4; i++)
          ptet[c][i] = index[ptet[c][i]];
      for (c=0; c < pnpyr; c++)
        for (i=0; i < 5; i++)
          ppyr[c][i] = index[ppyr[c][i]];
      for (c=0; c < pnpri; c++)
        for (i=0; i < 6; i++)
          ppri[c][i] = index[ppri[c][i]];
      for (c=0; c < pnhex; c++)
        for (i=0; i < 8; i++)
          phex[c][i] = index[phex[c][i]];
      for (c=0; c < pnply; c++)
        for (i=1; i <= ppoly[c][0][0]; i++)
          for (j=1; j <= ppoly[c][i][0]; j++)
            ppoly[c][i][j] = index[ppoly[c][i][j]];
      
      // now modify index to store the index of the processor the node is on
      for (i=0; i < nparts; i++)
        for (j=0; j < plist[i]->max; j++)
          index[plist[i]->list[j]] = j;

      // store nodes and map for current processor
      pnn = 0;
      for (i=0; i < nlist.max; i++)
      {
        n = nlist.list[i];
        pnode[i] = node[n];
        node_map[i][0] = n;
        node_map[i][1] = part[n];
        node_map[i][2] = index[n];
        //node_map[i][2] = plist[part[n]]->Index(n);
      }
      pnn = nlist.max;

      // create the parallel filename
      sprintf(buff,"%s",fname);
      if (strstr(fname,".sg") != NULL)
      {
        char *ptr = strstr(buff,".sg");
        if (ptr == NULL)
        {
          printf("\nSG suffix <.sg> not found in file name!");
          fflush(stdout);
          exit(0);
        } else
          *ptr = '\0';
        sprintf(pname,"%s_%d.sg",buff,p);
      } else if (strstr(fname,".cgns") != NULL)
      {
        char *ptr = strstr(buff,".cgns");
        if (ptr == NULL)
        {
          printf("\nCGNS suffix <.cgns> not found in file name!");
          fflush(stdout);
          exit(0);
        } else
          *ptr = '\0';
        sprintf(pname,"%s_%d.cgns",buff,p);
      }

      printf("\nWriting parallel file <%s>",pname);
      printf("\n Number of nodes = %d",pnn);
      if (pntet > 0) printf("\n Number of tetrahdedra = %d",pntet);
      if (pnpyr > 0) printf("\n Number of pyramids = %d",pnpyr);
      if (pnpri > 0) printf("\n Number of prisms = %d",pnpri);
      if (pnhex > 0) printf("\n Number of hexahedra = %d",pnhex);
      if (pnply > 0) printf("\n Number of polyhedra = %d",pnply);
      for (b=0; b < nb; b++)
      {
        if (pnt[b] > 0) printf("\n Boundary %d, number of triangles      = %d",b+1,pnt[b]);
        if (pnq[b] > 0) printf("\n Boundary %d, number of quadrilaterals = %d",b+1,pnq[b]);
        if (png[b] > 0) printf("\n Boundary %d, number of polygons       = %d",b+1,png[b]);
      }
      fflush(stdout);

      if (strstr(fname,".sg") != NULL)
      {
        P_SimGrid_write(1, pname, pnn, pnode, nb, b_name, nvol, vol_name, pnt, ptri, pnq, pquad, png, pngon,
                  pntet, ptet, ptet_vc, pnpyr, ppyr, ppyr_vc, pnpri, ppri, ppri_vc, pnhex, phex, phex_vc, pnply, ppoly, ppoly_vc,
                  node_map, tri_map, quad_map, ngon_map, tet_map, pyr_map, pri_map, hex_map, poly_map);
      } else if (strstr(fname,".cgns") != NULL)
      {
#ifdef HAVE_CGNS
        P_CGNS_write(1, pname, pnn, pnode, nb, b_name, pnt, ptri, pnq, pquad, png, pngon,
                    pntet, ptet, pnpyr, ppyr, pnpri, ppri, pnhex, phex, pnply, ppoly,
                    node_map, tri_map, quad_map, ngon_map, tet_map, pyr_map, pri_map, hex_map, poly_map);
#else
        fprintf(stderr,"\nCGNS not enabled!");
        fflush(stderr);
        exit(0);
#endif
      }

      nlist.Redimension(0);
      overlap.Redimension(0);

      for (b=0; b < nb; b++)
      {
        for (c=0; c < png[b]; c++)
        {
          free(pngon[b][c]);
          pngon[b][c] = 0;
        }
      }
      for (c=0; c < pnply; c++)
      {
        for (i=ppoly[c][0][0]; i >= 0; i--)
          free(ppoly[c][i]);
        free(ppoly[c]);
        ppoly[c] = 0;
      }
    }

    // clean up memory

    for (b=0; b < nb; b++)
    {
      if (nt[b] > 0)
      {
        for (c=0; c < nt[b]; c++)
        {
          free(ptri[b][c]);
        }
        free(ptri[b]);
        free(tri_map[b]);
      }
      if (nq[b] > 0)
      {
        for (c=0; c < nq[b]; c++)
        {
          free(pquad[b][c]);
        }
        free(pquad[b]);
        free(quad_map[b]);
      }
      if (ngon[b] > 0)
      {
        for (c=0; c < ngon[b]; c++)
        {
          free(pngon[b][c]);
        }
        free(pngon[b]);
        free(ngon_map[b]);
      }
    }
    free(ptri);
    free(tri_map);
    free(pquad);
    free(quad_map);
    free(pngon);
    free(ngon_map);
    free(pnt);
    free(pnq);
    free(png);
    if (ntet > 0)
    {
      for (c=0; c < ntet; c++)
      {
        free(ptet[c]);
      }
      free(ptet);
      free(ptet_vc);
      free(tet_map);
    }
    if (npyr > 0)
    {
      for (c=0; c < npyr; c++)
      {
        free(ppyr[c]);
      }
      free(ppyr);
      free(ppyr_vc);
      free(pyr_map);
    }
    if (npri > 0)
    {
      for (c=0; c < npri; c++)
      {
        free(ppri[c]);
      }
      free(ppri);
      free(ppri_vc);
      free(pri_map);
    }
    if (nhex > 0)
    {
      for (c=0; c < nhex; c++)
      {
        free(phex[c]);
      }
      free(phex);
      free(phex_vc);
      free(hex_map);
    }
    if (nply > 0)
    {
      free(ppoly);
      free(ppoly_vc);
      free(poly_map);
    }

    for (p=0; p < nparts; p++)
      delete plist[p];
    delete[] plist;

    for (n=0; n < nn; n++)
      free(node_map[n]);
    free(node_map);
    free(pnode);
    delete[] ntmp;
    delete[] index;
  }

  // free original memory
  for (b=0; b < nb; b++)
  {
    if (nt[b] > 0)
    {
      for (c=0; c < nt[b]; c++)
      {
        free(tri_conn[b][c]);
      }
      free(tri_conn[b]);
    }
    if (nq[b] > 0)
    {
      for (c=0; c < nq[b]; c++)
      {
        free(quad_conn[b][c]);
      }
      free(quad_conn[b]);
    }
    if (ngon[b] > 0)
    {
      for (c=0; c < ngon[b]; c++)
      {
        free(ngon_conn[b][c]);
      }
      free(ngon_conn[b]);
    }
  }
  free(tri_conn);
  free(quad_conn);
  free(ngon_conn);
  free(nt);
  free(nq);
  free(ngon);
  if (ntet > 0)
  {
    for (c=0; c < ntet; c++)
    {
      free(tet_conn[c]);
    }
    free(tet_conn);
    free(tet_vc);
  }
  if (npyr > 0)
  {
    for (c=0; c < npyr; c++)
    {
      free(pyr_conn[c]);
    }
    free(pyr_conn);
    free(pyr_vc);
  }
  if (npri > 0)
  {
    for (c=0; c < npri; c++)
    {
      free(pri_conn[c]);
    }
    free(pri_conn);
    free(pri_vc);
  }
  if (nhex > 0)
  {
    for (c=0; c < nhex; c++)
    {
      free(hex_conn[c]);
    }
    free(hex_conn);
    free(hex_vc);
  }
  if (nply > 0)
  {
    for (c=0; c < nply; c++)
    {
      for (i=poly_conn[c][0][0]; i >= 0; i--)
        free(poly_conn[c][i]);
      free(poly_conn[c]);
    }
    free(poly_conn);
    free(poly_vc);
  }

  delete[] part;

  return(error);
}

int recomp(char *fname, int nparts)
{
  int b, c, i, j, n, nn, nb, ntet, npyr, npri, nhex, nply, nvol;
  int *nt, *nq, *ngon, ***tri_conn, ***quad_conn, ***ngon_conn;
  int **tet_conn, **pyr_conn, **pri_conn, **hex_conn, ***poly_conn;
  int **node_map, **tri_map, **quad_map, **ngon_map, *tet_map, *pyr_map, *pri_map, *hex_map, *poly_map;
  int *tet_vc, *pyr_vc, *pri_vc, *hex_vc, *poly_vc;
  int pnn, pnb, pntet, pnpyr, pnpri, pnhex, pnply, pnvol;
  int *pnt, *pnq, *png, ***ptri, ***pquad, ***pngon, **ptet, **ppyr, **ppri, **phex, ***ppoly;
  int *ptet_vc, *ppyr_vc, *ppri_vc, *phex_vc, *ppoly_vc;
  char **b_name, **pb_name, **vol_name, **pvol_name;
  Point *node, *pnode;
  const int bdim = 132;
  char buff[bdim], pname[bdim];
  int error = 0;

  int p;

  // set all counts to 0
  nvol = nn = nb = ntet = npyr = npri = nhex = nply = 0;
  // set all pointers to null
  b_name = 0;
  vol_name = 0;
  node = 0;
  nt = 0;
  nq = 0;
  ngon = 0;
  tri_conn = 0;
  quad_conn = 0;
  ngon_conn = 0;
  tet_conn = 0;
  pyr_conn = 0;
  pri_conn = 0;
  hex_conn = 0;
  poly_conn = 0;
  tet_vc = 0;
  pyr_vc = 0;
  pri_vc = 0;
  hex_vc = 0;
  poly_vc = 0;

  // initialize parallel arrays
  pnvol = pnn = pnb = pntet = pnpyr = pnpri = pnhex = pnply = 0;
  pb_name = 0;
  pvol_name = 0;
  pnode = 0;
  pnt = 0;
  pnq = 0;
  png = 0;
  ptri = 0;
  pquad = 0;
  pngon = 0;
  ptet = 0;
  ppyr = 0;
  ppri = 0;
  phex = 0;
  ppoly = 0;
  node_map = 0;
  tri_map = 0;
  quad_map = 0;
  ngon_map = 0;
  tet_map = 0;
  pyr_map = 0;
  pri_map = 0;
  hex_map = 0;
  poly_map = 0;

  // read parallel files one at a time to determine totals
  for (p=0; p < nparts; p++)
  {
    // create the parallel filename
    sprintf(buff,"%s",fname);
    if (strstr(fname,".sg") != NULL)
    {
      char *ptr = strstr(buff,".sg");
      if (ptr == NULL)
      {
        printf("\nSG suffix <.sg> not found in file name!");
        fflush(stdout);
        exit(0);
      } else
        *ptr = '\0';
      sprintf(pname,"%s_%d.sg",buff,p);
    } else if (strstr(fname,".cgns") != NULL)
    {
      char *ptr = strstr(buff,".cgns");
      if (ptr == NULL)
      {
        printf("\nCGNS suffix <.cgns> not found in file name!");
        fflush(stdout);
        exit(0);
      } else
        *ptr = '\0';
      sprintf(pname,"%s_%d.cgns",buff,p);
    }

    // read parallel mesh file
    if (strstr(fname,".sg") != NULL)
    {
      if ((error = P_SimGrid_read(1, pname, pnn, &pnode, pnb, &pb_name, pnvol, &pvol_name, &pnt, &ptri, &pnq, &pquad, &png, &pngon,
                  pntet, &ptet, &ptet_vc, pnpyr, &ppyr, &ppyr_vc, pnpri, &ppri, &ppri_vc, pnhex, &phex, &phex_vc, pnply, &ppoly, &ppoly_vc,
                  &node_map, &tri_map, &quad_map, &ngon_map, &tet_map, &pyr_map, &pri_map, &hex_map, &poly_map)) != 0)
      {
        fprintf(stderr,"\nError reading parallel SG file <%s>.",pname);
        fflush(stderr);
        return(error);
        //exit(0);
      }
    } else if (strstr(fname,".cgns") != NULL)
    {
#ifdef HAVE_CGNS
      if ((error = P_CGNS_read(1, pname, pnn, &pnode, pnb, &pb_name, &pnt, &ptri, &pnq, &pquad, &png, &pngon,
                  pntet, &ptet, pnpyr, &ppyr, pnpri, &ppri, pnhex, &phex, pnply, &ppoly,
                  &node_map, &tri_map, &quad_map, &ngon_map, &tet_map, &pyr_map, &pri_map, &hex_map, &poly_map)) != 0)
      {
        fprintf(stderr,"\nError reading parallel CGNS file <%s>.",pname);
        fflush(stderr);
        return(error);
        //exit(0);
      }
#else
        fprintf(stderr,"\nCGNS not enabled!");
        fflush(stderr);
        exit(0);
#endif
      if (pntet > 0)
        ptet_vc = (int*)calloc(pntet,sizeof(int));
      if (pnpyr > 0)
        ppyr_vc = (int*)calloc(pnpyr,sizeof(int));
      if (pnpri > 0)
        ppri_vc = (int*)calloc(pnpri,sizeof(int));
      if (pnhex > 0)
        phex_vc = (int*)calloc(pnhex,sizeof(int));
      if (pnply > 0)
        ppoly_vc = (int*)calloc(pnply,sizeof(int));
    }
  
    printf("\nFile <%s> information:",pname);

    printf("\n  Number of nodes = %d",pnn);
    for (n=0; n < pnn; n++)
      nn = MAX(nn,node_map[n][0]+1);

    printf("\n  Number of volumes = %i",pnvol);
    if (nvol != pnvol)
    {
      if (vol_name == 0)
      {
        vol_name = (char**)malloc(MAX(nvol,pnvol)*sizeof(char*));
        for (i=0; i < MAX(nvol,pnvol); i++)
          vol_name[i] = (char*)malloc(33*sizeof(char));
      } else if (pnvol > nvol)
      {
        vol_name = (char**)realloc((void*)vol_name,pnvol*sizeof(char*));
        for (i=nvol; i < pnvol; i++)
          vol_name[i] = (char*)malloc(33*sizeof(char));
      }
      nvol = MAX(nvol,pnvol);
    }

    printf("\n  Number of boundaries = %d",pnb);
    if (nb != pnb)
    {
      if (b_name == 0)
      {
        b_name = (char**)malloc(MAX(nb,pnb)*sizeof(char*));
        for (i=0; i < MAX(nb,pnb); i++)
          b_name[i] = (char*)malloc(33*sizeof(char));
      } else if (pnb > nb)
      {
        b_name = (char**)realloc((void*)b_name,pnb*sizeof(char*));
        for (i=nb; i < pnb; i++)
          b_name[i] = (char*)malloc(33*sizeof(char));
      }
      if (nt == 0)
      {
        nt = (int*)malloc(MAX(nb,pnb)*sizeof(int));
        for (i=0; i < MAX(nb,pnb); i++)
          nt[i] = 0;
      } else if (pnb > nb)
      {
        nt = (int*)realloc((void*)nt,pnb*sizeof(int));
        for (i=nb; i < pnb; i++)
          nt[i] = 0;
      }
      if (nq == 0)
      {
        nq = (int*)malloc(MAX(nb,pnb)*sizeof(int));
        for (i=0; i < MAX(nb,pnb); i++)
          nq[i] = 0;
      } else if (pnb > nb)
      {
        nq = (int*)realloc((void*)nq,pnb*sizeof(int));
        for (i=nb; i < pnb; i++)
          nq[i] = 0;
      }
      if (ngon == 0)
      {
        ngon = (int*)malloc(MAX(nb,pnb)*sizeof(int));
        for (i=0; i < MAX(nb,pnb); i++)
          ngon[i] = 0;
      } else if (pnb > nb)
      {
        ngon = (int*)realloc((void*)ngon,pnb*sizeof(int));
        for (i=nb; i < pnb; i++)
          ngon[i] = 0;
      }
      nb = MAX(nb,pnb);
    }
      
    for (i=0; i < pnb; i++)
    {
      if (pnt[i] > 0)
        printf("\nBoundary %d, number of triangles      = %d",i+1,pnt[i]);
      for (j=0; j < pnt[i]; j++)
        nt[i] = MAX(nt[i],tri_map[i][j]+1);
      if (pnq[i] > 0)
        printf("\nBoundary %d, number of quadrilaterals = %d",i+1,pnq[i]);
      for (j=0; j < pnq[i]; j++)
        nq[i] = MAX(nq[i],quad_map[i][j]+1);
      if (png[i] > 0)
        printf("\nBoundary %d, number of polygons       = %d",i+1,png[i]);
      for (j=0; j < png[i]; j++)
        ngon[i] = MAX(ngon[i],ngon_map[i][j]+1);
    }
    if (pntet > 0)
    {
      printf("\n  Number of tetrahedra = %d",pntet);
      for (n=0; n < pntet; n++)
        ntet = MAX(ntet,tet_map[n]+1);
    }
    if (pnpyr > 0)
    {
      printf("\n  Number of pyramids   = %d",pnpyr);
      for (n=0; n < pnpyr; n++)
        npyr = MAX(npyr,pyr_map[n]+1);
    }
    if (pnpri > 0)
    {
      printf("\n  Number of prisms     = %d",pnpri);
      for (n=0; n < pnpri; n++)
        npri = MAX(npri,pri_map[n]+1);
    }
    if (pnhex > 0)
    {
      printf("\n  Number of hexahedra  = %d",pnhex);
      for (n=0; n < pnhex; n++)
        nhex = MAX(nhex,hex_map[n]+1);
    }
    if (pnply > 0)
    {
      printf("\n  Number of polyhedra  = %d",pnply);
      for (n=0; n < pnply; n++)
        nply = MAX(nply,poly_map[n]+1);
    }

    for (n=0; n < pnn; n++)
      free(node_map[n]);
    free(node_map);
    free(pnode);
    for (i=0; i < pnvol; i++)
      free(pvol_name[i]);
    free(pvol_name);
    for (i=0; i < pnb; i++)
    {
      if (pnt[i] > 0)
      {
        for (j=0; j < pnt[i]; j++)
          free(ptri[i][j]);
        free(ptri[i]);
        free(tri_map[i]);
      }
      if (pnq[i] > 0)
      {
        for (j=0; j < pnq[i]; j++)
          free(pquad[i][j]);
        free(pquad[i]);
        free(quad_map[i]);
      }
      if (png[i] > 0)
      {
        for (j=0; j < png[i]; j++)
          free(pngon[i][j]);
        free(pngon[i]);
        free(ngon_map[i]);
      }
      free(pb_name[i]);
    }
    free(pb_name);
    free(pnt);
    free(pnq);
    free(png);
    free(ptri);
    free(pquad);
    free(pngon);
    free(tri_map);
    free(quad_map);
    free(ngon_map);

    if (pntet > 0)
    {
      for (n=0; n < pntet; n++)
        free(ptet[n]);
      free(ptet);
      free(ptet_vc);
      free(tet_map);
    }
    if (pnpyr > 0)
    {
      for (n=0; n < pnpyr; n++)
        free(ppyr[n]);
      free(ppyr);
      free(ppyr_vc);
      free(pyr_map);
    }
    if (pnpri > 0)
    {
      for (n=0; n < pnpri; n++)
        free(ppri[n]);
      free(ppri);
      free(ppri_vc);
      free(pri_map);
    }
    if (pnhex > 0)
    {
      for (n=0; n < pnhex; n++)
        free(phex[n]);
      free(phex);
      free(phex_vc);
      free(hex_map);
    }
    if (pnply > 0)
    {
      for (n=0; n < pnply; n++)
      {
        for (i=ppoly[n][0][0]; i >= 0; i--)
          free(ppoly[n][i]);
        free(ppoly[n]);
      }
      free(ppoly);
      free(ppoly_vc);
      free(poly_map);
    }
    pnvol = pnn = pnb = pntet = pnpyr = pnpri = pnhex = pnply = 0;
    pb_name = 0;
    pvol_name = 0;
    pnode = 0;
    pnt = 0;
    pnq = 0;
    png = 0;
    ptri = 0;
    pquad = 0;
    pngon = 0;
    ptet = 0;
    ppyr = 0;
    ppri = 0;
    phex = 0;
    ppoly = 0;
    ptet_vc = 0;
    ppyr_vc = 0;
    ppri_vc = 0;
    phex_vc = 0;
    ppoly_vc = 0;
    node_map = 0;
    tri_map = 0;
    quad_map = 0;
    ngon_map = 0;
    tet_map = 0;
    pyr_map = 0;
    pri_map = 0;
    hex_map = 0;
    poly_map = 0;
  }
  printf("\n\n");

  printf("\nSerial file <%s> information:",fname);

  // allocate global arrays
  printf("\nTotal number of nodes = %d",nn);
  node = (Point*)malloc(nn*sizeof(Point));
  printf("\nTotal number of boundaries = %d",nb);
  tri_conn = (int***)malloc(nb*sizeof(int**));
  quad_conn = (int***)malloc(nb*sizeof(int**));
  ngon_conn = (int***)malloc(nb*sizeof(int**));
  for (i=0; i < nb; i++)
  {
    if (nt[i] > 0)
      printf("\nBoundary %d, total number of triangles      = %d",i+1,nt[i]);
    if (nt[i] == 0)
      tri_conn[i] = 0;
    else
    {
      tri_conn[i] = (int**)malloc(nt[i]*sizeof(int*));
      for (j=0; j < nt[i]; j++)
        tri_conn[i][j] = (int*)malloc(3*sizeof(int));
    }
    if (nq[i] > 0)
      printf("\nBoundary %d, total number of quadrilaterals = %d",i+1,nq[i]);
    if (nq[i] == 0)
      quad_conn[i] = 0;
    else
    {
      quad_conn[i] = (int**)malloc(nq[i]*sizeof(int*));
      for (j=0; j < nq[i]; j++)
        quad_conn[i][j] = (int*)malloc(4*sizeof(int));
    }
    if (ngon[i] > 0)
      printf("\nBoundary %d, total number of polygons       = %d",i+1,ngon[i]);
    if (ngon[i] == 0)
      ngon_conn[i] = 0;
    else
    {
      ngon_conn[i] = (int**)malloc(ngon[i]*sizeof(int*));
      for (j=0; j < ngon[i]; j++)
        ngon_conn[i][j] = 0;  // allocated per polygon later
    }
  }
  if (ntet > 0)
  {
    printf("\nTotal number of tetrahedra = %d",ntet);
    tet_vc = (int*)malloc(ntet*sizeof(int));
    tet_conn = (int**)malloc(ntet*sizeof(int*));
    for (i=0; i < ntet; i++)
      tet_conn[i] = (int*)malloc(4*sizeof(int));
  }
  if (npyr > 0)
  {
    printf("\nTotal number of pyramids   = %d",npyr);
    pyr_vc = (int*)malloc(npyr*sizeof(int));
    pyr_conn = (int**)malloc(npyr*sizeof(int*));
    for (i=0; i < npyr; i++)
      pyr_conn[i] = (int*)malloc(5*sizeof(int));
  }
  if (npri > 0)
  {
    printf("\nTotal number of prisms     = %d",npri);
    pri_vc = (int*)malloc(npri*sizeof(int));
    pri_conn = (int**)malloc(npri*sizeof(int*));
    for (i=0; i < npri; i++)
      pri_conn[i] = (int*)malloc(6*sizeof(int));
  }
  if (nhex > 0)
  {
    printf("\nTotal number of hexahedra  = %d",nhex);
    hex_vc = (int*)malloc(nhex*sizeof(int));
    hex_conn = (int**)malloc(nhex*sizeof(int*));
    for (i=0; i < nhex; i++)
      hex_conn[i] = (int*)malloc(8*sizeof(int));
  }
  if (nply > 0)
  {
    printf("\nTotal number of polyhedra  = %d",nply);
    poly_vc = (int*)malloc(nply*sizeof(int));
    poly_conn = (int***)malloc(nply*sizeof(int*));
    for (i=0; i < nply; i++)
      poly_conn[i] = 0; // allocated per element later
  }
  
  // reread parallel files and store data in global arrays
  for (p=0; p < nparts; p++)
  {
    // create the parallel filename
    sprintf(buff,"%s",fname);
    if (strstr(fname,".sg") != NULL)
    {
      char *ptr = strstr(buff,".sg");
      if (ptr == NULL)
      {
        printf("\nSG suffix <.sg> not found in file name!");
        fflush(stdout);
        exit(0);
      } else
        *ptr = '\0';
      sprintf(pname,"%s_%d.sg",buff,p);
    } else if (strstr(fname,".cgns") != NULL)
    {
      char *ptr = strstr(buff,".cgns");
      if (ptr == NULL)
      {
        printf("\nCGNS suffix <.cgns> not found in file name!");
        fflush(stdout);
        exit(0);
      } else
        *ptr = '\0';
      sprintf(pname,"%s_%d.cgns",buff,p);
    }

    // read serial mesh file
    if (strstr(fname,".sg") != NULL)
    {
      if ((error = P_SimGrid_read(1, pname, pnn, &pnode, pnb, &pb_name, pnvol, &pvol_name, &pnt, &ptri, &pnq, &pquad, &png, &pngon,
                  pntet, &ptet, &ptet_vc, pnpyr, &ppyr, &ppyr_vc, pnpri, &ppri, &ppri_vc, pnhex, &phex, &phex_vc, pnply, &ppoly, &ppoly_vc,
                  &node_map, &tri_map, &quad_map, &ngon_map, &tet_map, &pyr_map, &pri_map, &hex_map, &poly_map)) != 0)
      {
        fprintf(stderr,"\nError reading parallel SG file <%s>.",pname);
        fflush(stderr);
        return(error);
        //exit(0);
      }
    } else if (strstr(fname,".cgns") != NULL)
    {
#ifdef HAVE_CGNS
      if ((error = P_CGNS_read(1, pname, pnn, &pnode, pnb, &pb_name, &pnt, &ptri, &pnq, &pquad, &png, &pngon,
                  pntet, &ptet, pnpyr, &ppyr, pnpri, &ppri, pnhex, &phex, pnply, &ppoly,
                  &node_map, &tri_map, &quad_map, &ngon_map, &tet_map, &pyr_map, &pri_map, &hex_map, &poly_map)) != 0)
      {
        fprintf(stderr,"\nError reading parallel CGNS file <%s>.",pname);
        fflush(stderr);
        return(error);
        //exit(0);
      }
#else
        fprintf(stderr,"\nCGNS not enabled!");
        fflush(stderr);
        exit(0);
#endif
      if (pntet > 0)
        ptet_vc = (int*)calloc(pntet,sizeof(int));
      if (pnpyr > 0)
        ppyr_vc = (int*)calloc(pnpyr,sizeof(int));
      if (pnpri > 0)
        ppri_vc = (int*)calloc(pnpri,sizeof(int));
      if (pnhex > 0)
        phex_vc = (int*)calloc(pnhex,sizeof(int));
      if (pnply > 0)
        ppoly_vc = (int*)calloc(pnply,sizeof(int));
    }
  
    for (n=0; n < pnn; n++)
      node[node_map[n][0]] = pnode[n];

    for (b=0; b < pnvol; b++)
    {
      strcpy(vol_name[b],pvol_name[b]);
    }
    for (b=0; b < pnb; b++)
    {
      strcpy(b_name[b],pb_name[b]);
      for (i=0; i < pnt[b]; i++)
        for (j=0; j < 3; j++)
          tri_conn[b][tri_map[b][i]][j] = node_map[ptri[b][i][j]][0];
      for (i=0; i < pnq[b]; i++)
        for (j=0; j < 4; j++)
          quad_conn[b][quad_map[b][i]][j] = node_map[pquad[b][i][j]][0];
      for (i=0; i < png[b]; i++)
      {
        ngon_conn[b][ngon_map[b][i]] = (int*)realloc((void*)ngon_conn[b][ngon_map[b][i]],(pngon[b][i][0]+1)*sizeof(int));
        ngon_conn[b][ngon_map[b][i]][0] = pngon[b][i][0];
        for (j=1; j <= pngon[b][i][0]; j++)
          ngon_conn[b][ngon_map[b][i]][j] = node_map[pngon[b][i][j]][0];
      }
    }

    if (pntet > 0)
      for (c=0; c < pntet; c++)
      {
        for (i=0; i < 4; i++)
          tet_conn[tet_map[c]][i] = node_map[ptet[c][i]][0];
        tet_vc[tet_map[c]] = ptet_vc[c];
      }

    if (pnpyr > 0)
      for (c=0; c < pnpyr; c++)
      {
        for (i=0; i < 5; i++)
          pyr_conn[pyr_map[c]][i] = node_map[ppyr[c][i]][0];
        pyr_vc[pyr_map[c]] = ppyr_vc[c];
      }

    if (pnpri > 0)
      for (c=0; c < pnpri; c++)
      {
        for (i=0; i < 6; i++)
          pri_conn[pri_map[c]][i] = node_map[ppri[c][i]][0];
        pri_vc[pri_map[c]] = ppri_vc[c];
      }

    if (pnhex > 0)
      for (c=0; c < pnhex; c++)
      {
        for (i=0; i < 8; i++)
          hex_conn[hex_map[c]][i] = node_map[phex[c][i]][0];
        hex_vc[hex_map[c]] = phex_vc[c];
      }

    if (pnply > 0)
    {
      for (c=0; c < pnply; c++)
      {
        poly_conn[poly_map[c]] = (int**)realloc((void*)poly_conn[poly_map[c]],(ppoly[c][0][0]+1)*sizeof(int*));
        poly_conn[poly_map[c]][0] = (int*)malloc(sizeof(int));
        poly_conn[poly_map[c]][0][0] = ppoly[c][0][0];
        for (i=1; i <= ppoly[c][0][0]; i++)
        {
          poly_conn[poly_map[c]][i] = (int*)malloc((ppoly[c][i][0]+1)*sizeof(int));
          poly_conn[poly_map[c]][i][0] = ppoly[c][i][0];
          for (j=1; j <= ppoly[c][i][0]; j++)
            poly_conn[poly_map[c]][i][j] = node_map[ppoly[c][i][j]][0];
        }
        poly_vc[poly_map[c]] = ppoly_vc[c];
      }
    }
    for (n=0; n < pnn; n++)
      free(node_map[n]);
    free(node_map);
    free(pnode);
    for (i=0; i < pnb; i++)
    {
      if (pnt[i] > 0)
      {
        for (j=0; j < pnt[i]; j++)
          free(ptri[i][j]);
        free(ptri[i]);
        free(tri_map[i]);
      }
      if (pnq[i] > 0)
      {
        for (j=0; j < pnq[i]; j++)
          free(pquad[i][j]);
        free(pquad[i]);
        free(quad_map[i]);
      }
      if (png[i] > 0)
      {
        for (j=0; j < png[i]; j++)
          free(pngon[i][j]);
        free(pngon[i]);
        free(ngon_map[i]);
      }
      free(pb_name[i]);
    }
    free(pb_name);
    free(pnt);
    free(pnq);
    free(png);
    free(ptri);
    free(pquad);
    free(pngon);
    free(tri_map);
    free(quad_map);
    free(ngon_map);

    if (pntet > 0)
    {
      for (n=0; n < pntet; n++)
        free(ptet[n]);
      free(ptet);
      free(ptet_vc);
      free(tet_map);
    }
    if (pnpyr > 0)
    {
      for (n=0; n < pnpyr; n++)
        free(ppyr[n]);
      free(ppyr);
      free(ppyr_vc);
      free(pyr_map);
    }
    if (pnpri > 0)
    {
      for (n=0; n < pnpri; n++)
        free(ppri[n]);
      free(ppri);
      free(ppri_vc);
      free(pri_map);
    }
    if (pnhex > 0)
    {
      for (n=0; n < pnhex; n++)
        free(phex[n]);
      free(phex);
      free(phex_vc);
      free(hex_map);
    }
    if (pnply > 0)
    {
      for (n=0; n < pnply; n++)
      {
        for (i=ppoly[n][0][0]; i >= 0; i--)
          free(ppoly[n][i]);
        free(ppoly[n]);
      }
      free(ppoly);
      free(ppoly_vc);
      free(poly_map);
    }
    pnvol = pnn = pnb = pntet = pnpyr = pnpri = pnhex = pnply = 0;
    pb_name = 0;
    pvol_name = 0;
    pnode = 0;
    pnt = 0;
    pnq = 0;
    png = 0;
    ptri = 0;
    pquad = 0;
    pngon = 0;
    ptet = 0;
    ppyr = 0;
    ppri = 0;
    phex = 0;
    ppoly = 0;
    ptet_vc = 0;
    ppyr_vc = 0;
    ppri_vc = 0;
    phex_vc = 0;
    ppoly_vc = 0;
    node_map = 0;
    tri_map = 0;
    quad_map = 0;
    ngon_map = 0;
    tet_map = 0;
    pyr_map = 0;
    pri_map = 0;
    hex_map = 0;
    poly_map = 0;
  }

  printf("\nWriting %s",fname);
  fflush(stdout);

  // write serial mesh file
  if (strstr(fname,".sg") != NULL)
  {
    SimGrid_write(fname, nn, node, nb, b_name, nvol, vol_name, nt, tri_conn, nq, quad_conn, ngon, ngon_conn,
                ntet, tet_conn, tet_vc, npyr, pyr_conn, pyr_vc, npri, pri_conn, pri_vc, nhex, hex_conn, hex_vc, nply, poly_conn, poly_vc);
  } else if (strstr(fname,".cgns") != NULL)
  {
#ifdef HAVE_CGNS
    CGNS_write(fname, nn, node, nb, b_name, nt, tri_conn, nq, quad_conn, ngon, ngon_conn,
                ntet, tet_conn, npyr, pyr_conn, npri, pri_conn, nhex, hex_conn, nply, poly_conn);
#else
    fprintf(stderr,"\nCGNS not enabled!");
    fflush(stderr);
    exit(0);
#endif
  }

  // free up memory
  free(node);
  for (b=0; b < nb; b++)
  {
    if (nt[b] > 0)
    {
      for (i=0; i < nt[b]; i++)
        free(tri_conn[b][i]);
      free(tri_conn[b]);
    }
    if (nq[b] > 0)
    {
      for (i=0; i < nq[b]; i++)
        free(quad_conn[b][i]);
      free(quad_conn[b]);
    }
    if (ngon[b] > 0)
    {
      for (i=0; i < ngon[b]; i++)
        free(ngon_conn[b][i]);
      free(ngon_conn[b]);
    }
    free(b_name[b]);
  }
  free(b_name);
  free(nt);
  free(nq);
  free(ngon);
  free(tri_conn);
  free(quad_conn);
  free(ngon_conn);
  if (ntet > 0)
  {
    for (i=0; i < ntet; i++)
      free(tet_conn[i]);
    free(tet_conn);
    free(tet_vc);
  }
  if (npyr > 0)
  {
    for (i=0; i < npyr; i++)
      free(pyr_conn[i]);
    free(pyr_conn);
    free(pyr_vc);
  }
  if (npri > 0)
  {
    for (i=0; i < npri; i++)
      free(pri_conn[i]);
    free(pri_conn);
    free(pri_vc);
  }
  if (nhex > 0)
  {
    for (i=0; i < nhex; i++)
      free(hex_conn[i]);
    free(hex_conn);
    free(hex_vc);
  }
  if (nply > 0)
  {
    for (i=0; i < nply; i++)
    {
      for (j=poly_conn[i][0][0]; j >= 0; j--)
        free(poly_conn[i][j]);
      free(poly_conn[i]);
    }
    free(poly_conn);
    free(poly_vc);
  }

  return(error);

}

int main(int argcs, char* pArgs[])
{
  const int bdim = 132;
  char buff[bdim];
  char cname[bdim], sname[bdim], mname[bdim];
  FILE *in_f, *jou_f, *f_ptr;
  int volchk, wind, pflag;
  int error;
  double tol;

  printf("\n==========================================================================");
  printf("\n COPYRIGHT 2003-2012 THE UNIVERSITY OF TENNESSEE AT CHATTANOOGA           ");
  printf("\n                                                                          ");
  printf("\n                    RIGHTS IN DATA                                        ");
  printf("\n                                                                          ");
  printf("\n THIS SOFTWARE IS SUBMITTED WITH RESTRICTED RIGHTS UNDER GOVERNMENT       ");
  printf("\n    CONTRACTS. USE, REPRODUCTION, OR DISCLOSURE IS SUBJECT TO             ");
  printf("\n         RESTRICTIONS SET FORTH IN THESE CONTRACTS AND FEDERAL            ");
  printf("\n              RIGHTS IN DATA CONTRACT CLAUSES.                            ");
  printf("\n       ALL RIGHTS NOT RESERVED FOR THE GOVERNMENT ARE RETAINED BY         ");
  printf("\n              THE UNIVERSITY OF TENNESSEE AT CHATTANOOGA                  ");
  printf("\n                                                                          ");
  printf("\n Grid Conversion Utility (Conv)                                           ");
  printf("\n NOTE: This data includes the SimCenter convert code,                     ");
  printf("\n which was developed under private non-government funding.                ");
  printf("\n This software is submitted with limited rights to use, reproduce,        ");
  printf("\n and disclose this data for Government Purposes only.                     ");
  printf("\n Requests for access to the software for non-governmental purposes        ");
  printf("\n should be referrred to                                                   "); 
  printf("\n                                                                          ");
  printf("\n    Dr. Steve Karman                                                      "); 
  printf("\n    Steve-Karman@utc.edu                                                  "); 
  printf("\n    423-425-5492  or  423-425-5470                                        "); 
  printf("\n                                                                          ");
  printf("\n    SimCenter: National Center for Computational Engineering              "); 
  printf("\n    701 East M. L. King Boulevard                                         "); 
  printf("\n    Chattanooga, TN 37403                                                 "); 
  printf("\n==========================================================================\n");
  time_t tm;
  char *t_char;
  time(&tm);
  t_char = ctime(&tm);
  printf("\nRun started at %s",t_char);

  if ((jou_f=fopen("Conv.jou","w")) == NULL)
  {
    printf("\nCouldn't open file journal file");
    exit(0);
  }

  if (--argcs < 1)
  {
    printf("\nNo input file specified!");
    printf("\nUsing standard input!");
    in_f = stdin;
  } else if ((in_f=fopen(pArgs[argcs],"r")) == NULL)
  {
    printf("\nCouldn't open file <%s>\n",pArgs[argcs]);
    fprintf(jou_f,"\nCouldn't open file <%s>\n",pArgs[argcs]);
    fclose(jou_f);
    exit(0);
  }

  cname[0]='\0';
  sname[0]='\0';

  int mode, nparts, pmode, model;
  do
  {
    error = 0;

    printf("\n\nConv:");
    journal(in_f, stdout, jou_f, "#Exit            - 0\n"
                                 "#Convert file    - 1\n"
                                 "#Decomp file     - 2\n"
                                 "#Recomp files    - 3\n"
                                 "#Merge 2 files   - 4\n"
                                 "#Set volume name - 5\n"
                                 "#Choice >",mode);
    switch (mode)
    {
      case 1:
        printf("\nSupported File Formats                Input/Output  File Extension");
        printf("\n------------------------------------- ------------  --------------");
        printf("\nFieldview ASCII double precision         both       .crunch or .grd");
        printf("\nFieldview binary                         both       .uns or .fv");
        printf("\nNastran                                  both       .nas");
        printf("\nStarCD                                   both       .inp");
        printf("\nEnsight                                 output      .case");
#ifdef HAVE_CGNS
        printf("\nCGNS binary                              both       .cgns");
#endif
        printf("\nSimGrid binary                           both       .sg");
        printf("\nGeneric mesh ASCII double precision      both       .mesh3D");
        printf("\nCART3D boundary file                    output      .tri");
        printf("\nXpatch (FACET) boundary file            output      .facet\n");
        journal(in_f, stdout, jou_f, "#Enter input grid file name >",cname);
        journal(in_f, stdout, jou_f, "#Enter output grid file name >",sname);
        volchk;
        journal(in_f, stdout, jou_f, "#Check volume element winding? [0,1] >",volchk);
        wind;
        journal(in_f, stdout, jou_f, "#Check boundary element winding? [0,1] >",wind);
        pflag;
        journal(in_f, stdout, jou_f, "#Convert polyhedra to standard? [0,1] >",pflag);

        mname[0]='\0'; // don't set volume name

        error = convert(cname,sname,wind,volchk,pflag,mname);

        break;
      case 2:
        journal(in_f, stdout, jou_f, "#Enter input serial file name >",cname);
        journal(in_f, stdout, jou_f, "#Create partition file [ 1 ], apply partition file [ 2 ], both [ 3 ] >",pmode);
        if (pmode == 1 || pmode == 3)
        {
          journal(in_f, stdout, jou_f, "#Enter number of partitions >",nparts);
          journal(in_f, stdout, jou_f, "#Partitioning model [0 - metis, 1 - split tree] >",model);
          error = decomp(cname,nparts,model,1);
        }
        if (pmode == 2 || pmode == 3)
        {
          nparts = model = -1;
          error = decomp(cname,nparts,model,2);
        }

        break;
      case 3:
        journal(in_f, stdout, jou_f, "#Enter output serial file name >",sname);
        int flag;
        nparts = 0;
        // determine how many partition files exist
        do
        {
          // create part file name
          sprintf(buff,"%s",sname);
          if (strstr(sname,".sg") != NULL)
          {
            char *ptr = strstr(buff,".sg");
            if (ptr == NULL)
            {
              printf("\nSG suffix <.sg> not found in file name!");
              fflush(stdout);
              exit(0);
            } else
              *ptr = '\0';
            sprintf(cname,"%s_%d.sg",buff,nparts);
            // attempt to open part file
            if ((f_ptr=fopen(cname,"r")) != NULL)
            {
              nparts++;
              fclose(f_ptr);
              flag = 1;
              printf("\nSG file <%s> found.",cname);
            } else
            {
              flag = 0;
              printf("\nSG file <%s> not found.",cname);
            }
          } else
          {
            char *ptr = strstr(buff,".cgns");
            if (ptr == NULL)
            {
              printf("\nCGNS suffix <.cgns> not found in file name!");
              fflush(stdout);
              exit(0);
            } else
              *ptr = '\0';
            sprintf(cname,"%s_%d.cgns",buff,nparts);
            // attempt to open part file
            if ((f_ptr=fopen(cname,"r")) != NULL)
            {
              nparts++;
              fclose(f_ptr);
              flag = 1;
              printf("\nCGNS file <%s> found.",cname);
            } else
            {
              flag = 0;
              printf("\nCGNS file <%s> not found.",cname);
            }
          }
          fflush(stdout);
        } while (flag);

        if (nparts > 0)
        {
          printf("\nNumber of part files detected = %d",nparts);
          fflush(stdout);

          error = recomp(sname,nparts);
        } else
          error = 1;

        break;
      case 4:
        printf("\nSupported Input File Formats          File Extension");
        printf("\n------------------------------------- --------------");
        printf("\nFieldview ASCII double precision      .crunch or .grd");
        printf("\nFieldview binary                      .uns or .fv");
#ifdef HAVE_CGNS
        printf("\nCGNS binary                           .cgns");
#endif
        printf("\nSimGrid binary                        .sg");
        printf("\nGeneric mesh ASCII double precision   .mesh3D");
        journal(in_f, stdout, jou_f, "#Enter 1st grid file name >",cname);
        journal(in_f, stdout, jou_f, "#Enter 2nd grid file name >",sname);
        journal(in_f, stdout, jou_f, "#Enter output grid file name >",mname);
        journal(in_f, stdout, jou_f, "#Enter tolerance >",tol);

        error = merge(cname,sname,mname,tol);

        break;
      case 5:
        printf("\nSupported File Formats                File Extension");
        printf("\n------------------------------------- --------------");
        printf("\nStarCD                                .inp");
        printf("\nSimGrid binary                        .sg");
        journal(in_f, stdout, jou_f, "#Enter input grid file name >",cname);
        volchk=0;
        wind=0;
        pflag=0;
        journal(in_f, stdout, jou_f, "Enter volume name >",mname);

        error = convert(cname,cname,wind,volchk,pflag,mname);

        break;
      default:
        error = mode = 0;
        break;
    }

  } while (error != 0 || mode > 0);

  fclose(jou_f);

  if (in_f != stdin) fclose(in_f);

  time(&tm);
  t_char = ctime(&tm);
  printf("\nRun completed at %s",t_char);

  return 0;
}

