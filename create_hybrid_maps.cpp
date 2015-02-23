#ifdef PARALLEL
#include "mpi.h"
#include "Pmap.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "List.h"
#include "create_hybrid_maps.h"

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))

void create_hybrid_maps(FILE *out_f, int num_procs, int my_rank, int nn, int nb, int ntet, int npyr, int npri, int nhex, int nply,
                        int **tet_n, int **pyr_n, int **pri_n, int **hex_n, int ***poly_n,
                        int *nt, int ***t_n, int *nq, int ***q_n, int *ngon, int ***ngon_n,
                        int **nmap, int **tri_map, int **quad_map, int **ngon_map,
                        int *tet_map, int *pyr_map, int *pri_map, int *hex_map, int *poly_map)
{
  // THIS ASSUMES THE NODE OWNERSHIP IS CORRECT!!!
  int b, e, i, j, n, o;
#ifdef PARALLEL
  int k, l, m, p;
  MPI_Status single_status;
  int nreq_s, nreq_r, sposition, rposition;
  int *sendcnt;
  int *recvcnt;
  MPI_Request *srequest;
  MPI_Request *rrequest;
  MPI_Status *statuses;
  int *bdim;
  char **sbuff;
  char **rbuff;
  List **pelems;
#endif

  // first identify owning processors for each element as
  // lowest processor from any of the nodes of the element
  for (b=0; b < nb; b++)
  {
    for (i=0; i < nt[b]; i++)
    {
      o = num_procs;
      for (j=0; j < 3; j++)
      {
        n = t_n[b][i][j];
        o = MIN(o,nmap[n][1]);
      }
      if (o == my_rank)
        tri_map[b][i] = i;
      else
        tri_map[b][i] = -(o+1);
    }
    for (i=0; i < nq[b]; i++)
    {
      o = num_procs;
      for (j=0; j < 4; j++)
      {
        n = q_n[b][i][j];
        o = MIN(o,nmap[n][1]);
      }
      if (o == my_rank)
        quad_map[b][i] = i;
      else
        quad_map[b][i] = -(o+1);
    }
    for (i=0; i < ngon[b]; i++)
    {
      o = num_procs;
      for (j=1; j <= ngon_n[b][i][0]; j++)
      {
        n = ngon_n[b][i][j];
        o = MIN(o,nmap[n][1]);
      }
      if (o == my_rank)
        ngon_map[b][i] = i;
      else
        ngon_map[b][i] = -(o+1);
    }
  }

  for (e=0; e < ntet; e++)
  {
    o = num_procs;
    for (i=0; i < 4; i++)
    {
      n = tet_n[e][i];
      o = MIN(o,nmap[n][1]);
    }
    if (o == my_rank)
      tet_map[e] = e;
    else
      tet_map[e] = -(o+1);
  }
  for (e=0; e < npyr; e++)
  {
    o = num_procs;
    for (i=0; i < 5; i++)
    {
      n = pyr_n[e][i];
      o = MIN(o,nmap[n][1]);
    }
    if (o == my_rank)
      pyr_map[e] = e;
    else
      pyr_map[e] = -(o+1);
  }
  for (e=0; e < npri; e++)
  {
    o = num_procs;
    for (i=0; i < 6; i++)
    {
      n = pri_n[e][i];
      o = MIN(o,nmap[n][1]);
    }
    if (o == my_rank)
      pri_map[e] = e;
    else
      pri_map[e] = -(o+1);
  }
  for (e=0; e < nhex; e++)
  {
    o = num_procs;
    for (i=0; i < 8; i++)
    {
      n = hex_n[e][i];
      o = MIN(o,nmap[n][1]);
    }
    if (o == my_rank)
      hex_map[e] = e;
    else
      hex_map[e] = -(o+1);
  }
  for (e=0; e < nply; e++)
  {
    o = num_procs;
    for (i=1; i <= poly_n[e][0][0]; i++)
      for (j=1; j <= poly_n[e][i][0]; j++)
        o = MIN(o,nmap[poly_n[e][i][j]][1]);
    if (o == my_rank)
      poly_map[e] = e;
    else
      poly_map[e] = -(o+1);
  }
  
#ifdef PARALLEL
  // now take turns and set global element numbers
  if (num_procs > 1)
  {
    sendcnt = new int[num_procs];
    recvcnt = new int[num_procs];
    srequest = new MPI_Request[num_procs];
    rrequest = new MPI_Request[num_procs];
    statuses  = new MPI_Status[num_procs];
    // allocate space for message buffers
    bdim = new int[num_procs];
    for (p=0; p < num_procs; p++)
      bdim[p] = 0;
    sbuff = new char*[num_procs];
    rbuff = new char*[num_procs];
    for (p=0; p < num_procs; p++)
    {
      sbuff[p] = 0;
      rbuff[p] = 0;
    }
    pelems = new List*[num_procs];
    for (p=0; p < num_procs; p++)
      pelems[p] = new List();
  }

  MPI_Barrier(MPI_COMM_WORLD);
#endif

  int tet_total, pyr_total, pri_total, hex_total, poly_total;
  int *tri_total, *quad_total, *ngon_total;
  tri_total = new int[nb];
  quad_total = new int[nb];
  ngon_total = new int[nb];

  // Assign global element numbers
  if (my_rank == 0)
  {
    for (b=0; b < nb; b++)
    {
      tri_total[b] = quad_total[b] = ngon_total[b] = 0;
      for (i=0; i < nt[b]; i++)
        if (tri_map[b][i] >= 0)
          tri_map[b][i] = tri_total[b]++;
#ifdef PARALLEL
      for (i=1; i < num_procs; i++)
      {
        MPI_Send(&tri_total[b],1,MPI_INT, i, i, MPI_COMM_WORLD);
        MPI_Recv(&tri_total[b],1,MPI_INT,i,i,MPI_COMM_WORLD, &single_status);
      }
#endif
      for (i=0; i < nq[b]; i++)
        if (quad_map[b][i] >= 0)
          quad_map[b][i] = quad_total[b]++;
#ifdef PARALLEL
      for (i=1; i < num_procs; i++)
      {
        MPI_Send(&quad_total[b],1,MPI_INT, i, i, MPI_COMM_WORLD);
        MPI_Recv(&quad_total[b],1,MPI_INT,i,i,MPI_COMM_WORLD, &single_status);
      }
#endif
      for (i=0; i < ngon[b]; i++)
        if (ngon_map[b][i] >= 0)
          ngon_map[b][i] = ngon_total[b]++;
#ifdef PARALLEL
      for (i=1; i < num_procs; i++)
      {
        MPI_Send(&ngon_total[b],1,MPI_INT, i, i, MPI_COMM_WORLD);
        MPI_Recv(&ngon_total[b],1,MPI_INT,i,i,MPI_COMM_WORLD, &single_status);
      }
#endif
    }

    tet_total = pyr_total = pri_total = hex_total = poly_total = 0;
    for (e=0; e < ntet; e++)
      if (tet_map[e] >= 0)
        tet_map[e] = tet_total++;
#ifdef PARALLEL
    for (i=1; i < num_procs; i++)
    {
      MPI_Send(&tet_total,1,MPI_INT, i, i, MPI_COMM_WORLD);
      MPI_Recv(&tet_total,1,MPI_INT,i,i,MPI_COMM_WORLD, &single_status);
    }
#endif
    //fprintf(out_f,"\nCreate_Hybrid_Maps: Total number of tet elements = %d",tet_total);
    //fflush(out_f);
    for (e=0; e < npyr; e++)
      if (pyr_map[e] >= 0)
        pyr_map[e] = pyr_total++;
#ifdef PARALLEL
    for (i=1; i < num_procs; i++)
    {
      MPI_Send(&pyr_total,1,MPI_INT, i, i, MPI_COMM_WORLD);
      MPI_Recv(&pyr_total,1,MPI_INT,i,i,MPI_COMM_WORLD, &single_status);
    }
#endif
    //fprintf(out_f,"\nCreate_Hybrid_Maps: Total number of pyr elements = %d",pyr_total);
    //fflush(out_f);
    for (e=0; e < npri; e++)
      if (pri_map[e] >= 0)
        pri_map[e] = pri_total++;
#ifdef PARALLEL
    for (i=1; i < num_procs; i++)
    {
      MPI_Send(&pri_total,1,MPI_INT, i, i, MPI_COMM_WORLD);
      MPI_Recv(&pri_total,1,MPI_INT,i,i,MPI_COMM_WORLD, &single_status);
    }
#endif
    //fprintf(out_f,"\nCreate_Hybrid_Maps: Total number of pri elements = %d",pri_total);
    //fflush(out_f);
    for (e=0; e < nhex; e++)
      if (hex_map[e] >= 0)
        hex_map[e] = hex_total++;
#ifdef PARALLEL
    for (i=1; i < num_procs; i++)
    {
      MPI_Send(&hex_total,1,MPI_INT, i, i, MPI_COMM_WORLD);
      MPI_Recv(&hex_total,1,MPI_INT,i,i,MPI_COMM_WORLD, &single_status);
    }
#endif
    //fprintf(out_f,"\nCreate_Hybrid_Maps: Total number of hex elements = %d",hex_total);
    //fflush(out_f);
    for (e=0; e < nply; e++)
      if (poly_map[e] >= 0)
        poly_map[e] = poly_total++;
#ifdef PARALLEL
    for (i=1; i < num_procs; i++)
    {
      MPI_Send(&poly_total,1,MPI_INT, i, i, MPI_COMM_WORLD);
      MPI_Recv(&poly_total,1,MPI_INT,i,i,MPI_COMM_WORLD, &single_status);
    }
#endif
    //fprintf(out_f,"\nCreate_Hybrid_Maps: Total number of poly elements = %d",poly_total);
    //fflush(out_f);
  } else
  {
#ifdef PARALLEL
    for (b=0; b < nb; b++)
    {
      MPI_Recv(&tri_total[b],1,MPI_INT, 0, my_rank, MPI_COMM_WORLD, &single_status);
      for (i=0; i < nt[b]; i++)
        if (tri_map[b][i] >= 0)
          tri_map[b][i] = tri_total[b]++;
      MPI_Send(&tri_total[b],1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
      MPI_Recv(&quad_total[b],1,MPI_INT, 0, my_rank, MPI_COMM_WORLD, &single_status);
      for (i=0; i < nq[b]; i++)
        if (quad_map[b][i] >= 0)
          quad_map[b][i] = quad_total[b]++;
      MPI_Send(&quad_total[b],1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
      MPI_Recv(&ngon_total[b],1,MPI_INT, 0, my_rank, MPI_COMM_WORLD, &single_status);
      for (i=0; i < ngon[b]; i++)
        if (ngon_map[b][i] >= 0)
          ngon_map[b][i] = ngon_total[b]++;
      MPI_Send(&ngon_total[b],1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
    }
    MPI_Recv(&tet_total,1,MPI_INT, 0, my_rank, MPI_COMM_WORLD, &single_status);
    for (e=0; e < ntet; e++)
      if (tet_map[e] >= 0)
        tet_map[e] = tet_total++;
    MPI_Send(&tet_total,1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
    MPI_Recv(&pyr_total,1,MPI_INT, 0, my_rank, MPI_COMM_WORLD, &single_status);
    for (e=0; e < npyr; e++)
      if (pyr_map[e] >= 0)
        pyr_map[e] = pyr_total++;
    MPI_Send(&pyr_total,1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
    MPI_Recv(&pri_total,1,MPI_INT, 0, my_rank, MPI_COMM_WORLD, &single_status);
    for (e=0; e < npri; e++)
      if (pri_map[e] >= 0)
        pri_map[e] = pri_total++;
    MPI_Send(&pri_total,1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
    MPI_Recv(&hex_total,1,MPI_INT, 0, my_rank, MPI_COMM_WORLD, &single_status);
    for (e=0; e < nhex; e++)
      if (hex_map[e] >= 0)
        hex_map[e] = hex_total++;
    MPI_Send(&hex_total,1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
    MPI_Recv(&poly_total,1,MPI_INT, 0, my_rank, MPI_COMM_WORLD, &single_status);
    for (e=0; e < nply; e++)
      if (poly_map[e] >= 0)
        poly_map[e] = poly_total++;
    MPI_Send(&poly_total,1,MPI_INT,0,my_rank,MPI_COMM_WORLD);
#endif
  }

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);

  if (num_procs > 1)
  {
    // communicate and exchange the global number and indices not owned by me
    List **elist;
    int ne, ln;

    elist = new List*[nn];

    for (n=0; n < nn; n++)
      elist[n] = new List();

    int type, num_nodes;
    for (b=0; b < nb; b++)
    {
      MPI_Barrier(MPI_COMM_WORLD);

      for (type=0; type < 3; type++)
      {
        for (n=0; n < nn; n++)
          elist[n]->Redimension(0);
        for (p=0; p < num_procs; p++)
          pelems[p]->Redimension(0);

        switch (type)
        {
          case 0:
            num_nodes = 3;
            for (e=0; e < nt[b]; e++)
            {
              for (i=0; i < num_nodes; i++)
              {
                n = t_n[b][e][i];
                elist[n]->Add_To_List(e);
              }
              if ((p=-tri_map[b][e]-1) >= 0 && p != my_rank)
                pelems[p]->Add_To_List(e);
            }
            break;
          case 1:
            num_nodes = 4;
            for (e=0; e < nq[b]; e++)
            {
              for (i=0; i < num_nodes; i++)
              {
                n = q_n[b][e][i];
                elist[n]->Add_To_List(e);
              }
              if ((p=-quad_map[b][e]-1) >= 0 && p != my_rank)
                pelems[p]->Add_To_List(e);
            }
            break;
          case 4:
            num_nodes = 0;
            for (e=0; e < ngon[b]; e++)
            {
              for (j=1; j <= ngon_n[b][e][0]; j++)
              {
                n = ngon_n[b][e][j];
                elist[n]->Add_To_List(e);
              }
              if ((p=-ngon_map[b][e]-1) >= 0 && p != my_rank)
                pelems[p]->Add_To_List(e);
            }
            break;
        }

        List nodelist;
        // send buffer size
        for (p=0; p < num_procs; p++)
        {
          sendcnt[p] = 0;
          if (p != my_rank)
          {
            sendcnt[p] = sizeof(int);
            for (i=0; i < pelems[p]->max; i++)
            {
              e = pelems[p]->list[i];
              switch (type)
              {
                case 0:
                  num_nodes = 3;
                  break;
                case 1:
                  num_nodes = 4;
                  break;
                case 2:
                  nodelist.Redimension(0);
                  for (j=1; j <= ngon_n[b][e][0]; j++)
                    nodelist.Check_List(ngon_n[b][e][j]);
                  num_nodes = nodelist.max;
                  break;
              }
              sendcnt[p] += (num_nodes*3 + 2)*sizeof(int);
            }
          }
        }
    
        nreq_s = nreq_r = 0;
        for (p=0; p < num_procs; p++)
        {
          recvcnt[p] = 0;
          if (p != my_rank)
          {
            MPI_Isend(&sendcnt[p],1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
            MPI_Irecv(&recvcnt[p],1,MPI_INT,p,      p,MPI_COMM_WORLD,&rrequest[nreq_r]);
            nreq_s++;
            nreq_r++;
          }
        }
        MPI_Waitall(nreq_s,srequest,statuses);
        MPI_Waitall(nreq_r,rrequest,statuses);

        // size message buffers for each processor
        for (p=0; p < num_procs; p++)
        {
          bdim[p] = 0;
          if (p == my_rank || (sendcnt[p] == 0 && recvcnt[p] == 0))
            continue;

          bdim[p] = MAX(sendcnt[p],recvcnt[p]);
          if (sbuff[p] > 0) delete[] sbuff[p];
          if (rbuff[p] > 0) delete[] rbuff[p];
          if (bdim[p] > 0) sbuff[p] = new char[bdim[p]];
          if (bdim[p] > 0) rbuff[p] = new char[bdim[p]];
        }

        // package the list of element per processor
        for (p=0; p < num_procs; p++)
        {
          if (p == my_rank || sendcnt[p] == 0)
            continue;

          sposition=0;
          MPI_Pack(&(pelems[p]->max),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
          for (i=0; i < pelems[p]->max; i++)
          {
            e = pelems[p]->list[i];
            MPI_Pack(&e,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
            switch(type)
            {
              case 0:
                num_nodes = 3;
                MPI_Pack(&num_nodes,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                for (j=0; j < num_nodes; j++)
                {
                  n = t_n[b][e][j];
                  MPI_Pack(&(nmap[n][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(nmap[n][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(nmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                }
                break;
              case 1:
                num_nodes = 4;
                MPI_Pack(&num_nodes,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                for (j=0; j < num_nodes; j++)
                {
                  n = q_n[b][e][j];
                  MPI_Pack(&(nmap[n][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(nmap[n][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(nmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                }
                break;
              case 2:
                nodelist.Redimension(0);
                for (j=1; j <= ngon_n[b][e][0]; j++)
                  nodelist.Check_List(ngon_n[b][e][j]);
                num_nodes = nodelist.max;
                MPI_Pack(&num_nodes,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                for (j=0; j < nodelist.max; j++)
                {
                  n = nodelist.list[j];
                  MPI_Pack(&(nmap[n][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(nmap[n][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                  MPI_Pack(&(nmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                }
                break;
            }
          }
        }

        // exchange messages
        nreq_s = nreq_r = 0;
        for (p=0; p < num_procs; p++)
        {
          if (p == my_rank)
            continue;
          if (sendcnt[p] > 0)
          {
            MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
            nreq_s++;
          }
          if (recvcnt[p] > 0)
          {
            MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,      p,MPI_COMM_WORLD,&rrequest[nreq_r]);
            nreq_r++;
          }
        }
        MPI_Waitall(nreq_s,srequest,statuses);
        MPI_Waitall(nreq_r,rrequest,statuses);

        // unpack element list and re-package the current map for requested elements
        for (p=0; p < num_procs; p++)
        {
          if (p == my_rank)
            continue;

          sposition=rposition=0;
          Pmap *lmap;
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&ne,1,MPI_INT,MPI_COMM_WORLD);
          MPI_Pack(&ne,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
          for (i=0; i < ne; i++)
          {
            MPI_Unpack(rbuff[p],bdim[p],&rposition,&e,1,MPI_INT,MPI_COMM_WORLD);
            MPI_Pack(&e,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
            MPI_Unpack(rbuff[p],bdim[p],&rposition,&ln,1,MPI_INT,MPI_COMM_WORLD);
            lmap = new Pmap[ln];
            m = -1;
            for (n=0; n < ln; n++)
            {
              MPI_Unpack(rbuff[p],bdim[p],&rposition,&lmap[n].global,1,MPI_INT,MPI_COMM_WORLD);
              MPI_Unpack(rbuff[p],bdim[p],&rposition,&lmap[n].proc,1,MPI_INT,MPI_COMM_WORLD);
              MPI_Unpack(rbuff[p],bdim[p],&rposition,&lmap[n].index,1,MPI_INT,MPI_COMM_WORLD);
              // find a node to use with the local node-element hash table
              if (m < 0 && lmap[n].proc == my_rank)
                m = lmap[n].index;
            }
            if (m < 0)
            {
              if (my_rank == 0 && out_f > 0) fprintf(out_f,"\nCreate_Hybrid_Maps: Boundary elements: no local node identified!");
              if (my_rank == 0 && out_f > 0) fflush(out_f);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }
            int le = -1; // local element identity
            for (j=0; j < elist[m]->max && le < 0; j++)
            {
              e = elist[m]->list[j];
              if (type == 2)
              {
                nodelist.Redimension(0);
                for (k=1; k <= ngon_n[b][e][0]; k++)
                  nodelist.Check_List(ngon_n[b][e][k]);
                if (nodelist.max != ln)
                  continue;
              }
              bool ematch = true;
              for (k=0; k < ln && ematch; k++)
              {
                switch(type)
                {
                  case 0: n = t_n[b][e][k]; break;
                  case 1: n = q_n[b][e][k]; break;
                  case 2: n = nodelist.list[k]; break;
                }
                bool nmatch = false;
                for (l=0; l < ln && !nmatch; l++)
                  if (lmap[l].proc == nmap[n][1] && lmap[l].index == nmap[n][2])
                    nmatch = true;
                ematch = nmatch;
              }
              if (ematch)
                le = e;
            }
            if (type == 2) nodelist.Redimension(0);

            if (le < 0)
            {
              if (my_rank == 0 && out_f > 0) fprintf(out_f,"\nCreate_Hybrid_Maps: Boundary Elements: no local element identified!");
              if (my_rank == 0 && out_f > 0) fflush(out_f);
              MPI_Abort(MPI_COMM_WORLD,0);
              exit(0);
            }

            switch(type)
            {
              case 0: MPI_Pack(&tri_map[b][le],1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD); break;
              case 1: MPI_Pack(&quad_map[b][le],1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD); break;
              case 2: MPI_Pack(&ngon_map[b][le],1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD); break;
            }
  
            delete[] lmap;
          }
        }

        // exchange messages
        nreq_s = nreq_r = 0;
        for (p=0; p < num_procs; p++)
        {
          if (p == my_rank)
            continue;
          if (recvcnt[p] > 0)
          {
            MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
            nreq_s++;
          }
          if (sendcnt[p] > 0)
          {
            MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,      p,MPI_COMM_WORLD,&rrequest[nreq_r]);
            nreq_r++;
          }
        }
        MPI_Waitall(nreq_s,srequest,statuses);
        MPI_Waitall(nreq_r,rrequest,statuses);

        // unpack returned messages
        for (p=0; p < num_procs; p++)
        {
          if (p == my_rank)
            continue;

          sposition=rposition=0;
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&ne,1,MPI_INT,MPI_COMM_WORLD);
          for (m=0; m < ne; m++)
          {
            MPI_Unpack(rbuff[p],bdim[p],&rposition,&e,1,MPI_INT,MPI_COMM_WORLD);
            MPI_Unpack(rbuff[p],bdim[p],&rposition,&i,1,MPI_INT,MPI_COMM_WORLD);
            switch(type)
            {
              case 0: tri_map[b][e] = i; break;
              case 1: quad_map[b][e] = i; break;
              case 2: ngon_map[b][e] = i; break;
            }
          }
        }
      }
    }

    for (type=0; type < 5; type++)
    {

      MPI_Barrier(MPI_COMM_WORLD);

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: type = %d",type);
      //fflush(out_f);

      for (n=0; n < nn; n++)
        elist[n]->Redimension(0);
      for (p=0; p < num_procs; p++)
        pelems[p]->Redimension(0);

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: elist and pelems reset");
      //fflush(out_f);

      // create node to element hash table
      // create list of elements to send to each processor
      switch (type)
      {
        case 0:
          num_nodes = 4;
          for (e=0; e < ntet; e++)
          {
            for (i=0; i < num_nodes; i++)
            {
              n = tet_n[e][i];
              elist[n]->Add_To_List(e);
            }
            if ((p=-tet_map[e]-1) >= 0 && p != my_rank)
              pelems[p]->Add_To_List(e);
          }
          break;
        case 1:
          num_nodes = 5;
          for (e=0; e < npyr; e++)
          {
            for (i=0; i < num_nodes; i++)
            {
              n = pyr_n[e][i];
              elist[n]->Add_To_List(e);
            }
            if ((p=-pyr_map[e]-1) >= 0 && p != my_rank)
              pelems[p]->Add_To_List(e);
          }
          break;
        case 2:
          num_nodes = 6;
          for (e=0; e < npri; e++)
          {
            for (i=0; i < num_nodes; i++)
            {
              n = pri_n[e][i];
              elist[n]->Add_To_List(e);
            }
            if ((p=-pri_map[e]-1) >= 0 && p != my_rank)
              pelems[p]->Add_To_List(e);
          }
          break;
        case 3:
          num_nodes = 8;
          for (e=0; e < nhex; e++)
          {
            for (i=0; i < num_nodes; i++)
            {
              n = hex_n[e][i];
              elist[n]->Add_To_List(e);
            }
            if ((p=-hex_map[e]-1) >= 0 && p != my_rank)
              pelems[p]->Add_To_List(e);
          }
          break;
        case 4:
          num_nodes = 0;
          for (e=0; e < nply; e++)
          {
            for (i=1; i <= poly_n[e][0][0]; i++)
              for (j=1; j <= poly_n[e][i][0]; j++)
              {
                n = poly_n[e][i][j];
                elist[n]->Add_To_List(e);
              }
            if ((p=-poly_map[e]-1) >= 0 && p != my_rank)
              pelems[p]->Add_To_List(e);
          }
          break;
      }

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: pelems list made");
      //fflush(out_f);

      List nodelist;
      // send buffer size
      for (p=0; p < num_procs; p++)
      {
        sendcnt[p] = 0;
        if (p != my_rank)
        {
          sendcnt[p] = sizeof(int);
          for (i=0; i < pelems[p]->max; i++)
          {
            e = pelems[p]->list[i];
            switch (type)
            {
              case 0:
                num_nodes = 4;
                break;
              case 1:
                num_nodes = 5;
                break;
              case 2:
                num_nodes = 6;
                break;
              case 3:
                num_nodes = 8;
                break;
              case 4:
                nodelist.Redimension(0);
                for (j=1; j <= poly_n[e][0][0]; j++)
                  for (k=1; k <= poly_n[e][j][0]; k++)
                    nodelist.Check_List(poly_n[e][j][k]);
                num_nodes = nodelist.max;
                break;
            }
            sendcnt[p] += (num_nodes*3 + 2)*sizeof(int);
          }
        }
      }
    
      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: sendcnt made");
      //fflush(out_f);

      nreq_s = nreq_r = 0;
      for (p=0; p < num_procs; p++)
      {
        recvcnt[p] = 0;
        if (p != my_rank)
        {
          MPI_Isend(&sendcnt[p],1,MPI_INT,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
          MPI_Irecv(&recvcnt[p],1,MPI_INT,p,      p,MPI_COMM_WORLD,&rrequest[nreq_r]);
          nreq_s++;
          nreq_r++;
        }
      }
      MPI_Waitall(nreq_s,srequest,statuses);
      MPI_Waitall(nreq_r,rrequest,statuses);

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: sendcnt exchanged");
      //fflush(out_f);

      // size message buffers for each processor
      for (p=0; p < num_procs; p++)
      {
        bdim[p] = 0;
        if (p == my_rank || (sendcnt[p] == 0 && recvcnt[p] == 0))
          continue;

        bdim[p] = MAX(sendcnt[p],recvcnt[p]);
        if (sbuff[p] > 0) delete[] sbuff[p];
        if (rbuff[p] > 0) delete[] rbuff[p];
        if (bdim[p] > 0) sbuff[p] = new char[bdim[p]];
        if (bdim[p] > 0) rbuff[p] = new char[bdim[p]];
      }

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: buffers resized");
      //fflush(out_f);

      // package the list of element per processor
      for (p=0; p < num_procs; p++)
      {
        if (p == my_rank || sendcnt[p] == 0)
          continue;

        sposition=0;
        MPI_Pack(&(pelems[p]->max),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        for (i=0; i < pelems[p]->max; i++)
        {
          e = pelems[p]->list[i];
          MPI_Pack(&e,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
          switch(type)
          {
            case 0:
              num_nodes = 4;
              MPI_Pack(&num_nodes,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              for (j=0; j < num_nodes; j++)
              {
                n = tet_n[e][j];
                MPI_Pack(&(nmap[n][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              }
              break;
            case 1:
              num_nodes = 5;
              MPI_Pack(&num_nodes,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              for (j=0; j < num_nodes; j++)
              {
                n = pyr_n[e][j];
                MPI_Pack(&(nmap[n][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              }
              break;
            case 2:
              num_nodes = 6;
              MPI_Pack(&num_nodes,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              for (j=0; j < num_nodes; j++)
              {
                n = pri_n[e][j];
                MPI_Pack(&(nmap[n][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              }
              break;
            case 3:
              num_nodes = 8;
              MPI_Pack(&num_nodes,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              for (j=0; j < num_nodes; j++)
              {
                n = hex_n[e][j];
                MPI_Pack(&(nmap[n][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              }
              break;
            case 4:
              nodelist.Redimension(0);
              for (j=1; j <= poly_n[e][0][0]; j++)
                for (k=1; k <= poly_n[e][j][0]; k++)
                  nodelist.Check_List(poly_n[e][j][k]);
              num_nodes = nodelist.max;
              MPI_Pack(&num_nodes,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              for (j=0; j < nodelist.max; j++)
              {
                n = nodelist.list[j];
                MPI_Pack(&(nmap[n][0]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][1]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
                MPI_Pack(&(nmap[n][2]),1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
              }
              break;
          }
        }
      }

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: messages packed");
      //fflush(out_f);

      // exchange messages
      nreq_s = nreq_r = 0;
      for (p=0; p < num_procs; p++)
      {
        if (p == my_rank)
          continue;
        if (sendcnt[p] > 0)
        {
          MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
          nreq_s++;
        }
        if (recvcnt[p] > 0)
        {
          MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,      p,MPI_COMM_WORLD,&rrequest[nreq_r]);
          nreq_r++;
        }
      }
      MPI_Waitall(nreq_s,srequest,statuses);
      MPI_Waitall(nreq_r,rrequest,statuses);

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: messages exhanged");
      //fflush(out_f);

      // unpack element list and re-package the current map for requested elements
      for (p=0; p < num_procs; p++)
      {
        if (p == my_rank)
          continue;

        sposition=rposition=0;
        Pmap *lmap;
        MPI_Unpack(rbuff[p],bdim[p],&rposition,&ne,1,MPI_INT,MPI_COMM_WORLD);
        MPI_Pack(&ne,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
        for (i=0; i < ne; i++)
        {
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&e,1,MPI_INT,MPI_COMM_WORLD);
          MPI_Pack(&e,1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD);
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&ln,1,MPI_INT,MPI_COMM_WORLD);
          lmap = new Pmap[ln];
          m = -1;
          for (n=0; n < ln; n++)
          {
            MPI_Unpack(rbuff[p],bdim[p],&rposition,&lmap[n].global,1,MPI_INT,MPI_COMM_WORLD);
            MPI_Unpack(rbuff[p],bdim[p],&rposition,&lmap[n].proc,1,MPI_INT,MPI_COMM_WORLD);
            MPI_Unpack(rbuff[p],bdim[p],&rposition,&lmap[n].index,1,MPI_INT,MPI_COMM_WORLD);
            // find a node to use with the local node-element hash table
            if (m < 0 && lmap[n].proc == my_rank)
              m = lmap[n].index;
          }
          if (m < 0)
          {
            if (my_rank == 0 && out_f > 0) fprintf(out_f,"\nCreate_Hybrid_Maps: no local node identified!");
            if (my_rank == 0 && out_f > 0) fflush(out_f);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          }

          int le = -1; // local element identity
          for (j=0; j < elist[m]->max && le < 0; j++)
          {
            e = elist[m]->list[j];
            if (type == 4)
            {
              nodelist.Redimension(0);
              for (j=1; j <= poly_n[e][0][0]; j++)
                for (k=1; k <= poly_n[e][j][0]; k++)
                  nodelist.Check_List(poly_n[e][j][k]);
              if (nodelist.max != ln)
                continue;
            }
            bool ematch = true;
            for (k=0; k < ln && ematch; k++)
            {
              switch(type)
              {
                case 0: n = tet_n[e][k]; break;
                case 1: n = pyr_n[e][k]; break;
                case 2: n = pri_n[e][k]; break;
                case 3: n = hex_n[e][k]; break;
                case 4: n = nodelist.list[k]; break;
              }
              bool nmatch = false;
              for (l=0; l < ln && !nmatch; l++)
                if (lmap[l].proc == nmap[n][1] && lmap[l].index == nmap[n][2])
                  nmatch = true;
              ematch = nmatch;
            }
            if (ematch)
              le = e;
          }
          if (type == 4) nodelist.Redimension(0);

          if (le < 0)
          {
            if (my_rank == 0 && out_f > 0) fprintf(out_f,"\nCreate_Hybrid_Maps: no local element identified!");
            if (my_rank == 0 && out_f > 0) fflush(out_f);
            MPI_Abort(MPI_COMM_WORLD,0);
            exit(0);
          }

          switch(type)
          {
            case 0: MPI_Pack(&tet_map[le],1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD); break;
            case 1: MPI_Pack(&pyr_map[le],1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD); break;
            case 2: MPI_Pack(&pri_map[le],1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD); break;
            case 3: MPI_Pack(&hex_map[le],1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD); break;
            case 4: MPI_Pack(&poly_map[le],1,MPI_INT,sbuff[p],bdim[p],&sposition,MPI_COMM_WORLD); break;
          }

          delete[] lmap;
        }
      }

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: messages repacked");
      //fflush(out_f);

      // exchange messages
      nreq_s = nreq_r = 0;
      for (p=0; p < num_procs; p++)
      {
        if (p == my_rank)
          continue;
        if (recvcnt[p] > 0)
        {
          MPI_Isend(sbuff[p],bdim[p],MPI_CHAR,p,my_rank,MPI_COMM_WORLD,&srequest[nreq_s]);
          nreq_s++;
        }
        if (sendcnt[p] > 0)
        {
          MPI_Irecv(rbuff[p],bdim[p],MPI_CHAR,p,      p,MPI_COMM_WORLD,&rrequest[nreq_r]);
          nreq_r++;
        }
      }
      MPI_Waitall(nreq_s,srequest,statuses);
      MPI_Waitall(nreq_r,rrequest,statuses);

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: messages re-exchanged");
      //fflush(out_f);

      // unpack returned messages
      for (p=0; p < num_procs; p++)
      {
        if (p == my_rank)
          continue;

        sposition=rposition=0;
        MPI_Unpack(rbuff[p],bdim[p],&rposition,&ne,1,MPI_INT,MPI_COMM_WORLD);
        for (m=0; m < ne; m++)
        {
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&e,1,MPI_INT,MPI_COMM_WORLD);
          MPI_Unpack(rbuff[p],bdim[p],&rposition,&i,1,MPI_INT,MPI_COMM_WORLD);
          switch(type)
          {
            case 0: tet_map[e] = i; break;
            case 1: pyr_map[e] = i; break;
            case 2: pri_map[e] = i; break;
            case 3: hex_map[e] = i; break;
            case 4: poly_map[e] = i; break;
          }
        }
      }

      //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: messages unpacked and saved");
      //fflush(out_f);

    }

    // release communication memory
    for (p=0; p < num_procs; p++)
    {
      if (sbuff[p] > 0) delete[] sbuff[p];
      if (rbuff[p] > 0) delete[] rbuff[p];
      delete pelems[p];
    }
    delete[] pelems;
    delete[] sendcnt;
    delete[] recvcnt;
    delete[] srequest;
    delete[] rrequest;
    delete[] statuses;
    // allocate space for message buffers
    delete[] bdim;
    delete[] sbuff;
    delete[] rbuff;
    for (n=0; n < nn; n++)
      delete elist[n];
    delete[] elist;
  }
#endif

  //if (my_rank == 0) fprintf(out_f,"\nCREATE_HYBRID_MAPS: done");
  //fflush(out_f);

  return;
}
