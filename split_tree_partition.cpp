#include <stdio.h>
#include <sys/param.h>
#include "Point.h"
#include "List.h"
#include "Util.h"
#include "split_tree_partition.h"

class split_obj
{
  public:
  split_obj() { vlist = 0; }

  Point lo, hi;
  List *vlist;
};

int split_tree_partition(int nv, Point vert[], int nparts, int part[], int mode, FILE *outf)
{
  int dir, i, j, k, m, n, nsplit;
  double dx, dy, dz, v;
  split_obj *split;
  int flag = 1;

  if (nparts > nv)
  {
    fprintf(stderr,"\nSPLIT_TREE_PARTITION: Number of partitions exceeds number of vertices!\n");
    fflush(stderr);
    return(0);
  } else if (nparts == nv)
  {
    for (n=0; n < nv; n++)
      part[n] = n;
    return(flag);
  }

  // create all split objects
  split = new split_obj[nparts];
  for (i=0; i < nparts; i++)
    split[i].vlist = new List();

  // initialize first object
  split[0].lo = vert[0];
  split[0].hi = vert[0];
  split[0].vlist->Redimension(nv);
  split[0].vlist->max = 0;
  for (n=0; n < nv; n++)
  {
    split[0].vlist->Add_To_List(n);
    for (i=0; i < 3; i++)
    {
      split[0].lo[i] = MIN(split[0].lo[i],vert[n][i]);
      split[0].hi[i] = MAX(split[0].hi[i],vert[n][i]);
    }
  }
  nsplit = 1;

  // perform splits on largest object until # of splits reached
  while (nsplit < nparts)
  {
    if (outf) fprintf(outf,"\nCurrent number of partitions = %d",nsplit);
    if (outf) fflush(outf);

    // find split with largest number of vertices
    j = 0;
    m = 0;
    for (i=0; i < nsplit; i++)
    {
      if (split[i].vlist->max > m)
      {
        m = split[i].vlist->max;
        j = i;
      }
    }

    // determine split direction
    dx = split[j].hi[0]-split[j].lo[0];
    dy = split[j].hi[1]-split[j].lo[1];
    dz = split[j].hi[2]-split[j].lo[2];
    if (dx >= dy && dx >= dz)
      dir = 0;
    else if (dy >= dx && dy >= dz)
      dir = 1;
    else
      dir = 2;

    int *sn = new int[split[j].vlist->max];
    double *sv = new double[split[j].vlist->max];
    List tmp;
    tmp.Redimension(split[j].vlist->max);
    tmp.max = 0;
    for (i=0; i < split[j].vlist->max; i++)
    {
      n=split[j].vlist->list[i];
      v=vert[n][dir];
      sn[i] = i;
      sv[i] = v;
      tmp.Add_To_List(n);
    }
    int lv = split[j].vlist->max;
    sort_lt(lv,sn,sv);

    double mid = (split[j].hi[dir] + split[j].lo[dir])*0.5;
    if (mode == 0)
    {
      k=lv/2;
      v=sv[sn[k]];
    } else
    {
      k=0;
      v=sv[sn[k]];
      while (v < mid && k < lv-2)
        v = sv[sn[++k]];
    }
    if (outf && dir == 0) fprintf(outf,"\nSplitting partition %d at X = %lg",j,v);
    if (outf && dir == 1) fprintf(outf,"\nSplitting partition %d at Y = %lg",j,v);
    if (outf && dir == 2) fprintf(outf,"\nSplitting partition %d at Z = %lg",j,v);
    if (outf) fflush(outf);
    
    // set new split highs and lows
    split[nsplit].lo = split[j].lo;
    split[nsplit].hi = split[j].hi;
    split[j].hi[dir] = v;
    split[nsplit].lo[dir] = v;

    // store vertex lists
    split[nsplit].vlist->Redimension(split[j].vlist->max);
    split[nsplit].vlist->max = 0;
    split[j].vlist->max = 0;
    for (i=0; i < k; i++)
    {
      n=tmp.list[sn[i]];
      split[j].vlist->Add_To_List(n);
    }
    for (i=k; i < lv; i++)
    {
      n=tmp.list[sn[i]];
      split[nsplit].vlist->Add_To_List(n);
    }
    split[j].vlist->Redimension(split[j].vlist->max);
    split[nsplit].vlist->Redimension(split[nsplit].vlist->max);
    delete[] sn;
    delete[] sv;
    tmp.Redimension(0);

    nsplit++;
  }
  if (outf)
  {
    for (i=0; i < nparts; i++)
      fprintf(outf,"\nPartition %d has %d vertices",i,split[i].vlist->max);
    fflush(outf);
  }

  // record partition for each vertex
  for (n=0; n < nv; n++)
    part[n] = 0;

  for (i=0; i < nparts; i++)
  {
    for (j=0; j < split[i].vlist->max; j++)
    {
      n = split[i].vlist->list[j];
      part[n] = i;
    }
  }
  
  // free up splits
  for (i=0; i < nparts; i++)
    delete split[i].vlist;
  delete[] split;

  return(flag);

}
