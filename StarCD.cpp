#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "Point.h"
#include "convert.h"
#include "Util.h"

int GetNumberLines(const char *file)
{
  FILE *fp = fopen(file, "r"); 
  int ch, lines = 0, prev = '\n' /* so empty files have no lines */;

  if(fp)
  {
    while((ch = fgetc(fp)) != EOF) /* Read all chars in the file. */
    {
      if(ch == '\n')
        ++lines; /* Bump the counter for every newline. */

      prev = ch; /* Keep a copy to later test whether... */
    }

    fclose(fp);

    if(prev != '\n') /* ...the last line did not end in a newline. */
      ++lines; /* If so, add one more to the total. */

    //printf("lines = %d\n", lines);

    return lines;
  }
  else
  {
    printf("ERROR: (%s) failed to open <%s> (file <%s>, line %d)\n", __func__, file, __FILE__, __LINE__);
    return -1;
  }
}

int StarCD_read(char *fname, int &nnodes, Point **nodes, int &nb, char ***b_name,
                int **ntri, int ****tri_conn,
                int **nquad, int ****quad_conn,
                int **ngon, int ****gon_conn,
                int &nvol, char ***vol_name, 
                int &ntet, int ***tet_conn, int **tet_vc,
                int &npyr, int ***pyr_conn, int **pyr_vc,
                int &npri, int ***pri_conn, int **pri_vc,
                int &nhex, int ***hex_conn, int **hex_vc,
                int &nply, int ****poly_conn, int **poly_vc)
{
  FILE *fp;
  int i, j;
  int ncells, nbface;
  int temp;

  char buff[50];
  char basepath[50];
  char *pnt = NULL;
  char inp_file[50];
  char cel_file[50];
  char vrt_file[50];
  char bnd_file[50];

  memset(basepath, '\0', 50);

  if((pnt = strrchr(fname, '/'))) 
    strncpy(basepath, fname, pnt-fname + 1);

  printf("\n\nReading STAR-CD file <%s> (basepath is <%s>)\n\n", fname, basepath);

  //-----------------------Reading StarCD files--------------------------//

  strcpy(inp_file, fname);
  strcpy(cel_file, basepath);
  strcpy(vrt_file, basepath);
  strcpy(bnd_file, basepath);

  printf("Reading input file <%s>\n\n", inp_file);

  int nlines = GetNumberLines(inp_file); 
  //fscanf(fp,"%d",&nnodes);

  if ((fp=fopen(inp_file, "r")) == NULL)
  {
    printf("Error opening input file <%s>\n", inp_file); 
    fflush(stdout);
    return -1;
  }

  nb = 0;
  nvol = 0;

  char strbuff[1000];
  char *ptr;
  int *tag = (int*)calloc(nlines, sizeof(int));
  char **vol_names = (char **) calloc(nlines, sizeof(char*));
  for(i = 0; i < nlines; i++) 
    vol_names[i] = (char *) calloc(50, sizeof(char));

  for(i = 0; i < nlines; i++)
  {
    
    fgets(strbuff, 1000, fp);

    if(strncmp(strbuff, "CTAB", 4) == 0) 
    {
      sscanf(strbuff, "%s %d %s  %d   %d   %d   %d   %d   %d %s", buff, &tag[nvol], buff, &temp, &temp, &temp, &temp, &temp, &temp, vol_names[nvol]);   
      //printf("nvol = %d, volume name is <%s> (tag = %d)\n", nvol + 1, vol_names[nvol], tag[nvol]);
      nvol++;
    }
    else if(strncmp(strbuff, "RDEF", 4) == 0) 
    {
      nb++;
      //ptr = strtok(strbuff, ",");

      //while(ptr != NULL)
      //{
      //  printf("%s\n", ptr);
      //  ptr = strtok(NULL, ",");
      //}
    }
    else if(strncmp(strbuff, "cread", 5) == 0) 
    {
      ptr = strtok(strbuff, ",");
      ptr = strtok(NULL, ",");
      strcat(cel_file, ptr);
      printf("\tCell file is     <%s>\n", cel_file);
    }
    else if(strncmp(strbuff, "vread", 5) == 0) 
    {
      ptr = strtok(strbuff, ",");
      ptr = strtok(NULL, ",");
      strcat(vrt_file, ptr);
      printf("\tVertex file is   <%s>\n", vrt_file);
    }
    else if(strncmp(strbuff, "bread", 5) == 0) 
    {
      ptr = strtok(strbuff, ",");
      ptr = strtok(NULL, ",");
      strcat(bnd_file, ptr);
      printf("\tBoundary file is <%s>\n", bnd_file);
    }
    else
    {
      printf("WARNING: %s read a line it doesn't know what to do with (file = <%s>, line = %d)\n", __func__, __FILE__, __LINE__);
    }
  }
  printf("\n");
  
  int unique_nvol = nvol;

  for(i = 0; i < nvol; i++) 
  {
    for(j = 0; j < i; j++) 
    {
      if(tag[j] == tag[i]) 
      {
        unique_nvol--;
        tag[j] *= -1;
      }
    }
  }

  printf("Number of unique volume conditions is %d (removed %d duplicates)\n", unique_nvol, nvol - unique_nvol);
  printf("Number of boundary condtions is %d\n", nb);

  (*vol_name) = (char **) calloc(unique_nvol, sizeof(char*));
  for(i = 0; i < unique_nvol; i++) 
    (*vol_name)[i] = (char *) calloc(50, sizeof(char));

  temp = 0;

  for(i = 0; i < nvol; i++) 
    if(tag[i] >= 0)
      memcpy((*vol_name)[temp++], vol_names[i], 50*sizeof(char));
  
  printf("\nVolume condtions:\n");
  for(i = 0; i < unique_nvol; i++) 
    printf("\tVolume %d: name = <%s>\n", i + 1, (*vol_name)[i]);

  printf("\n");

  nvol = unique_nvol;

  fclose(fp);

  free(tag);
  for(i = 0; i < nlines; i++) 
    free(vol_names[i]);
  free(vol_names);

  // Reading Coordinates

  printf("Reading coordinates <%s>\n", vrt_file);

  nnodes = GetNumberLines(vrt_file); 
  //fscanf(fp,"%d",&nnodes);

  if ((fp=fopen(vrt_file, "r")) == NULL)
  {
    printf("Error opening coordinates file <%s>\n", vrt_file); 
    fflush(stdout);
    return -1;
  }

  printf("\n\tNumber of nodes = %d \n\n", nnodes);
  fflush(stdout);

  (*nodes) = (Point *) calloc(nnodes, sizeof(Point));

  for(i = 0; i < nnodes; i++)
    fscanf(fp,"%d   %lf   %lf   %lf ",&temp, &(*nodes)[i][0], &(*nodes)[i][1], &(*nodes)[i][2]);   

  fclose(fp);

  //
  // Reading Cell to Node Connectivity:
  // For hex cells all 8 nodes are going to be different
  // For prism cells node numbers (3 and 4) and (7 and 8) are going to be the same 
  // For pyramid cells node numbers (5 and 6) are going to be the same
  // For tet cells node numbers (3 and 4) and (5, 6, 7 and 8) are going to be the same
  //

  printf("Reading cell connectivity <%s>\n", cel_file);

  //fscanf(fp,"%d",&ncells);
  ncells = GetNumberLines(cel_file);

  if ((fp = fopen(cel_file, "r")) == NULL)
  {
    printf("Error opening cell connectivity file <%s>\n", cel_file);
    fflush(stdout);
    exit(EXIT_FAILURE);
  }

  printf("\n\tNumber of cells = %d \n\n", ncells);
  fflush(stdout);

  int **c2n = (int**) calloc(ncells, sizeof(int *));
  for(i = 0; i < ncells; i++)
    c2n[i] = (int *)calloc(8, sizeof(int));

  int *vol_tag = (int *)calloc(ncells, sizeof(int)); //volume tag
  int* cell_size = (int *)calloc(ncells, sizeof(int));

  ntet = 0;
  npyr = 0;
  npri = 0;
  nhex = 0;
  nply = 0;

  for(i = 0; i < ncells; i++)
  {
    fscanf(fp, "%d" ,&temp); //cell number

    for(j = 0; j < 8; j++)
      fscanf(fp, "%d", &c2n[i][j]); //connectivity

    fscanf(fp, "%d %d", &vol_tag[i], &temp); //volume tag and something else

    if ((c2n[i][2] == c2n[i][3]) && (c2n[i][4] == c2n[i][5])) 
    {
      ntet++;     // no of tet cells
      c2n[i][3] = c2n[i][4];
      c2n[i][4] = c2n[i][5];

      cell_size[i] = 4;
    }
    else if ((c2n[i][2] != c2n[i][3]) && (c2n[i][4] == c2n[i][5]))
    {
      npyr++;     // no of pyramid cells
      cell_size[i] = 5;
    }
    else if ((c2n[i][2] == c2n[i][3]) && (c2n[i][4] != c2n[i][5]) && (c2n[i][6] == c2n[i][7])) 
    {
      npri++;       // no of prism cells
      c2n[i][3] = c2n[i][4];
      c2n[i][4] = c2n[i][5];
      c2n[i][5] = c2n[i][6];

      cell_size[i] = 6;
    }
    else  
    {
      nhex++;        // no of hex cells
      cell_size[i] = 8;
    }

    //
    //c index cell nodes
    //
    for(j = 0; j < 8; j++)
      c2n[i][j]--;
  }

  printf("\tNumber of tetrahedral cells = %d\n", ntet);
  printf("\tNumber of pyramid cells     = %d\n", npyr);
  printf("\tNumber of prism cells       = %d\n", npri);
  printf("\tNumber of hexahedral cells  = %d\n\n", nhex);

  fclose(fp);

  // Reading Boundary Information
  // For quad faces all 4 nodes are going to be different
  // For triangular faces node numbers (3 and 4) are going to be same

  printf("Reading boundary information <%s>\n", bnd_file);

  nbface = GetNumberLines(bnd_file);
  //fscanf(fp, "%d", &nbface);

  if ((fp = fopen(bnd_file, "r")) == NULL)
  {
    printf("Error opening the boundary file <%s>\n", bnd_file);
    fflush(stdout); 
    return -1;
  }

  printf("\n\tNumber of boundary faces = %d\n", nbface);
  fflush(stdout);

  int **ibface = (int**) calloc(nbface, sizeof(int *));
  char **ibname = (char **) malloc(nbface*sizeof(char *));

  for(i = 0; i < nbface; i++)
  {
    ibface[i] = (int *) calloc(5, sizeof(int));
    ibname[i] = (char *) malloc(50*sizeof(char));
  }

  int nbface_quad = 0;
  int nbface_tri  = 0;

  for(i = 0; i < nbface; i++)
  {
    fscanf(fp, "%d  %d  %d  %d  %d  %d  %d  %s", &temp, &ibface[i][0],
           &ibface[i][1], &ibface[i][2], &ibface[i][3], &ibface[i][4],
           &temp, ibname[i]);
    //printf("nodes = (%d, %d, %d, %d)\n", ibface[i][0], ibface[i][1], ibface[i][2], ibface[i][3]);
    //printf("tag = %d, boundary name = %s\n", ibface[i][4], ibname[i]);

    if (ibface[i][2] == ibface[i][3])
      nbface_tri++;   // no of triangular faces
    else
      nbface_quad++;  // no of quad faces

    //
    //c index boundary nodes and tags
    //
    for(j = 0; j < 5; j++)
      ibface[i][j]--;

  }

  (*b_name) = (char **) malloc(nb*sizeof(char *));

  for(i = 0; i < nb; i++) 
    (*b_name)[i] = (char *) malloc(50*sizeof(char));

  fclose(fp);

  //
  //Convert data structures for to be compatible with structures expected by
  //convert
  //

  (*ntri) = (int*) calloc(nb, sizeof(int));
  (*nquad) = (int*) calloc(nb, sizeof(int));
  (*ngon) = (int*) calloc(nb, sizeof(int));

  for(i = 0; i < nb; i++) 
  {
    (*ntri)[i] = 0;
    (*nquad)[i] = 0;
    (*ngon)[i] = 0;
  }

  for(i = 0; i < nbface; i++) 
  {
    if (ibface[i][2] == ibface[i][3])
      (*ntri)[ibface[i][4]]++;   // no of triangular faces
    else
      (*nquad)[ibface[i][4]]++;   // no of quadrilateral faces

    strcpy((*b_name)[ibface[i][4]], ibname[i]);
  }

  (*tri_conn) = (int***) calloc(nb, sizeof(int**));
  (*quad_conn) = (int***) calloc(nb, sizeof(int**));
  (*gon_conn) = (int***) calloc(nb, sizeof(int**));

  for(i = 0; i < nb; i++) 
  {
    if((*ntri)[i] != 0) 
      (*tri_conn)[i] = (int**)calloc((*ntri)[i], sizeof(int*));
    else
      (*tri_conn)[i] =  NULL;

    if((*nquad)[i] != 0) 
      (*quad_conn)[i] = (int**)calloc((*nquad)[i], sizeof(int*));
    else
      (*quad_conn)[i] = NULL;

    (*gon_conn)[i] = NULL;
  }

  for(i = 0; i < nb; i++) 
  {
    for(j = 0; j < (*ntri)[i]; j++) 
      (*tri_conn)[i][j] = (int*) calloc(3, sizeof(int));

    for(j = 0; j < (*nquad)[i]; j++) 
      (*quad_conn)[i][j] = (int*) calloc(4, sizeof(int));
  }
  
  int face, tri_cnt[nb], quad_cnt[nb];
  for(i = 0; i < nb; i++) 
  {
    tri_cnt[i] = 0;
    quad_cnt[i] = 0;
  }

  int itri, iquad;
  for(i = 0; i < nbface; i++) 
  {
    face = ibface[i][4]; 
    if (ibface[i][2] == ibface[i][3])
    {
      itri = tri_cnt[face]++;
      (*tri_conn)[face][itri][0] = ibface[i][0];
      (*tri_conn)[face][itri][1] = ibface[i][1];
      (*tri_conn)[face][itri][2] = ibface[i][2];
    }
    else
    {
      iquad = quad_cnt[face]++;
      (*quad_conn)[face][iquad][0] = ibface[i][0];
      (*quad_conn)[face][iquad][1] = ibface[i][1];
      (*quad_conn)[face][iquad][2] = ibface[i][2];
      (*quad_conn)[face][iquad][3] = ibface[i][3];
    }
  }

  itri = iquad = 0;
  for(i = 0; i < nb; i++) 
  {
    itri += tri_cnt[i];
    iquad += quad_cnt[i];
  }
  assert(itri == nbface_tri && iquad == nbface_quad);

  (*tet_conn) = (int **)calloc(ntet, sizeof(int*));
  (*tet_vc) = (int *)calloc(ntet, sizeof(int));
  (*pyr_conn) = (int **)calloc(npyr, sizeof(int*));
  (*pyr_vc) = (int *)calloc(npyr, sizeof(int));
  (*pri_conn) = (int **)calloc(npri, sizeof(int*));
  (*pri_vc) = (int *)calloc(npri, sizeof(int));
  (*hex_conn) = (int **)calloc(nhex, sizeof(int*));
  (*hex_vc) = (int *)calloc(nhex, sizeof(int));
  (*poly_conn) = NULL;
  (*poly_vc) = NULL;

  for(i = 0; i < ntet; i++) 
    (*tet_conn)[i] = (int *)calloc(4, sizeof(int));
  for(i = 0; i < npyr; i++) 
    (*pyr_conn)[i] = (int *)calloc(5, sizeof(int));
  for(i = 0; i < npri; i++) 
    (*pri_conn)[i] = (int *)calloc(6, sizeof(int));
  for(i = 0; i < nhex; i++) 
    (*hex_conn)[i] = (int *)calloc(8, sizeof(int));

  int tet_cnt = 0, pyramid_cnt = 0, prism_cnt = 0, hex_cnt = 0;

  for(i = 0; i < ncells; i++) 
  {
    if (cell_size[i] == 4) 
    {
      for(j = 0; j < 4; j++) 
        (*tet_conn)[tet_cnt][j] = c2n[i][j];

      (*tet_vc)[tet_cnt] = vol_tag[i];

      tet_cnt++;
    }
    else if (cell_size[i] == 5)
    {
      for(j = 0; j < 5; j++) 
        (*pyr_conn)[pyramid_cnt][j] = c2n[i][j];

      (*pyr_vc)[pyramid_cnt] = vol_tag[i];

      pyramid_cnt++;
    }
    else if (cell_size[i] == 6) 
    {
      for(j = 0; j < 6; j++) 
        (*pri_conn)[prism_cnt][j] = c2n[i][j];

      (*pri_vc)[prism_cnt] = vol_tag[i];

      prism_cnt++;
    }
    else  
    {
      for(j = 0; j < 8; j++) 
        (*hex_conn)[hex_cnt][j] = c2n[i][j];

      (*hex_vc)[hex_cnt] = vol_tag[i];

      hex_cnt++;
    }
  }

  assert(tet_cnt == ntet && pyramid_cnt == npyr && 
         prism_cnt == npri && hex_cnt == nhex);
  
  for(i = 0; i < ncells; i++) 
    free(c2n[i]);
  free(c2n);

  free(vol_tag);

  free(cell_size);

  for(i = 0; i < nbface; i++) 
  {
    free(ibface[i]);
    free(ibname[i]);
  }
  free(ibface);
  free(ibname);

  printf("\nReading <%s> complete\n", fname);

  return 0;
}

int StarCD_write(char *fname, const int nnodes, Point *nodes,
                 const int nb, char **b_name,
                 const int * const ntri, int ***tri_conn,
                 const int * const nquad, int ***quad_conn,
                 const int * const ngon, int ***gon_conn,
                 const int nvol, char **vol_name, 
                 const int ntet, int **tet_conn, const int * const tet_vc,
                 const int npyr, int **pyr_conn, const int * const pyr_vc,
                 const int npri, int **pri_conn, const int * const pri_vc,
                 const int nhex, int **hex_conn, const int * const hex_vc,
                 const int nply, int ***poly_conn, const int * const poly_vc)
{
  int i, j, n;
  char basename[50];
  char inp_file[50];
  char cel_file[50];
  char vrt_file[50];
  char bnd_file[50];
  char *pnt = NULL;
  FILE *fp = NULL;

  memset(basename, '\0', 50);

  if((pnt = strrchr(fname, '.'))) 
    strncpy(basename, fname, pnt-fname);

  printf("\n\nWriting Star-CD file <%s> (basename is <%s>)\n", fname, basename);

  sprintf(inp_file, "%s.inp", basename);
  sprintf(cel_file, "%s.cel", basename);
  sprintf(vrt_file, "%s.vrt", basename);
  sprintf(bnd_file, "%s.bnd", basename);

  //
  //write mesh info to .inp file
  //
  printf("\nWriting input file to <%s>\n", inp_file);

  if((fp = fopen(inp_file, "w")) == NULL)
  {
    printf("Error opening input file <%s>\n", inp_file);
    return -1;
  }

  printf("\n\tNumber of volume tags = %d\n", nvol);
  //
  //write volume condition data, I have no idea what the six integers are
  //supposed to be
  //
  for(i = 0; i < nvol; i++) 
  {
    fprintf(fp, "CTAB %8d %s %8d %8d %8d %8d %8d %8d %s\n", i + 1, "VOLUME_TYPE", 2, 0, 1, 1, 0, 0, vol_name[i]);
  }

  printf("\n\tNumber of boundaries = %d\n", nb);
  //
  //write boundary names
  //
  for(i = 0; i < nb; i++) 
    fprintf(fp, "RDEF,%d,%s\n", i + 1, b_name[i]);

  //
  //information about .cel file
  //
  fprintf(fp, "cread,%s,%s,,,,%s\n", cel_file, "icvo", "coded");

  //
  //information about .vrt file
  //
  fprintf(fp, "vread,%s,%s,,,,%s\n", vrt_file, "icvo", "coded");

  //
  //information about .bnd file
  //
  fprintf(fp, "bread,%s,%d,%s,,%s,%s,%d,%d\n", bnd_file, 0, "all", "add", "code", 0, 0);

  fclose(fp);

  //
  //write node coordinates to .vrt file
  //
  printf("\nWriting coordinates to <%s>\n", vrt_file);

  printf("\n\tNumber of nodes = %d\n", nnodes);

  if((fp = fopen(vrt_file, "w")) == NULL)
  {
    printf("Error opening coordinates file <%s>\n", vrt_file);
    return -1;
  }

  for(n = 0; n < nnodes; n++) 
  {
    //printf("%lf %lf %lf\n", coords[n][0], coords[n][1], coords[n][2]);
    fprintf(fp, "%-6d   %-18.10e  %-18.10e  %-18.10e\n", n + 1, nodes[n][0], nodes[n][1], nodes[n][2]);
  }

  fclose(fp);

  //
  //write cell to node connectivity to .cel file
  //
  printf("\nWriting elements to <%s>\n", cel_file);
  printf("\n\tNumber of cells = %d\n\n", ntet + npyr + npri + nhex);

  if((fp = fopen(cel_file, "w")) == NULL)
  {
    printf("Error opening cell to node connectivity file <%s>\n", cel_file);
    return -1;
  }

  printf("\tNumber of tetrahedra = %d\n", ntet);
  printf("\tNumber of pyramids   = %d\n", npyr);
  printf("\tNumber of prisms     = %d\n", npri);
  printf("\tNumber of hexahedra  = %d\n", nhex);

  char cel_format_string[] = "%9d       %8d %8d %8d %8d %8d %8d %8d %8d %8d %4d\n"; 

  for(n = 0; n < ntet; n++) 
  {
    //
    //one based node numbers
    //
    for(i = 0; i < 4; i++) 
      tet_conn[n][i] += 1;

    //fprintf(fp, "%9d       %8d %8d %8d %8d %8d %8d %8d %8d %8d %4d\n", n + 1, tet_conn[n][0], tet_conn[n][1], tet_conn[n][2], tet_conn[n][2], tet_conn[n][3], tet_conn[n][3], tet_conn[n][3], tet_conn[n][3], tet_vc[n], 1);
    if(tet_vc != NULL) 
      fprintf(fp, cel_format_string, n + 1, tet_conn[n][0], tet_conn[n][1], tet_conn[n][2], tet_conn[n][2], tet_conn[n][3], tet_conn[n][3], tet_conn[n][3], tet_conn[n][3], tet_vc[n], 1);
    else
      fprintf(fp, cel_format_string, n + 1, tet_conn[n][0], tet_conn[n][1], tet_conn[n][2], tet_conn[n][2], tet_conn[n][3], tet_conn[n][3], tet_conn[n][3], tet_conn[n][3], -1, 1);
  }

  int nelem = ntet;

  for(n = 0; n < npyr; n++) 
  {
    //
    //one based node numbers
    //
    for(i = 0; i < 5; i++) 
      pyr_conn[n][i] += 1;

    //fprintf(fp, "%9d       %8d %8d %8d %8d %8d %8d %8d %8d %8d %4d\n", nelem + n + 1, pyr_conn[n][0], pyr_conn[n][1], pyr_conn[n][2], pyr_conn[n][3], pyr_conn[n][4], pyr_conn[n][4], pyr_conn[n][4], pyr_conn[n][4], pyr_vc[n], 1);
    if(pyr_vc != NULL) 
      fprintf(fp, cel_format_string, nelem + n + 1, pyr_conn[n][0], pyr_conn[n][1], pyr_conn[n][2], pyr_conn[n][3], pyr_conn[n][4], pyr_conn[n][4], pyr_conn[n][4], pyr_conn[n][4], pyr_vc[n], 1);
    else
      fprintf(fp, cel_format_string, nelem + n + 1, pyr_conn[n][0], pyr_conn[n][1], pyr_conn[n][2], pyr_conn[n][3], pyr_conn[n][4], pyr_conn[n][4], pyr_conn[n][4], pyr_conn[n][4], -1, 1);
  }

  nelem += npyr;

  for(n = 0; n < npri; n++) 
  {
    //
    //one based node numbers
    //
    for(i = 0; i < 6; i++) 
      pri_conn[n][i] += 1;

    //fprintf(fp, "%9d       %8d %8d %8d %8d %8d %8d %8d %8d %8d %4d\n", nelem + n + 1, pri_conn[n][0], pri_conn[n][1], pri_conn[n][2], pri_conn[n][2], pri_conn[n][3], pri_conn[n][4], pri_conn[n][5], pri_conn[n][5], pri_vc[n], 1);
    if(pri_vc != NULL) 
      fprintf(fp, cel_format_string, nelem + n + 1, pri_conn[n][0], pri_conn[n][1], pri_conn[n][2], pri_conn[n][2], pri_conn[n][3], pri_conn[n][4], pri_conn[n][5], pri_conn[n][5], pri_vc[n], 1);
    else
      fprintf(fp, cel_format_string, nelem + n + 1, pri_conn[n][0], pri_conn[n][1], pri_conn[n][2], pri_conn[n][2], pri_conn[n][3], pri_conn[n][4], pri_conn[n][5], pri_conn[n][5], -1, 1);

  }

  nelem += npri;

  for(n = 0; n < nhex; n++) 
  {
    //
    //one based node numbers
    //
    for(i = 0; i < 8; i++) 
      hex_conn[n][i] += 1;

    //fprintf(fp, "%9d       %8d %8d %8d %8d %8d %8d %8d %8d %8d %4d\n", nelem + n + 1, hex_conn[n][0], hex_conn[n][1], hex_conn[n][2], hex_conn[n][3], hex_conn[n][4], hex_conn[n][5], hex_conn[n][6], hex_conn[n][7], hex_vc[n], 1);
    if(hex_vc != NULL) 
      fprintf(fp, cel_format_string, nelem + n + 1, hex_conn[n][0], hex_conn[n][1], hex_conn[n][2], hex_conn[n][3], hex_conn[n][4], hex_conn[n][5], hex_conn[n][6], hex_conn[n][7], hex_vc[n], 1);
    else
      fprintf(fp, cel_format_string, nelem + n + 1, hex_conn[n][0], hex_conn[n][1], hex_conn[n][2], hex_conn[n][3], hex_conn[n][4], hex_conn[n][5], hex_conn[n][6], hex_conn[n][7], -1, 1);
  }

  fclose(fp);

  //
  //write boundary element information to .bnd file
  //
  printf("\nWriting boundary elements to <%s>\n", bnd_file);

  if((fp = fopen(bnd_file, "w")) == NULL)
  {
    printf("Error opening boundary file <%s>\n", bnd_file);
    return -1;
  }

  int face_elem = 0;

  printf("\n\tNumber of boundaries = %d\n", nb);

  char bnd_format_string[] = "%8d       %8d %8d %8d %8d %6d %6d      %s\n";

  for(n = 0; n < nb; n++) 
  {
    printf("\n\tBoundary %d is %s\n", n + 1, b_name[n]);
    printf("\t\tNumber of triangles = %d\n", ntri[n]);

    //
    //Triangle faces
    //
    for(i = 0; i < ntri[n]; i++) 
    {
      //
      //one based node numbers
      //
      for(j = 0; j < 3; j++) 
        tri_conn[n][i][j] += 1;

      //fprintf(fp, "%8d       %8d %8d %8d %8d %6d %6d      %s\n", ++face_elem, tri_conn[n][i][0], tri_conn[n][i][1], tri_conn[n][i][2], tri_conn[n][i][2], n + 1, 0, b_name[n]);
      fprintf(fp, bnd_format_string, ++face_elem, tri_conn[n][i][0], tri_conn[n][i][1], tri_conn[n][i][2], tri_conn[n][i][2], n + 1, 0, b_name[n]);
    }

    printf("\t\tNumber of quads = %d\n", nquad[n]);
    //
    //Quad faces
    //
    for(i = 0; i < nquad[n]; i++) 
    {
      //
      //one based node numbers
      //
      for(j = 0; j < 4; j++) 
        quad_conn[n][i][j] += 1;

      //fprintf(fp, "%8d       %8d %8d %8d %8d %6d %6d      %s\n", ++face_elem, quad_conn[n][i][0], quad_conn[n][i][1], quad_conn[n][i][2], quad_conn[n][i][3], n + 1, 0, b_name[n]);
      fprintf(fp, bnd_format_string, ++face_elem, quad_conn[n][i][0], quad_conn[n][i][1], quad_conn[n][i][2], quad_conn[n][i][3], n + 1, 0, b_name[n]);
    }
  }

  fclose(fp);

  printf("\nWriting <%s> complete\n", fname);

  return 0;
}
