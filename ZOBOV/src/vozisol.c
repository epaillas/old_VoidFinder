#include <ctype.h>
#include "voz.h"

#define DL for (d=0;d<3;d++)
#define PRINTFREQ 1000
/* print out particle volume every PRINTFREQ particles */

int delaunadj (coordT *points, int nvp, int nvpbuf, int nvpall, PARTADJ **adjs);
int vorvol (coordT *deladjs, coordT *points, pointT *intpoints, int numpoints,  realT *vol);
int posread_isol(char *posfile, float ***p, float fact, int *np, int *npreal);

int main(int argc, char *argv[]) {
  int exitcode;
  int i, j, k, np, npreal;
  float **r;
  coordT rtemp[3], *parts;
  coordT deladjs[3*MAXVERVER], points[3*MAXVERVER];
  pointT intpoints[3*MAXVERVER];
  FILE *pos, *adj, *vol;
  char posfile[256], adjfile[256], volfile[256], *prefix;
  char asciiadjfile[261], asciivolfile[261];
  char asciiyn;
  PARTADJ *adjs;
  double *vols;
  float volley;
  float predict, xmin,xmax,ymin,ymax,zmin,zmax;
  
  float width, width2, totwidth, totwidth2, bf, s, g;
  float border, boxsize;
  float c[3];
  int hasfewadj;
  int numborder, numfewadj, numinside, nout,maxver;
  double totalvol;

  char d;

  if (argc != 4) {
    printf("Wrong number of arguments.\n");
    printf("arg1: file prefix\n");
    printf("arg2: boxsize\n");
    printf("arg3: Output vol, adj data in ascii? (y,n for yes,no)\n");
    exit(0);
  }
  prefix = argv[1];
  sprintf(posfile,"%s.pos",prefix);
  if (sscanf(argv[2],"%f",&boxsize) != 1) {
    printf("That's no boxsize; try again.\n");
    exit(0);
  }
  asciiyn = tolower(argv[3][0]);
  printf("%c\n",asciiyn);

  maxver = MAXVERVER;
  /* Boxsize should be the range in r, yielding a range 0-1 */

  posread_isol(posfile,&r,1./boxsize,&np,&npreal);
  printf("%d particles; %d of them real\n",np,npreal);fflush(stdout);
  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=0; i<np;i++) {
    if (r[i][0]<xmin) xmin = r[i][0]; if (r[i][0]>xmax) xmax = r[i][0];
    if (r[i][1]<ymin) ymin = r[i][1]; if (r[i][1]>ymax) ymax = r[i][1];
    if (r[i][2]<zmin) zmin = r[i][2]; if (r[i][2]>zmax) zmax = r[i][2];
  }
  printf("np: %d, x: %f,%f; y: %f,%f; z: %f,%f\n",np,xmin,xmax, ymin,ymax, zmin,zmax); fflush(stdout);

  adjs = (PARTADJ *)malloc(np*sizeof(PARTADJ));
  if (adjs == NULL) {
    printf("Unable to allocate adjs\n");
    exit(0);
  }

  parts = (coordT *)malloc(3*np*sizeof(coordT));
  
  for (i=0; i<np; i++) {
    parts[3*i] = r[i][0];
    parts[3*i+1] = r[i][1];
    parts[3*i+2] = r[i][2];

    if (r[i][0] < xmin) xmin = r[i][0];
    if (r[i][0] > xmax) xmax = r[i][0];
    if (r[i][1] < ymin) ymin = r[i][1];
    if (r[i][1] > ymax) ymax = r[i][1];
    if (r[i][2] < zmin) zmin = r[i][2];
    if (r[i][2] > zmax) zmax = r[i][2];
  }
  for (i=0;i<np;i++) free(r[i]);
  free(r);
  
  /* Do tesselation*/
  printf("File read.  Tessellating ...\n"); fflush(stdout);
  exitcode = delaunadj(parts, np, np, np, &adjs);
  
  /* Now calculate volumes*/
  printf("Now finding volumes ...\n"); fflush(stdout);
  vols = (double *)malloc(npreal*sizeof(double));
  
  numborder = 0;
  numfewadj = 0;
  printf("npreal = %d\n",npreal); fflush(stdout);

  for (i=0; i<npreal; i++) { /* Just the original particles
			     Assign adjacency coordinate array*/
    hasfewadj = 0;
    vols[i] = 0.;
    /* Volumes */
    for (j = 0; j < adjs[i].nadj; j++) {
      if ((adjs[i].adj[j] >= npreal) || (adjs[i].adj[j] < 0)) {
	vols[i] = 1e-30; /*Border point*/
	/* Get rid of the adjacency outside the boundary */
	adjs[i].nadj = j;
	if (adjs[i].adj[j] < 0) {
	  printf("Negative adj: %d %d %d\n",i,j,adjs[i].adj[j]);
	}
      }
    }
    if (adjs[i].nadj < 4) {
      if ((vols[i] == 0.)&&(numborder == 0)) {
	hasfewadj = 1;
	numfewadj ++; /* Only if it initially had too few */
      }
      vols[i] = 1e-30;
      if (adjs[i].nadj < 0) {
	printf("#adj(%d)=%d; on the boundary. Expect warning in jo?o?.\n",i,adjs[i].nadj);
      }
    }
    if (vols[i] == 0.) { /* non-border point */
      for (j = 0; j < adjs[i].nadj; j++) {
	DL {
	  deladjs[3*j + d] = parts[3*adjs[i].adj[j]+d] - parts[3*i+d];
	}
      }
      exitcode = vorvol(deladjs, points, intpoints, adjs[i].nadj, &(vols[i]));
      /*if (vols[i] > 1e10) {
	printf("Bigvol! %d: %d adjs, %g\n",i,adjs[i].nadj,vols[i]);
	}*/
      vols[i] *= (double)np;
    } else {
      if (hasfewadj == 0) /* only if it's not already counted in numfewadj */
	numborder ++;
    }
    if ((i % PRINTFREQ) == 0){
      printf("%d: %d, %g\n",i,adjs[i].nadj,vols[i]);fflush(stdout);
    }
  }

  totalvol = 0.;
  numinside=0;
  for (i=0;i<npreal; i++) {
    if (vols[i] > 1.01e-30) {
      totalvol += vols[i];
      numinside++;
    }
  }
  printf("numborder = %d, numfewadj = %d, numinside = %d\n",numborder,numfewadj, numinside);
  printf("Total %d out of %d\n",numborder+numfewadj+numinside,npreal);

  if (numborder+numfewadj+numinside != npreal) {
    printf("It doesn't add up!\n");
  }
  printf("Average volume = %g\n",totalvol/(double)numinside);
  
  /* Now the output! */

  sprintf(adjfile,"%s.adj",prefix);
  sprintf(asciiadjfile,"%s.ascii.adj",prefix);
  sprintf(volfile,"%s.vol",prefix);
  sprintf(asciivolfile,"%s.ascii.vol",prefix);

  printf("Outputting to %s, %s\n\n",adjfile,volfile);

  adj = fopen(adjfile,"w");
  if (adj == NULL) {
    printf("Unable to open %s\n",adjfile);
    exit(0);
  }
  fwrite(&npreal,1, sizeof(int),adj);

  /* Adjacencies: first the numbers of adjacencies, 
     and the number we're actually going to write per particle */
  for (i=0;i<npreal;i++) {
    if (adjs[i].nadj < 0) {
      adjs[i].nadj = 0; /* In a weird case of no boundary, it could be -1 */
    }
    fwrite(&adjs[i].nadj,1,sizeof(int),adj);
  }

  /* Now the lists of adjacencies (without double counting) */
  for (i=0;i<npreal;i++) {
    nout = 0;
    for (j=0;j<adjs[i].nadj; j++) if (adjs[i].adj[j] > i) nout++;
    fwrite(&nout,1,sizeof(int),adj);      
    for (j=0;j<adjs[i].nadj; j++) 
      if (adjs[i].adj[j] > i) 
	fwrite(&(adjs[i].adj[j]),1,sizeof(int),adj);
  }

  fclose(adj);

  /* Volumes */
  vol = fopen(volfile,"w");
  if (vol == NULL) {
    printf("Unable to open %s\n",volfile);
    exit(0);
  }
  fwrite(&npreal,1, sizeof(int),vol);
  for (i=0;i<npreal;i++) {
    volley = (float)vols[i];
    fwrite(&volley,sizeof(float),1,vol);
  }

  /* Ascii copies, if you want */
  if (asciiyn == 'y') {
    adj = fopen(asciiadjfile,"w");
    if (adj == NULL) {
      printf("Unable to open %s\n",adjfile);
      exit(0);
    }
    fprintf(adj,"%d\n",npreal);

    /* Adjacencies: first the numbers of adjacencies, 
       and the number we're actually going to write per particle */
    for (i=0;i<npreal;i++) {
      if (adjs[i].nadj < 0) {
	adjs[i].nadj = 0; /* In a weird case of no boundary, it could be -1 */
      }
      fprintf(adj,"%d %d\n",i,adjs[i].nadj);
    }
    
    /* Now the lists of adjacencies (without double counting) */
    for (i=0;i<npreal;i++) {
      nout = 0;
      for (j=0;j<adjs[i].nadj; j++) if (adjs[i].adj[j] > i) nout++;
      fprintf(adj,"%d %d: ",i,nout);
      for (j=0;j<adjs[i].nadj; j++) 
	if (adjs[i].adj[j] > i) 
	  fprintf(adj,"%d ",adjs[i].adj[j]);
      
      fprintf(adj,"\n");
    }
    fclose(adj);

    /* Volumes */
    vol = fopen(asciivolfile,"w");
    if (vol == NULL) {
      printf("Unable to open %s\n",volfile);
      exit(0);
    }
    fprintf(vol,"%d\n",npreal);
    for (i=0;i<npreal;i++) {
      volley = (float)vols[i];
      fprintf(vol,"%f\n",volley);
    }
    fclose(vol);
  }
  return(0);
}

