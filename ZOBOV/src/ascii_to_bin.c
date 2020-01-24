/* 
Takes as an input an ASCII file where
the first three columns correspond to
the particle positions.

Writes as an output an unformatted
binary file that is readable
by ZOBOV

E. Paillas <epaillas@astro.puc.cl>
*/

#include <stdio.h>
#include <stdlib.h>
#include "user.h"

int main(int argc, char *argv[]){

  FILE *fin = fopen(argv[1],"r");
  FILE *fout = fopen(argv[2],"wb");
  char ch;
  int i, np;
  realT *x_buf, *y_buf, *z_buf;
  realT x, y, z;

  if (argc < 3)  
    { 
        printf("Some arguments are missing.\n");
        printf("1) input_ascii\n");
        printf("2) output_bin\n");
        return 0; 
    } 

  printf("Input file: %s \n", argv[1]);
  printf("Output file: %s \n", argv[2]);

  /* Count the number of lines in the input file */
  np = 0;
  while ((ch=getc(fin)) != EOF) {
    if (ch == '\n') { ++np; }
  }
  rewind(fin);
  printf("Number of particles: %d \n\n", np);

  /* Allocate memory for position arrays */
  x_buf = (realT *)malloc(np*sizeof(realT));
  y_buf = (realT *)malloc(np*sizeof(realT));
  z_buf = (realT *)malloc(np*sizeof(realT));

  /* Write positions to arrays */
  for (i=0 ; i < np; i++){
    fscanf(fin, "%lf %lf %lf", &x_buf[i], &y_buf[i], &z_buf[i]);
  }
  printf("First particle positions:\n%lf %lf %lf \n", x_buf[1], y_buf[1], z_buf[1]);

  /* Write to output binary file. */
  fwrite(&np,sizeof(int),1,fout);
  fwrite(x_buf,sizeof(realT),np,fout);
  fwrite(y_buf,sizeof(realT),np,fout);
  fwrite(z_buf,sizeof(realT),np,fout);

  fclose(fin);
  fclose(fout);
  return 0;
}

