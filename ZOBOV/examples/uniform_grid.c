#include <stdio.h>
#include <stdlib.h>
#include "user.h"

int main(void){

  srand(0);

  FILE *fout = fopen("uniform_grid.raw","wb");
  int i, j, k, counter = 0, np = 4*4*4;
  realT *x_buf, *y_buf, *z_buf;
  realT x, y, z;
  fwrite(&np,sizeof(int),1,fout);

  x_buf = (realT *)malloc(np*sizeof(realT));
  y_buf = (realT *)malloc(np*sizeof(realT));
  z_buf = (realT *)malloc(np*sizeof(realT));

  for(i = 0; i < 4; i++){
    x = ((realT)i)/4.0+1.0/8.0;
    for(j = 0; j < 4; j++){
      y = ((realT)j)/4.0+1.0/8.0;
      for(k = 0; k < 4; k++){
	z = ((realT)k)/4.0+1.0/8.0;
	x_buf[counter] = x;
	y_buf[counter] = y;
	z_buf[counter++] = z;
      }
    }    
  }

  fwrite(x_buf,sizeof(realT),np,fout);
  fwrite(y_buf,sizeof(realT),np,fout);
  fwrite(z_buf,sizeof(realT),np,fout);

  free(y_buf);
  free(z_buf);

  fclose(fout);

  /* velocities */
  fout = fopen("uniform_grid_vel.raw","wb"); 
  fwrite(&np,sizeof(int),1,fout);
  for(i = 0; i < np; i++){
    x_buf[i] = 0.0;
  }

  fwrite(x_buf,sizeof(realT),np,fout);
  fwrite(x_buf,sizeof(realT),np,fout);
  fwrite(x_buf,sizeof(realT),np,fout);
  fclose(fout);

  return 0;
}
