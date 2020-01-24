#include <stdio.h>
#include <stdlib.h>
#include "user.h"

int main(void){
  
  srand(1);

  FILE *fout = fopen("single_halo.raw","wb");
  int i, j, k, counter = 0, np = 4*4*4 + 120;
  realT *x_buf, *y_buf, *z_buf;
  realT x, y, z;
  fwrite(&np,sizeof(int),1,fout);

  x_buf = (realT *)malloc(np*sizeof(realT));
  y_buf = (realT *)malloc(np*sizeof(realT));
  z_buf = (realT *)malloc(np*sizeof(realT));

  /* uniform grid */
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

  /* blob around (0.25, 0.25, 0.25) */
  for(j = 1; j < 7; j++){
    for(i = 0; i < 20; i++){
      x = 0.25 + j*0.02*(1.0 - 2.0*rand()/((realT)RAND_MAX));
      y = 0.25 + j*0.02*(1.0 - 2.0*rand()/((realT)RAND_MAX));
      z = 0.25 + j*0.02*(1.0 - 2.0*rand()/((realT)RAND_MAX));
      x_buf[counter] = x;
      y_buf[counter] = y;
      z_buf[counter++] = z;
    }
  }

  fwrite(x_buf,sizeof(realT),np,fout);
  fwrite(y_buf,sizeof(realT),np,fout);
  fwrite(z_buf,sizeof(realT),np,fout);


  fclose(fout);

  /* velocities */
  fout = fopen("single_halo_vel.raw","wb"); 
  fwrite(&np,sizeof(int),1,fout);
  for(i = 0; i < np; i++){
    x_buf[i] = 0.0;
    y_buf[i] = 0.0;
    z_buf[i] = 0.0;
  }

  /* crank up the guard particles' kinetic energy */
  for(i = 0; i < 64; i++){
    x_buf[i] = 100.0*(1.0 - 2.0*rand()/((realT)RAND_MAX));
    y_buf[i] = 100.0*(1.0 - 2.0*rand()/((realT)RAND_MAX));
    z_buf[i] = 100.0*(1.0 - 2.0*rand()/((realT)RAND_MAX));
  }

  fwrite(x_buf,sizeof(realT),np,fout);
  fwrite(y_buf,sizeof(realT),np,fout);
  fwrite(z_buf,sizeof(realT),np,fout);

  fclose(fout);

  free(x_buf);
  free(y_buf);
  free(z_buf);

  return 0;
}
