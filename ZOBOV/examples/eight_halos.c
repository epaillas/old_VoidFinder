#include <stdio.h>
#include <stdlib.h>
#include "user.h"

int main(void){

  srand(8);

  FILE *fout = fopen("eight_halos.raw","wb");
  int i, j, k, l, m, counter = 0, np = 16*16*16 + 8*120;
  realT *x_buf, *y_buf, *z_buf;
  realT x, y, z;

  fwrite(&np,sizeof(int),1,fout);

  x_buf = (realT *)malloc(np*sizeof(realT));
  y_buf = (realT *)malloc(np*sizeof(realT));
  z_buf = (realT *)malloc(np*sizeof(realT));

  /* uniform grid */
  for(i = 0; i < 16; i++){
    x = ((realT)i)/16.0+1.0/32.0;
    for(j = 0; j < 16; j++){
      y = ((realT)j)/16.0+1.0/32.0;
      for(k = 0; k < 16; k++){
	z = ((realT)k)/16.0+1.0/32.0;
	x_buf[counter] = x;
	y_buf[counter] = y;
	z_buf[counter++] = z;
      }
    }    
  }

  /* blob around (0.25, 0.25, 0.25) */
  for(i = 0; i < 2; i++){
    for(j = 0; j < 2; j++){
      for(k = 0; k < 2; k++){
	for(l = 1; l < 7; l++){
	  for(m = 0; m < 20; m++){
	    x = i*0.5 + 0.25 + l*0.02*(1.0 - 2.0*rand()/((realT)RAND_MAX));
	    y = j*0.5 + 0.25 + l*0.02*(1.0 - 2.0*rand()/((realT)RAND_MAX));
	    z = k*0.5 + 0.25 + l*0.02*(1.0 - 2.0*rand()/((realT)RAND_MAX));
	    x_buf[counter] = x;
	    y_buf[counter] = y;
	    z_buf[counter++] = z;
	  }
	}
      }
    }
  }

  fwrite(x_buf,sizeof(realT),np,fout);
  fwrite(y_buf,sizeof(realT),np,fout);
  fwrite(z_buf,sizeof(realT),np,fout);

  fclose(fout);

  /* velocities */
  fout = fopen("eight_halos_vel.raw","wb"); 
  fwrite(&np,sizeof(int),1,fout);

  for(i = 0; i < np; i++){
    x_buf[i] = 0.0;
    y_buf[i] = 0.0;
    z_buf[i] = 0.0;
  }

  /* crank up the guard particles' kinetic energy */
  for(i = 0; i < 16*16*16; i++){
    x_buf[i] = 100.0*(1.0 - 2.0*rand()/((realT)RAND_MAX));
    y_buf[i] = 100.0*(1.0 - 2.0*rand()/((realT)RAND_MAX));
    z_buf[i] = 100.0*(1.0 - 2.0*rand()/((realT)RAND_MAX));
  }

  fwrite(x_buf,sizeof(realT),np,fout);
  fwrite(y_buf,sizeof(realT),np,fout);
  fwrite(z_buf,sizeof(realT),np,fout);

  fclose(fout);

  return 0;
}
