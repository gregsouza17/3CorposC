#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eqnl.h"
#include "2dstruc.h"

#define massa_1 1
#define massa_2 0.00095
#define  K  4*M_PI*M_PI*massa_1*massa_2
#define M_total  massa_1 + massa_2
#define m_reduzida massa_1*massa_2/M_total
  
double sinal(double x);
double forca_x(double x, double y);
double forca_y(double x, double y);
double ext_x(double x, double y);
double ext_y(double x, double y);


void main(void){
  int i, N = 100000;
  struct vec r1,v1,r2,v2,r,v,R_CM , V_CM,rnext,vnext;
  struct vec K1,K2,K3,K4;
  double dt = 0.00001,  normr;
  double **dist_rela, **pos_CM;
  
    dist_rela = malloc(4*sizeof(double *));
    pos_CM = malloc(4*sizeof(double *));
    for(i=0; i<4; i++){
      dist_rela[i] = malloc( N*sizeof(double));
      pos_CM[i] = malloc(N*sizeof(double));
    }

  
  //Condições Iniciais de posição
  r1.x = r1.y = r2.x= 0.0; 
  r2.y = 5.2;
  //Condições iniciais de velocidade
  v1.x = v1.y = v2.y = 0.0; 
  v2.x = 2.64;
  
  //centro de massa
  R_CM.x = (r1.x*massa_1 + r2.x*massa_2)/M_total;
  R_CM.y = (r1.y*massa_1 + r2.y*massa_2)/M_total;
  V_CM.x = (v1.x*massa_1 + v2.x*massa_2)/M_total;
  V_CM.y = (v1.y*massa_1 + v2.y*massa_2)/M_total;

  //Indo pro referencial do centro de massa
  r1.x = r1.x-R_CM.x;
  r1.y = r1.y - R_CM.y;
  r2.x = r2.x-R_CM.x;
  r2.y = r2.y - R_CM.y;

  
  v1.x = v1.x-V_CM.x;
  v1.y = v1.y - V_CM.y;
  v2.x = v2.x-V_CM.x;
  v2.y = v2.y - V_CM.y;

  
  r.x = r1.x - r2.x; 
  r.y = r1.y - r2.y; 
  v.x = v1.x - v2.x;
  v.y = v1.y - v2.y;
  

  

  dist_rela = Newton_2D(forca_x,forca_y, r, v, N, dt);
  

  for(i=0; i<N; i++){
    if(i%10==0){
      printf("%f %f %f\n", i*dt, R_CM.x - massa_1*dist_rela[0][i]/M_total
	     , R_CM.y - massa_1*dist_rela[1][i]/ M_total);
    }          
  }
  free(dist_rela);
  free(pos_CM);
}


double sinal(double x){
  if(x>0)
    return -1;
  else
    return 1;
}


double forca_x(double x, double y){
  double cos_angulo=1, dist;
  dist = sqrt(x*x + y*y);
  
  if(dist!=0)
    cos_angulo = fabs(x)/dist;
  

  
  return sinal(x)* K*cos_angulo/(dist*dist*m_reduzida);
  //return sinal(x)*dist*cos_angulo;
}


double forca_y(double x, double y){
  double cos_angulo, dist, sen_angulo=1;
  dist = sqrt(x*x + y*y);
  
  if(dist!=0){
    cos_angulo = x/dist;
    sen_angulo = sqrt(1 - cos_angulo*cos_angulo);
  }
  
  return sinal(y)*K*sen_angulo/(dist*dist*m_reduzida);
  
  //return sinal(y)*dist*sen_angulo;
}

