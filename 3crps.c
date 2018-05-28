#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eqnl.h"
#include "2dstruc.h"


//Constantes;

#define massa_1 1
#define massa_2 2.43e-6
#define massa_3 3e-6

#define G 4*M_PI*M_PI

#define massa_t massa_1 + massa_2+massa_3

//Casting das forças
double sinal(double x);

double forca_x(double x, double y);

double forca_y(double x, double y);

main(){
  int i,j, N = 1000000;

  
  struct vec r[3],v[3],RCM,VCM;
  double dt = 0.000005;
  
  //Condições iniciais.
  r[0].x= 0;
  r[1].y= 0.7;
  r[2].y = -1;
  r[0].y=r[1].x=r[2].x = 0;
  
  v[0].y = 0;
  v[1].x = 7;
  v[2].x = -5.75;
  v[0].x = v[1].y = v[2].y = 0;

  //Mudando para referencial do centro de Massa
  RCM.x = r[0].x + r[1].x + r[2].x; RCM.x = RCM.x/massa_t;
  RCM.y = r[0].y + r[1].y + r[2].y; RCM.y = RCM.y/massa_t;
  VCM.x = v[0].x + v[1].x + v[2].x; VCM.x = VCM.x/massa_t;
  VCM.y = v[0].y + v[1].y + v[2].y; VCM.y = VCM.y/massa_t;
  


  for(i=0; i<3; i++){
  r[i].x = r[i].x - RCM.x;
  r[i].y = r[i].y - RCM.y;
  v[i].x = v[i].x - VCM.x;
  v[i].y = v[i].y - VCM.y;
  
  }
  
  //Transformação para posiçao relativa;
  struct vec s[3], vs[3];
  
  s[0].x = r[2].x - r[1].x;
  s[0].y = r[2].y - r[1].y;
  vs[0].x = v[2].x - v[1].x;
  vs[0].y = v[2].y - v[1].y;
 
  
  s[1].x = r[0].x - r[2].x;
  s[1].y = r[0].y - r[2].y;
  vs[1].x = v[0].x - v[2].x;
  vs[1].y = v[0].y - v[2].y;

  
  s[2].x = r[1].x - r[0].x;
  s[2].y = r[1].y - r[0].y;
  vs[2].x = v[1].x - v[0].x;
  vs[2].y = v[1].y - v[0].y;

  double ***solucao;
  //Alocação do solucao
  solucao = malloc(3*sizeof(double **));
  for(i=0; i<3; i++){
    solucao[i] = malloc(4*sizeof(double *));
    for(j=0; j<4; j++){
      solucao[i][j] = malloc(N*sizeof(double));
    }
  }

  //Resolve s1,s2,s3 e joga em solucao
  for(i=0;i<3;i++){
    solucao[i] = Newton_2D(forca_x,forca_y, r[i], v[i], N, dt);
  }

  //Destransforma e imprime 
  for(j=0; j<N; j++){
    
    r[0].x = massa_2*solucao[2][0][j]-massa_3*solucao[1][0][j];
    r[0].x =  r[0].x/massa_t;
    r[0].y = massa_2*solucao[2][1][j]-massa_3*solucao[1][1][j];
    r[0].y =  r[0].y/massa_t;

    r[1].x = massa_3*solucao[0][0][j] - massa_1*solucao[2][0][j];
    r[1].x = r[1].x/massa_t;
    r[1].y = massa_3*solucao[0][1][j] - massa_1*solucao[2][1][j];
    r[1].y = r[1].y/massa_t;

    r[2].x = massa_1*solucao[1][0][j] - massa_2*solucao[0][0][j];
    r[2].x = r[2].x/massa_t;
    r[2].y = massa_1*solucao[1][1][j] - massa_2*solucao[0][1][j];
    r[2].y = r[2].y/massa_t;

    if(j%10 ==0){
      for(i=0; i<3; i++){
	
	printf("%f %f\t",r[i].x, r[i].y);
      }
      printf("\t %f \n", j*dt);
    }

}
  free(solucao);
}

double sinal(double x){
  if(x>0)
    return -1;
  else
    return 1;
}


double forca_x(double x, double y){
  double K, cos_angulo=0, norm;
  K = G*massa_t;
  norm = sqrt(x*x + y*y);
  if(norm > 0.01)
    cos_angulo = fabs(x/norm);
  
  return sinal(x)*cos_angulo*K/(norm*norm);
  
}


double forca_y(double x, double y){
  double K, cos_angulo, norm,sen_angulo = 0;
  K = G*massa_t;
  norm = sqrt(x*x + y*y);
  if(norm > 0.01){
    cos_angulo = fabs(x/norm);
    sen_angulo = sqrt(1 - cos_angulo*cos_angulo);
  }
  
  return sinal(y)*sen_angulo*K/(norm*norm);
}
