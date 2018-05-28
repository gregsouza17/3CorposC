#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct vec {
  double x;
  double y;
} vec;


double norm(struct vec vet){
  return sqrt(vet.x*vet.x + vet.y*vet.y);
}


double **Newton_2D(double (*forca_x)(), double (*forca_y)(),
		   struct vec r_inic, struct vec v_inic, int N, double dt){
  
  int i;
  struct vec r,rnext,v,vnext,velo1,velo2,velo3,velo4, pos1,pos2,pos3,pos4;
  
  double **solucao,t,normr;
  
  //Aloca espaço de solucao
  solucao = malloc(4*sizeof(double *));
  
  for(i=0; i<4; i++){
    solucao[i] = malloc( N*sizeof(double));
    //printf("%d\n",i);
  }
  

  //Aplica condições iniciais
  r.x = r_inic.x;
  r.y = r_inic.y;
  v.x = v_inic.x;
  v.y = v_inic.y;
  
  
   i=0;
   
  while(i<N){
  	//Here i'm setting the solution, the first two rows are de x and y position, and the third and fourth are the velocitys
    solucao[0][i] = r.x;
    solucao[1][i] = r.y;
    solucao[2][i] = v.x;
    solucao[3][i] = v.y;
 
      
    t = i*dt;
    
    
    velo1.x = dt*forca_x(r.x,r.y);
    pos1.x = dt*v.x ;
    
    velo1.y = dt*  forca_y(r.x,r.y);
    pos1.y = v.y*dt;
    

    velo2.x = dt*forca_x(r.x +  0.5*pos1.x, r.y  + 0.5*pos1.y  );
    pos2.x = (v.x + velo1.x/2)*dt/2;

    
    velo2.y = dt*forca_y(r.x +  0.5*pos1.x, r.y  + 0.5*pos1.y  );
    pos2.y = (v.y +  velo1.y/2    )*dt;
 
    
    velo3.x =  dt*forca_x(r.x +  0.5*pos2.x, r.y  + 0.5*pos2.y);
    pos3.x = (v.x +  velo2.x/2    )  *dt;
    
    velo3.y =  dt*forca_y(r.x +  0.5*pos2.x, r.y  + 0.5*pos2.y);
    pos3.y = (v.y +  velo2.y/2     )  *dt;
   
    
    velo4.x =  dt*forca_x(r.x + pos3.x, r.y + pos3.y);
    pos4.x = (v.x +  velo3.x     )*dt;
    
    velo4.y =  dt*forca_y(r.x + pos3.x, r.y + pos3.y);
    pos4.y = (v.y +   velo3.y      )*dt;
  
    
    vnext.x  = v.x  + (velo1.x + velo2.x + velo3.x + velo4.x)/6.0 ;
    rnext.x  = r.x + (pos1.x + pos2.x + pos3.x + pos4.x)/6.0 ;
    
    vnext.y  = v.y  + (velo1.y + velo2.y + velo3.y + velo4.y)/6.0 ;
    rnext.y  = r.y + (pos1.y + pos2.y + pos3.y + pos4.y)/6.0 ;

    // if(i%1000 == 0)
    //printf("%e %e %e\n", r.x, r.y, t);
      //scanf("*");

    r.x = rnext.x;
    r.y = rnext.y;
    
    v.x = vnext.x;
    v.y = vnext.y;

    
    i++;

   
      // printf("%f %f %f\n",i*dt,solucao[1][i],solucao[2][i]);
  
    
  }
  return solucao;

  free(solucao);
}



