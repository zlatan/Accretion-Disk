#include <stdio.h>
#include <math.h>


double deltag( double qy, double qz)
{
  return -M_PI/2 - atan2(M_PI,2*sqrt(qy*qy + qz*qz) );
}

double deltau( double qy, double qz)
{

return -M_PI/2;
}



int main(void)
{


double qz=0;
double ID=0;
double DQ;
int i,N,N1,iter,Niter=20;

N=100;
N1=N+1;

double D[N1][N1],D_i[N1][N1],phi[N1][N1],phi_i[N1][N1],c[N1][N1],s[N1][N1],q[N1][N1],dg[N1][N1];
double du[N1][N1],sg[N1][N1],cg[N1][N1],su[N1][N1],cu[N1][N1],sug[N1][N1];
double phinew[N1][N1], Dnew[N1][N1], X[N1][N1],Y[N1][N1],sig,siu,G[N1][N1];
double Q,P,XQ,YQ,x,y,aN,rN=1.,rNnew,err1,rN2=1.,zN2,zN;
double d_reg,err2;


double C_D=0.0;

int iy,iz;
int jy,jz;
int ky,kz;
double Ly=3,Lz=3;
double DQy,DQz,Q23;
double Qy,Qz,Ky,Kz;
int my=31,mz=31;
int Ny,Nz;



Ny=2*my+1;
Nz=2*mz+1;
DQy=2*Ly/(Ny-1);
DQz=2*Lz/(Nz-1);

/*for(XQ=0;XQ<=10;XQ=XQ+0.1)
{
 for(YQ=0;YQ<=10;YQ=YQ+0.1)
 {
  printf("%f %f %f\n",XQ,YQ,deltag(XQ,YQ) );
  }
}
*/


 for(iy=1;iy<=Ny;iy++)
  for(iz=1;iz<=Nz;iz++)
	{
        Qy=-Ly+(iy-1)*DQy;
        Qz=-Lz+(iz-1)*DQz;
	dg[iy][iz]=deltag(Qy,Qz);

        sg[iy][iz]=sin(dg[iy][iz]);
        cg[iy][iz]=cos(dg[iy][iz]);
        du[iy][iz]=deltau(Qy,Qz);
        su[iy][iz]=sin(du[iy][iz]);
        cu[iy][iz]=cos(du[iy][iz]);
        sug[iy][iz]=sin(du[iy][iz]-dg[iy][iz]);

	D[iy][iz]=rand();
	phi[iy][iz]=2.*M_PI*rand();
        c[iy][iz]=cos(phi[iy][iz]);
        s[iy][iz]=sin(phi[iy][iz]);
	}


for(iter=0;iter<=Niter;iter++)
{

for(iy=1;iy<=Ny;iy++)
{
  Qy=-Ly+(iy-1)*DQy;
  //printf("Qy=%f iy=%i\n",Qy,iy);
  if(iy==my+1)  continue; 
  if(Qy==0.)  continue; 
  Q23=pow(fabs(Qy),2./3);
 
for(iz=1;iz<=Nz;iz++)
 {
  Qz=-Lz+(iz-1)*DQz;
   X[iy][iz]=0.0; 
   Y[iy][iz]=0.0;
  for(jy=1;jy<=Ny;jy++)
   {
   Ky=-Ly+(jy-1)*DQy;
   ky=iy-jy+my+1;
   if(ky<0 && ky>Ny) continue; 
   if (Ky < -Ly && Ky > Ly ) continue; 
   for(jz=1;jz<=Nz;jz++)
     {
   Kz=-Lz+(jz-1)*DQz;
   kz=iz-jz+mz+1;
   if(kz<0 && kz>Nz) continue; 
   if (Kz < -Lz && Kz > Lz ) continue; 
   X[iy][iz]+=(Ky+Kz)*D[jy][jz]*D[ky][kz]*(c[jy][jz]*c[ky][kz]+s[jy][jz]*s[ky][kz])*DQy*DQz/Q23;
   Y[iy][iz]+=(Ky+Kz)*D[jy][jz]*D[ky][kz]*(s[jy][jz]*c[ky][kz]-c[jy][jz]*s[ky][kz])*DQy*DQz/Q23;
     }// next jz
   }// next jy
 }// next iz
}// next iy

for(iy=1;iy<=Ny;iy++)
{
 for(iz=1;iz<=Nz;iz++)
 {      
           D_i[iy][iz]=sqrt(X[iy][iz]*X[iy][iz]+Y[iy][iz]*Y[iy][iz]);// Dnew_initial
           phi_i[iy][iz]=atan2(Y[iy][iz],X[iy][iz]);
           sig=sin(phi_i[iy][iz]-dg[iy][iz]);
           siu=sin(phi_i[iy][iz]-du[iy][iz]);
           
	   x=sig*cos(du[iy][iz])+sin(phi_i[iy][iz]-du[iy][iz])*cos(dg[iy][iz]);
           y=sin(phi_i[iy][iz]-dg[iy][iz])*sin(du[iy][iz])+siu*sin(dg[iy][iz]);
           phi[iy][iz]=atan2(y,x);       
	   aN=(sig*su[iy][iz]+siu*sg[iy][iz])*(sig*su[iy][iz]+siu*sg[iy][iz])
             +(sig*cu[iy][iz]+siu*cg[iy][iz])*(sig*cu[iy][iz]+siu*cg[iy][iz]); // a Numerator
           //G[iy][iz]=sqrt(aN)/fabs(sug[iy][iz]); 
	 Qy=-Ly+(iy-1)*DQy;
	 Qz=-Lz+(iz-1)*DQz;
	   G[iy][iz]=(1)/(Qy*Qy + Qz*Qz + 1);
	//  G[iy][iz]=1/(Qy*Qy + Qz*Qz +1);
	   Dnew[iy][iz]=G[iy][iz]*D_i[iy][iz];
           //printf("%f\n",siu);
	   //if ( iz == 22 ) Dnew[iy][22]=0.;
           //printf("%i %i %f\n",iy,iz,Dnew[iy][iz]);               
	  

 }// next iz
}// next iy

        
  	rN2=0.;// real norma D	
        for(iy=1;iy<=Ny;iy++)
        {
	 for(iz=1;iz<=Nz;iz++)
         {
         rN2+=Dnew[iy][iz]*Dnew[iy][iz]*DQy*DQz;
	 }
	 }
        

        rN=sqrt(rN2); 

	 for(iy=1;iy<=Ny;iy++)
         {
	 for(iz=1;iz<=Nz;iz++)
              {
              Dnew[iy][iz]/=rN;
              }// Normalization. Now |Dnew|=1
          }
        

	 for(iy=1;iy<=Ny;iy++)
         {
	 for(iz=1;iz<=Nz;iz++) {
                        D[iy][iz]=Dnew[iy][iz];
                        c[iy][iz]=cos( phi[iy][iz]);
                        s[iy][iz]=sin( phi[iy][iz]);
				}// getting old
	}
                        

                        // zN=? 
                        if(iter==Niter-1) {zN=rN; 
					for(iy=1;iy<=Ny;iy++)
				         {
	 					for(iz=1;iz<=Nz;iz++)
						{ 
						D[iy][iz]/=zN; D_i[iy][iz]=D[iy][iz]/G[iy][iz]; 
				 
		//printf("%i %i %g\n",iy,iz,D[iy][iz]);

						}// |D|=1/zN
					}
					    }


} // next iter
			
 
for(iy=1;iy<=Ny;iy++)
         {
	Qy=-Ly+(iy-1)*DQy;
	 for(iz=1;iz<=Nz;iz++) 
  	 {
          C_D+=fabs(Qy)*D[iy][iz]*D[iy][iz];
          C_D=C_D*DQy*DQz/4;
	  printf("%i %i %g\n",iy,iz,fabs(Qy)*D[iy][iz]*D[iy][iz]);
	  //printf("%i %i %g\n",iy,iz,G[iy][iz]);
	  }
	}
//printf("%g\n",C_D);



return 0;
}
