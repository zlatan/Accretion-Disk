#include <stdio.h>
#include <math.h>


double Q0=M_PI/2;
int D=3;

const double wc=0.;
const double theta=M_PI/2;

void force(double f[D],double y[D], double t, double K_y, double K_z)
{

f[1] = y[2];
f[2] = -( (K_y*K_y + K_z*K_z) -(1/(1+2*t*t +t*t*t*t)))*y[1];
}


double faza(double K_y, double K_z,double initialg, double initialu)
{

double a,b,y0[D],h2,k1[D],k2[D],k3[D],k4[D],dy[D],h,y[D],yc[D],t;
int N,n,d,i=0;
double Q,xi,tau,p[10],tau_a[10];
double cc=0,cs=0,ss=0,cp=0,sp=0,dd=0,d2=0;
double c,det,aa,bb,f,tt,siika;

Q=sqrt(K_y*K_y+K_z*K_z);

a=0;//*Q/K_y;
b=100;//*Q/K_y;
N=1000;
h=(b-a)/(N-1);
h2=h/2;

y0[1]=initialg;//
y0[2]=initialu;//

t=a;
xi=a*K_y/Q;
for(d=1;d<=D;d++) {y[d]=y0[d];}

for(n=1;n<N;n++)
{
t=a+h*(n-1);// t_n
force(k1,y,t,K_y,K_z);
for(d=1;d<=D;d++) yc[d]=y[d]+h2*k1[d];
force(k2,yc,t+h2,K_y,K_z); 
for(d=1;d<=D;d++) yc[d]=y[d]+h2*k2[d];
force(k3,yc,t+h2,K_y,K_z); 
for(d=1;d<=D;d++) yc[d]=y[d]+h*k3[d];
force(k4,yc,t+h,K_y,K_z);

for(d=1;d<=D;d++)  y[d]+=h*(k1[d]+k4[d]+(k2[d]+k3[d])*2)/6; 

xi=(t+h)*K_y/Q;
tau=t+h;


 if ( i < 7 ) 
			{			
		if ( t>30.0 )
			{						
			p[i]=y[1];				
			tau_a[i] = t;
			i++;
			}
		}


}


for(i=0;i<7;i++)
{
tt=tau_a[i]*K_z;
c=cos(tt); 
siika=sin(tt);
cc=cc + c*c;
cs=cs + c*siika;
ss=ss + siika*siika;
sp=sp + siika*p[i];
cp=cp + c*p[i];
}

det=cc*ss-cs*cs;
aa=(ss*cp-cs*sp)/det;
bb=(cs*cp-cc*sp)/det;
dd=sqrt(aa*aa+bb*bb);
f=atan2(bb,aa);

 //printf("=%g\n",f);

return f;

}


double deltag( double qy, double qz)
{
return ( faza(qy,qz,0.0,1.0) );
}

double deltau( double qy, double qz)
{
return ( faza(qy,qz,1.0,0.0) );
}



int main(void)
{

double qqq=0;
double ID=0;
double DQ,Qmax=5.;
int N,N1,i,j,ja,k,ka,iter,Niter=41;
	N=1000;
	N1=N+1;

double D[N1],D_i[N1],phi[N1],phi_i[N1],c[N1],s[N1],q[N1],dg[N1],du[N1],sg[N1],cg[N1],su[N1],cu[N1],sug[N1];
double phinew[N1], Dnew[N1];
double Q,P,X,Y,XQ,YQ,x,y,sig,siu,aN,rN=1.,rNnew,err1,rN2=1.,G[N1],zN2,zN;
DQ=Qmax/N;
double d_reg,err2;

for (qqq=-13.7;qqq<13.7;qqq=qqq+0.08)
{
 for(i=1;i<=N;i++)
	{
	q[i]=DQ*i;
	dg[i]=deltag(q[i],qqq);
        sg[i]=sin(dg[i]);
        cg[i]=cos(dg[i]);
        du[i]=deltau(q[i],qqq);
        su[i]=sin(du[i]);
        cu[i]=cos(du[i]);
        sug[i]=sin(du[i]-dg[i]);
	}

 for(i=1;i<=N;i++)
	{
	D[i]=rand();
	phi[i]=2.*M_PI*rand();
        c[i]=cos(phi[i]);
        s[i]=sin(phi[i]);
	}
        
for(iter=0;iter<=Niter;iter++)
{
	for(i=1;i<=N;i++)
           {// tova e cycler po Q
	   Q=q[i];  
         XQ=0.; YQ=0.;
	 for(j=-N;j<=N;j++)
           {// tova e cycle po P in (-infty,+infty)
            ja=abs(j);
            P=DQ*j;
            k=i-j; ka=abs(k);//QP=Q-P=DQ*(i-j)=DQ*k 
            if(ka>N) continue;
            if(ka==0) continue; // i,ja,ka in [1,N]]
            if(ja==0) continue; // i,j,k in [-N, -1] U [1, N]
            XQ+=P*D[ja]*D[ka]*(c[ja]*c[ka]+s[ja]*s[ka]);//*DQ
            YQ+=P*D[ja]*D[ka]*(s[ja]*c[ka]-c[ja]*s[ka]);//*DQ
           }// next j
           d_reg=pow(fabs(Q),2./3);
           //d_reg=1.;
           X=(XQ*DQ)/d_reg;//*DQ
           Y=(YQ*DQ)/d_reg;//*DQ
           D_i[i]=sqrt(X*X+Y*Y);// Dnew_initial
           phi_i[i]=atan2(Y,X);
           sig=sin(phi_i[i]-dg[i]);
           siu=sin(phi_i[i]-du[i]);
           
	   x=sig*cos(du[i])+sin(phi_i[i]-du[i])*cos(dg[i]);
           y=sin(phi_i[i]-dg[i])*sin(du[i])+siu*sin(dg[i]);
           phi[i]=atan2(y,x);       
	   aN=(sig*su[i]+siu*sg[i])*(sig*su[i]+siu*sg[i])
             +(sig*cu[i]+siu*cg[i])*(sig*cu[i]+siu*cg[i]); // a Numerator
           G[i]=sqrt(aN)/fabs(sug[i]); 
	   Dnew[i]=G[i]*D_i[i];
	}// next i  
	
        rN2=0.;// real norma D	
        for(i=1;i<=N;i++) {rN2+=Dnew[i]*Dnew[i]*DQ;}
        rNnew=sqrt(rN2); 
        if(iter==Niter){err2=(rNnew-1/zN)*2/(rNnew+1/zN);}
        err1=(rNnew-rN)*2./(rNnew+rN);
        rN=rNnew;

        for(i=1;i<=N;i++) {Dnew[i]/=rN;}// Normalization. Now |Dnew|=1
        for(i=1;i<=N;i++) {
                        D[i]=Dnew[i];
                        c[i]=cos( phi[i]);
                        s[i]=sin( phi[i]);
				}// getting old
                        
                        // zN=? 
                        if(iter==Niter-1) {zN=rN; for(i=1;i<=N;i++) D[i]/=zN; D_i[i]=D[i]/G[i]; }// |D|=1/zN
} // next iter
			

       double rID=0;
       for(i=1;i<=N;i++)
       {
       rID+=D[i]*D[i]*q[i];
       //printf("%g %g\n",q[i],D[i]*q[i]);
   //printf("%g %g\n",q[i],D[i]);       
	} 
       rID=(rID*DQ)/(zN*zN);
      //printf("Niter=%i rID=%g\n",Niter,rID);
ID += rID*0.08;
printf("qqq=%f -> I_D=%f\n",qqq,ID);
}

/*
	double res,Qy,Qz;

	for (Qy=-2.141;Qy<2.2;Qy+=0.01)
	   for (Qz=-3.141;Qz<3.2;Qz+=0.1)
		{
		res=inegral(Qy,0.);    	     	    
		printf ("%g %g \n",Qy,res);  
		} 


*/

return 0;
}
