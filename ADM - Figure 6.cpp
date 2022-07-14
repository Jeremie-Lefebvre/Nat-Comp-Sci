
#include <iostream>
#include <string.h>
#include <fstream>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <conio.h>
#include <stdlib.h>    
#include <stdio.h>

///////////////RAN2PARAMETERS///////////////////////

#define IM1   2147483563
#define IM2   2147483399
#define AM    (1.0/IM1)
#define IMM1  (IM1-1)
#define IA1   40014
#define IA2   40692
#define IQ1   53668
#define IQ2   52774
#define IR1   12211
#define IR2   3791
#define NTAB  32
#define NDIV  (1+IMM1/NTAB)
#define EPS   1.2e-7
#define RNMX  (1.0 - EPS)
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define Pi 3.14159

///////////////PARAMETERS///////////////////////

#define T 0//200000 //Time of simulation in msec
#define N 100//number of neurons
#define dt 0.01 //time step for Euler Maruyama integration
#define w_o 0.21//synaptic connectivity weight (criticality) 
#define beta 25 //activation function non-linear gain
#define rho 0.15//connection probability
#define Omega 1000//10.0 //network spatial scale in mm
#define h 0.10//activation function threshold (inflexion point)
#define c_min 0.1 //minimal conduction velocity in m/s
#define c_max 1000//maximal conduction velocity in m/s
#define Q 5//parameters grid size

//////////////////////FUNCTION DEFINITIONS///////////////////////
double f(double u);//neuron response function
double f_prime(double u);//derivative of f
double f_doubleprime(double u);//2nd derivative of f
float ran2(long *idum);//random number generator - initial cond.

//////////////////////OBJECT DEFINITIONS///////////////////////
double x[N];//neuron position along x axis
double y[N];//neuron position along y axis
double z[N];//neuron position along z axis

double c[N][N];//conduction velocity matrix
double w[N][N];//synaptic weight matrix
double l[N][N];//axonal tract lengths matrix
int tau[N][N];//conduction delays matrix
double Velocity[Q];//conduction velocity matrix

double xi[T];//noise array
double u[T];//mean field dynamics

double u_1[T];//linearized dynamics fized point phi_1
double u_3[T];//linearized dynamics fized point phi_3

double MEAN_FIRING_RATE[Q];//mean field firing rate(theoretical)
double MEAN_VARIANCE[Q];//mean field variance(theoretical)
double MEAN_CORRELATION[Q]; //mean field correlation (theoretical)

//variance of the linearized mean activity (numerical)
double Variance_1[Q];
double Variance_3[Q];

////theoretical variances  calculated from : "Kramers-Moyal expansion for stochastic differential equations 
//with single and multiple delays: Applications to financial physics and neurophysics"
//Physics Letters A 360, 552-562 (2007).
double TD_FRANK_Variance_1[Q];
double TD_FRANK_Variance_3[Q];

using namespace std;


int main()
{
	
				for (int q=0;q<Q;q++)
				{
	
					Velocity[q] = pow(10,q-2);//0.1+c_max/((double)Q)*q;
					
					
					
					
					
//////////////////////DEFINING CONNECTIVITY, DELAYS AND CONDUCTION VELOCITIES///////////////////////
					MEAN_FIRING_RATE[q]=0;
					MEAN_VARIANCE[q]=0;
					for (int i=0;i<N;i++)
					{
						long seed1 = (long)21+i*187198+56*i+12*1+1+q*485;
						long seed2 = (long)21+i*56+56*i+11*1+13783*q;
						long seed3 = (long)21+i*789+56*i+8*1+1535*q+q;
		                          	  					
						x[i] = ran2(&seed1)*Omega;
						y[i] = ran2(&seed2)*Omega;
						z[i] = ran2(&seed3)*Omega;
					}
					
					for(int i=0;i<N;i++)
					{
						for (int j=0;j<N;j++)
						{
							w[i][j] = 0;
							l[i][j] = fabs(sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j])));
							c[i][j] = Velocity[q];
							tau[i][j] = (int) floor(l[i][j]/(c[i][j]));
							
						}
						
					}
				
	


					
					
//////////////////////EULER MARUYAMA INTEGRATION SCHEME///////////////////////
					double fixedPoint_1=0.0335;
					double fixedPoint_3=0.19;
		

//////////////////////MEAN VALUE AND VARIANCE OF THE MEAN ACTIVITY///////////////////////
		
   				    double Gamma=0;
   				    double s=0.1;
			
        			for (int i=0;i<N;i++)
        			{
        				for (int j=0;j<N;j++)
        				{
        					Gamma = Gamma+1/((double)N*N)*exp(-tau[i][j]*dt*s);
						}
					}
					

					
					double mean_firing_rate=0;
					double mean_variance =0;
					double S=w_o*w_o/((double)N*(1-w_o*f_prime(fixedPoint_1)*Gamma));   
					cout<<S<<endl;
					cout<<f_prime(fixedPoint_1)*f_prime(fixedPoint_1)*w_o*w_o*f(fixedPoint_1)<<endl;
					
					mean_firing_rate =   0.2537354531e-1 * sqrt(-(0.2383736622e1 * S + 0.6684252738e0) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.3805155286e-1 * S - 0.1923133216e-2) / S) / (0.1e1 / 0.3141592654e1 * sqrt(-(0.2383736622e1 * S + 0.6684252738e0) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.3805155286e-1 * S - 0.1923133216e-2) / S) / 0.2e1 + 0.1e1 / 0.3141592654e1 * sqrt(-(0.2101433283e2 * S + 0.1162972721e2) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.7165986182e-1 * S + 0.1904648304e-2) / S) / 0.2e1) + 0.1439796044e0 * sqrt(-(0.2101433283e2 * S + 0.1162972721e2) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.7165986182e-1 * S + 0.1904648304e-2) / S) / (0.1e1 / 0.3141592654e1 * sqrt(-(0.2383736622e1 * S + 0.6684252738e0) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.3805155286e-1 * S - 0.1923133216e-2) / S) / 0.2e1 + 0.1e1 / 0.3141592654e1 * sqrt(-(0.2101433283e2 * S + 0.1162972721e2) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.7165986182e-1 * S + 0.1904648304e-2) / S) / 0.2e1);

					
					mean_variance =  0.1406739725e-1 * sqrt(-(0.2383736622e1 * S + 0.6684252738e0) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.3805155286e-1 * S - 0.1923133216e-2) / S) * pow(0.1e1 / 0.3141592654e1 * sqrt(-(0.2383736622e1 * S + 0.6684252738e0) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.3805155286e-1 * S - 0.1923133216e-2) / S) / 0.2e1 + 0.1e1 / 0.3141592654e1 * sqrt(-(0.2101433283e2 * S + 0.1162972721e2) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.7165986182e-1 * S + 0.1904648304e-2) / S) / 0.2e1, -0.2e1) * sqrt(-(0.2101433283e2 * S + 0.1162972721e2) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.7165986182e-1 * S + 0.1904648304e-2) / S);

					
													
					MEAN_FIRING_RATE[q] = 100*mean_firing_rate;
				
					MEAN_VARIANCE[q] =100*100* mean_variance;
					MEAN_CORRELATION[q] =mean_variance/(mean_variance+0.5*f_prime(fixedPoint_1)*f_prime(fixedPoint_1)*w_o*w_o*f(fixedPoint_1));
					cout<<" "<<Variance_1[q]<<" "<<TD_FRANK_Variance_1[q]<<" "<<Variance_3[q]<<" "<<TD_FRANK_Variance_3[q]<<" "<< MEAN_FIRING_RATE[q]<<"	"<<MEAN_VARIANCE[q]<<"	"<<MEAN_CORRELATION[q]<<endl;
					cout<<100*100*f_prime(fixedPoint_1)*f_prime(fixedPoint_1)*Variance_1[q]<<endl;
				}// q loop
		
					
			
					
/////////////////////OUTPUT ROUTINE/////////////////////////////					
					
		ofstream outfile;
      	
 	
    
    	outfile.open("ADM - Mean Firing Rate.txt", ios::out);
    	for(int q=0;q<Q;q++)
    	{
	   		outfile<<Velocity[q]<<"	"<<MEAN_FIRING_RATE[q]<<endl;
    	}  
    	outfile.close(); 
    	
    	    	outfile.open("ADM - Mean Firing Rate Variance.txt", ios::out);
    	for(int q=0;q<Q;q++)
    	{
	   		outfile<<Velocity[q]<<"	"<<MEAN_VARIANCE[q]<<endl;
    	}  
    	outfile.close(); 
    	
   	
    	outfile.open("ADM - Variance of v at phi_1.txt", ios::out);
    	for(int q=0;q<Q;q++)
    	{
	   		outfile<<Velocity[q]<<"	"<<Variance_1[q]<<"	"<<TD_FRANK_Variance_1[q]<<endl;
    	}  
    	outfile.close();
    	
    	
    
    	
    	
					
		cout<<"Simulations complete..."<<endl;
      
return 0;    
}


///////////////FUNCTION DEFINITIONS///////////////////////

float ran2(long *idum)
{
  int j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0) {                             /* initialize */
    if (-(*idum) < 1)                           /* prevent idum == 0 */
      *idum = 1;
    else
      *idum = -(*idum);                         /* make idum positive */
    idum2 = (*idum);
    for (j = NTAB + 7; j >= 0; j--) {           /* load the shuffle table */
      k = (*idum) / IQ1;
      *idum = IA1 * (*idum - k*IQ1) - k*IR1;
      if (*idum < 0)
        *idum += IM1;
      if (j < NTAB)
        iv[j] = *idum;
    }
    iy = iv[0];
  }

  k = (*idum) / IQ1;
  *idum = IA1 * (*idum - k*IQ1) - k*IR1;
  if (*idum < 0)
    *idum += IM1;
  k = idum2/IQ2;
  idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;
  if (idum2 < 0)
    idum2 += IM2;
  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if (iy < 1)
    iy += IMM1;
  if ((temp = AM * iy) > RNMX)
    return RNMX;                                /* avoid endpoint */
  else
    return temp;
}



double f(double u)
{
	double output;
	output = 1/(1+exp(-beta*(u-h)));
	return output;
	
}

double f_prime(double u)
{
	double output;
	output = 1/(1+exp(-beta*(u-h)))*1/(1+exp(-beta*(u-h)))*exp(-beta*(u-h))*beta;
	return output;
	
}

double f_doubleprime(double u)
{
	double output;
	output=0.2e1 * pow(0.1e1 + exp(-beta * (u - h)), -0.3e1) * beta * beta * pow(exp(-beta * (u - h)), 0.2e1) - pow(0.1e1 + exp(-beta * (u - h)), -0.2e1) * beta * beta * exp(-beta * (u - h));
	//output = (beta*beta* exp(beta*(h + u))*(exp(beta*h) - exp(beta*u)))/( (exp(beta*h) + exp(beta*u)) * (exp(beta*h) + exp(beta*u)) * (exp(beta*h) + exp(beta*u)) );
	return output;
	
}





