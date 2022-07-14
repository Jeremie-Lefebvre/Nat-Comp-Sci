
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


//R_A_N_2_P_A_R_A_M_E_T_E_R_S___________________________________________________

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


//_P_A_R_A_M_E_T_E_R_S__________________________________________________________

#define T 1000*100//Time of simulation in msec
#define T_window 1000 //sliding time window for simulations
#define Delta_T_rate 50//time wondow for firing rate
#define Sample 5 //Number of nodes selected for plotting purposes
#define N 100//Number of neurons/populations
#define Trials 20//number of trials
#define S 5//discretization for parameter values 
#define dt 0.1 //time step for integration
#define Dt 0.001// Value in sec of 1 dt e.g. Dt=0.001 means that 1 dt equals 1 msec
#define w_o 0.21//0.21
#define temporal_window 10 //Temporal window for spike coherence measures
#define beta 25 //activation function non-linear gain
#define rho 0.15//Connexion probability
#define Omega 10.0 //Network spatial dimension in mm
#define h 0.10//activation function threshold (inflexion point) 
#define c_min 0.1//minimal conduction velocity
#define c_max 100//maximal conduction velocity

#define alpha_retraction 0.0001 //rate of myelin retraction 
#define alpha_formation 0.001//rate of myelin formation

#define t_stim_on T/4//stimulus onset time
#define t_stim_off T//stimulus offset time

#define T_plastic_on 0//plasticity onset time (post transient)
#define T_plastic_off T//plasticity offset time

#define t_delay_statistics_initial 20//initial time in which delay statistics are computed
#define t_delay_statistics_final T-T_window-100//final time in which delay statistics are computed
#define delay_graining 50//number of bins for delay statistics
#define speed_graining 50//number of bins for conduction velocity statistics
#define max_delay 500//maximal delay for delay statistics
#define max_speed 40//maximum conduction velocity for conduction velocity statistics

double f(double u, double threshold);//activation function of neurons
double f_prime(double u, double threshold);//derivative of the activation function of neurons
double Heaviside(double input);//heaviside function
float ran2(long *idum);//random number generator - initial cond.
void shuffle_indices(int Array[], int size);//shuffle an array with int  entries
void shuffle(double Array[], int size);//shuffle an array with double  entries
void sort(double Array[], int size);//sort array 
void four1(double data[], unsigned long nn, int isign);//FFT
double coherence(double X1[], double X2[], int window, int total_time);//spike coherence calculations

double u_e[N][T_window];//membrane voltage of neurons
double mean_u_e[T];//mean network membrane woltage in time 
double mean_rate[T];//mean network firing rate in time
double x[N];//spatial position of neurons in x 
double y[N];//spatial position of neurons in y 
double z[N];//spatial position of neurons in z 
double c[N][N];//conduction velocity matrix
double c_o[N][N];//baseline conduction velocity matrix
double w[N][N];//synaptic weight matrix
double l[N][N];//axonal lenght matrix
double eta_c[T];//noise array for stimulus
double eta_i[N][T];//noise array for stimulus
double OU_process[N][T];//OU process stimulus
double Xi_e[N][T_window];//localized noisy input
int tau[N][N];//conduction delays
int tau_o[N][N]; //baseline conduction delays
double X_e[N][T_window];//individual spike trains in time
double r_e[N][T_window];//individual firing rates in time
double I[N];//stimulus (spatial)
double CV_in_time[Sample][T];//selected conduction velocities for display purposes
double Delay_in_time[Sample][T];//selected conduction delays for display purposes
double  Delay_distribution_final[delay_graining];//delay distribution after learning
double  Delay_distribution_initial[delay_graining];//delay distribution before learning
double Speed_distribution_final[speed_graining];//conduction velocity distribution after learning
double Speed_distribution_initial[speed_graining];//conduction velocity distribution before learning

double AMPLITUDE[S];//Tonic input amplitude
double NORMALIZED_MYELINATION[S][S];//normalized myelination
double 	RATE[S][S];//mean network firing rate

using namespace std;


int main()
{
	
srand (time(NULL));


			cout<<"Importing data..."<<endl;
		ifstream sfile("LENGTHS10.txt", ios::out);
         for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  			sfile>>l[i][j];
         		
			}
                         
         }  
        
        
        ifstream ufile("CVS10.txt", ios::out);
         for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  			ufile>>c_o[i][j];
        
			}
                         
         }  

// for(int i=0;i<N;i++)
//         {
//         	for (int j=0;j<N;j++)
//         	{
//         	
//         			
//         		  			c_o[i][j]=c_min;
//        
//			}
//                         
//         }  





		

//					for (int i=0;i<N;i++)
//					{
//						long seed1=(long)21+i*187198+56*i+12*1+1;
//						long seed2=(long)21+i*56+56*i+11*1+1;
//						long seed3=(long)21+i*789+56*i+8*1+1;
//		                          	  					
//						x[i] = ran2(&seed1)*Omega;
//						y[i] = ran2(&seed2)*Omega;
//						z[i] = ran2(&seed3)*Omega;
//
//					}
					double mean_w=0;
					for(int i=0;i<N;i++)
					{
						for (int j=0;j<N;j++)
						{
						//	l[i][j]= fabs(sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j])));
//							c[i][j]=c_min;//c_o;
							if(i==j)
							{
								w[i][j]=0;
							}
							else
							{
								if (fabs(i-j)<N)
								{
									long seed=(long)21+i*j*187198+56*j+12*1+1;
		                          	  
									if(ran2(&seed)<rho)
									{
											w[i][j]=w_o/rho;
											
									}
								}
								else
								{
									w[i][j]=0;
								}
								
							}
							mean_w=mean_w+1/((double)N*N)*w[i][j];
		 				
		 					//cout<<tau[i][j]<<endl;
						}
						
					}
					w[1][2]=w_o/rho;
					
					double c_mean_old=0;
					for (int i=0;i<N;i++)
					{
						for (int j=0;j<N;j++)
						{
							c_mean_old = c_mean_old+1/((double)N*N)*c_o[i][j];
						}
					}
				
					
					int sampled=0;
					int sampled_index_i[Sample];
					int sampled_index_j[Sample];
				
					for (int i=0;i<N;i++)
					{
							for (int j=0;j<N;j++)
							{
											
											if(w[i][j]>0&&sampled<Sample)
											{
												sampled_index_i[sampled] =i;
												sampled_index_j[sampled] =j;
												sampled++;
												
											}

						}				
					}
				
				
				
////////////////////EULER SCHEME//////////////////////////////////////////////////////////////////////////				
					
					//initial history
					
	for (int s1=0;s1<S;s1++)
	{
		
		for (int s2=0;s2<1;s2++)
		{
					AMPLITUDE[s1] = -0.1+0.2*s1/((double)S-1);
								
					cout<<AMPLITUDE[s1]<<endl;
				
					
					NORMALIZED_MYELINATION[s1][s2] = 0;
					RATE[s1][s2] = 0;
					for (int trial=0;trial<Trials;trial++)
					{
						cout<<s1<<"	"<<s2<<"	"<<trial<<endl;
						
						
						
					for (int i=0;i<N;i++)
					{
						for (int j=0;j<N;j++)
						{
							c[i][j] = c_o[i][j];
							tau[i][j] = (int) floor(l[i][j]/(c_o[i][j]));
		 					tau_o[i][j]=tau[i][j];	
						}	
						
					}	
					
					
					
					
		        	for (int s=0;s<T_window-1;s++)
					{
					
						for(int n=0;n<N;n++)
		        		{
		                            long d=rand()%105;      
		                            long noisyseed1=(long)21*s+n*n*187198+56*d+12*1+trial;
		                            long noisyseed2=(long)69*s*n+11*n+s+45*d+2*trial+1; 
		            				Xi_e[n][s]= sqrt(-2*log(ran2(&noisyseed1)))*cos(2*3.14159*ran2(&noisyseed2));
         							u_e[n][s]=Xi_e[n][s];	
         							r_e[n][s] = 0;
						}
						for(int i=0;i<Sample;i++)
						{
							CV_in_time[i][s]=c[sampled_index_i[i]][sampled_index_j[i]];
							Delay_in_time[i][s]=tau[sampled_index_i[i]][sampled_index_j[i]];
						}
						mean_rate[T_window+ s] =0;
						
						
					
		           }
					
					
					
					
					//Sliding window
					for (int s=0;s<T-T_window;s++)
					{
						//cout<<s<<endl;
						
						//update sliding window
						for (int t=0;t<T_window-1;t++)
						{
							for (int i=0;i<N;i++)
							{
									u_e[i][t] = u_e[i][t+1]	;
									X_e[i][t] = X_e[i][t+1];
									r_e[i][t] = r_e[i][t+1];
								
									Xi_e[i][t] = Xi_e[i][t+1]; 
							}	
						}
						
						
						//define input value
					
						if (s>t_stim_on&&s<t_stim_off)
						{
							for (int k=0;k<N;k++)
							{
							
								I[k] = AMPLITUDE[s1];
							}
						}
						else
						{
							for (int k=0;k<N;k++)
							{
								I[k]=0;
							}
						}
					
							
						double sum_ee; double sum_theo;
						for (int i=0;i<N;i++)
						{
						
							sum_ee=0;
							sum_theo=0;
							for (int j=0;j<N;j++)
							{								
									sum_ee=sum_ee+1/((double)N)*w[i][j]*X_e[j][T_window-1-tau[i][j]];
								
	 								
							}

						//	u_theo[T_window-1] = u_theo[T_window-1]+dt*alpha_u*(-1*u_theo[T_window-1]+sum_theo+I);
							u_e[i][T_window-1]=u_e[i][T_window-1]+dt*(-1*u_e[i][T_window-1]+sum_ee+I[i]);
					
										
										   // Poisson neurons
										 
										  long seed_q =rand()%101;  
						   				  long seed_e = (long) 12+i*i+s+15+s*seed_q+56*i+i+99*s+2;              
			                              double p_e =ran2(&seed_e);
			                              double p_fire_e = ( 1-exp(-f(u_e[i][T_window-1],h)*dt)); 
			                              if (p_e<p_fire_e)//spike occurs
			                               {                   
			                                               X_e[i][T_window-1] =1/dt; 
			                                               
			        										
			                               }
			                               else//no spike
			                               {
			                                               X_e[i][T_window-1]=0; 
			                                               			                               
										   }
										   double sum=0;
										   for (int t_index=0;t_index<Delta_T_rate;t_index++)
										   {
										   			sum=sum+1/((double)Delta_T_rate*Dt)*X_e[i][T_window-1-t_index]*dt;
										   			
										   }
										   r_e[i][T_window-1]=sum;

			                              
						}
						
						if(s>T_plastic_on&&s<T_plastic_off)
						{
								
						for (int i=0;i<N;i++)
						{
							for (int j=0;j<N;j++)
							{
								if (i==j)
								{
									c[i][j]=c_min;//
								}
								else
								{
									
										c[i][j]=c[i][j]+dt*((c_min-c[i][j])*alpha_retraction+0.3*X_e[j][T_window-1]*(l[i][j]/(c[i][j]))*alpha_formation);//kappa*(r-1));
									
								}
								 if(c[i][j]<c_min)
								 {
								 	c[i][j]=c_min;
								 }
								 else
								 {
								 	
								 }
								tau[i][j] = (int) floor(l[i][j]/(c[i][j]));
								
							}
						}
						}
						
						for(int i=0;i<Sample;i++)
						{
							CV_in_time[i][T_window+s]=c[sampled_index_i[i]][sampled_index_j[i]];
							Delay_in_time[i][T_window+s]=tau[sampled_index_i[i]][sampled_index_j[i]];
						}
						
						mean_rate[T_window+ s]=0;
						for (int i=0;i<N;i++)
						{
									mean_u_e[T_window+ s] = mean_u_e[T_window+ s]+1/((double)N)*u_e[i][T_window-1];
									mean_rate[T_window+ s] = mean_rate[T_window+ s]+1/((double)N)*r_e[i][T_window-1];
						}
						

					}// s loop
				

/////////////OUTPUT//////////////////////////////////////////////////////////////////////////	
	
		
			
			double ave=0;double var=0;
			for (int s=(int)T/2;s<T;s++)
			{
				ave = ave +1/((double)T/2)*mean_rate[s];
			}
			for (int s=(int)T/2;s<T;s++)
			{
				var = var +1/((double)T/2)*(mean_rate[s]-ave)*(mean_rate[s]-ave);
			}
		
			double c_mean_new=0;
		
			
			for (int i=0;i<N;i++)
			{
				for (int j=0;j<N;j++)
				{
					c_mean_new = c_mean_new+1/((double)N*N)*c[i][j];
				}
			}
			NORMALIZED_MYELINATION[s1][s2] = NORMALIZED_MYELINATION[s1][s2]+1/((double)Trials)*c_mean_new/c_mean_old;
			RATE[s1][s2] = RATE[s1][s2]+1/((double)Trials)*ave;
			
			}//trial loop
	}// s2 loop
}//s1 loop
			
			
			
			
        ofstream outfile;
        
        
     	  
     	 
        outfile.open("ADM - tonic inputs - Normalized myelination vs corr vs var.txt", ios::out);
         for(int s1=0;s1<S;s1++)
         {
         	for (int s2=0;s2<1;s2++)
         	{
			 
         	
         		  			outfile<<AMPLITUDE[s1]<<"	"<<NORMALIZED_MYELINATION[s1][s2]<<"	"<<RATE[s1][s2]<<endl;
         }
                         
         }  
        outfile.close(); 

        
        outfile.open("Hebbian myelination  - long runs -  Conduction velocities in time.txt", ios::out);
         for(int t=0;t<T-T_window;t=t+200)
         {
         	
         		  			outfile<<t<<"	"<<CV_in_time[0][t]<<"	"<<CV_in_time[1][t]<<"	"<<CV_in_time[2][t]<<"	"<<CV_in_time[3][t]<<"	"<<CV_in_time[4][t]<<endl;
         	
                         
         }  
        outfile.close(); 
        
         outfile.open("Hebbian myelination  - long runs -  Conduction delays in time.txt", ios::out);
         for(int t=0;t<T-T_window;t=t+200)
         {
         	
         		  			outfile<<t<<"	"<<Delay_in_time[0][t]<<"	"<<Delay_in_time[1][t]<<"	"<<Delay_in_time[2][t]<<"	"<<Delay_in_time[3][t]<<"	"<<Delay_in_time[4][t]<<endl;
         	
                         
         }  
        outfile.close(); 
        
        
        
          outfile.open("Hebbian myelination  - long runs -  Mean activity and Firing rate.txt", ios::out);
         for(int t=0;t<T-T_window;t=t+200)
         {
         	
         		  			outfile<<t<<"	"<<mean_u_e[t]<<"	"<<mean_rate[t]<<"	"<<endl;
         	
                         
         }  
        outfile.close(); 
        
        
        
      	outfile.open("Hebbian myelination  - long runs - Initial and Final Delay distribution.txt", ios::out);
		
		
		
				for (int k=0;k<delay_graining;k++)
				{
				 
				 	
		         					outfile<<(max_delay)/((double) delay_graining)*k<<"	"<<Delay_distribution_initial[k]<<"	"<<Delay_distribution_final[k]<<endl;
		         		 
		     	}
		     	
   
        outfile.close(); 
        
        	outfile.open("Hebbian myelination  - long runs -  Initial and Final Conduction Velocity distribution.txt", ios::out);
		
		
		
				for (int k=0;k<speed_graining;k++)
				{
				 
				 	
		         					outfile<<(max_speed)/((double) speed_graining)*k<<"	"<<Speed_distribution_initial[k]<<"	"<<Speed_distribution_final[k]<<endl;
		         		 
		     	}
		     	
   
        outfile.close(); 
        
       
        
       
        
    
        outfile.open("Hebbian myelination  -  Connectivity.txt", ios::out);
         for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  			outfile<<i<<"	"<<j<<"	"<<c[i][j]<<"	"<<tau[i][j]<<"	"<<w[i][j]<<endl;
         		
			}
                         
         }  
        outfile.close(); 
        
      outfile.open("LENGTHS.txt", ios::out);
         for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  			outfile<<l[i][j]<<endl;
         		
			}
                         
         }  
        outfile.close(); 
        
         outfile.open("CVS.txt", ios::out);
         for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  			outfile<<c[i][j]<<endl;
         		
			}
                         
         }  
        outfile.close(); 
      
        
        
 cout<<"Simulations complete..."<<endl;
      
return 0;    
}

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



double f(double u,double threshold)
{
	double output;
	output = 1/(1+exp(-beta*(u-threshold)));
	return output;
	
}

double Heaviside(double input)
{
	double output=0;
	if (input>0)
	{
		output=1;
	}
	else{}
	return output;
}


void sort(double Array[],int size)
{
	
	for (int q=0;q<size;q++)
	{
	
		for (int i=0;i<size;i++)
		{
				double temp;
				if (Array[i]>Array[i+1])
				{
				temp = Array[i+1];
				Array[i+1]= Array[i];
				Array[i] =temp;
				
				}
		
		}
	}
	
	
}



void shuffle_indices(int Array[], int size)
{
   int temporary;
   int randomNum;
   int last;
 
   for (last = size; last > 1; last=last-1)
   {
      randomNum = (int) floor(rand() % last);
      temporary = Array[randomNum];
      Array[randomNum] = Array[last - 1];
      Array[last - 1] = temporary;
   }
}

void shuffle(double Array[], int size)
{
   double temporary;
   int randomNum;
   int last;
 
   for (last = size; last > 1; last=last-1)
   {
      randomNum = rand() % last;
      temporary = Array[randomNum];
      Array[randomNum] = Array[last - 1];
      Array[last - 1] = temporary;
   }
}


/******************************************************************************/
void four1(double data[], unsigned long nn, int isign)
/*******************************************************************************
Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as
1; or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform,
if isign is input as -1.  data is a complex array of length nn or, equivalently,
a real array of length 2*nn.  nn MUST be an integer power of 2 (this is not
checked for!).
*******************************************************************************/
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) { /* This is the bit-reversal section of the routine. */
		if (j > i) {
			SWAP(data[j],data[i]); /* Exchange the two complex numbers. */
			SWAP(data[j+1],data[i+1]);
		}
		m=nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}

	mmax=2;
	while (n > mmax) { /* Outer loop executed log2 nn times. */
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax); /* Initialize the trigonometric recurrence. */
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) { /* Here are the two nested inner loops. */
			for (i=m;i<=n;i+=istep) {
				j=i+mmax; /* This is the Danielson-Lanczos formula. */
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr; /* Trigonometric recurrence. */
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
   
double coherence(double X1[], double X2[], int window)
{
		int Bins = T/((int) window);
		double Y1[Bins];
		double Y2[Bins];
		//create binned spike train representations
		for (int k=0;k<Bins;k++)
		{
			Y1[k]=0;
			Y2[k]=0;
			for (int s=0;s<window;s++)
			{
				if (X1[k*window+s]>0)
				{
					Y1[k]=1;
				}
				if (X2[k*window+s]>0)
				{
					Y2[k]=1;
				}
				
			}
					
		}
		//compute correlation
		double sumY1Y2=0;
		double sumY1=0;
		double sumY2=0;
		for (int k=0;k<Bins;k++)
		{
			sumY1Y2 = sumY1Y2+Y1[k]*Y2[k];
			sumY1 = sumY1+Y1[k];
			sumY2 = sumY2+Y2[k];
		}
		
		double output;
		if (sumY1>0&&sumY2>0)
		{
				output = sumY1Y2/sqrt(sumY1*sumY2);
		}
		else{ output =0;}		 
		return output;
		
		

}




