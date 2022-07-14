
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

#define T 1024*16//Time of simulation in msec
#define T_window 1000 //Time window for firing rate calculations
#define T_min 1000//minimal time post transient to initiate statistic calculations 
#define S 7//number of conditions - no conduction velocities, correlation time constants, spatial correlations 
#define Q 20//number of trials
#define N 100 //Number of neurons/populations
#define dt 0.1 //time step for integration
#define Dt 0.001// Value in sec of 1 dt e.g. Dt=0.001 means that 1 dt equals 1 msec
#define w_o 0.21//mean synaptic coupling
#define I_o 0.0// Input amplitude
#define temporal_window 10 //Temporal window for spike coherence measures
#define beta 25 //activation function non-linear gain
#define rho 0.15//Connexion probability
#define Omega 10.0 //Network spatial dimension in mm
#define h 0.10//activation function threshold (inflexion point) 
#define c_min 0.1//minimal conduction velocity
#define c_max 100//maximal conduction velocity

#define alpha_retraction 0.0001 //rate of myelin retraction 
#define alpha_formation 0.001//rate of myelin formation

#define t_stim_on 0//stimulus onset time
#define t_stim_off T//stimulus offset time

#define T_plastic_on T//plasticity onset time (post transient)
#define T_plastic_off T//plasticity offset time

//stimulated nodes
#define node_1 0
#define node_2 30

//receording nodes
#define node_3 70
#define node_4 100

double f(double u, double threshold);//activation function of neurons
double f_prime(double u, double threshold);//derivative of the activation function of neurons
double Heaviside(double input);//heaviside function
float ran2(long *idum);//random number generator - initial cond.
void shuffle_indices(int Array[], int size);//shuffle an array with int  entries
void shuffle(double Array[], int size);//shuffle an array with double  entries
void sort(double Array[], int size);//sort array 
void four1(double data[], unsigned long nn, int isign);//FFT
double coherence(double X1[], double X2[], int window, int total_time);//spike coherence calculations



double u_e[N][T];//membrane voltage of neurons
double mean_u_e[T];//mean network membrane woltage in time 
double mean_coherence[T];//mean spike coherence in time
double stability_metric[T];//mean stability metric in time
double mean_correlation[T]; //mean cond firing rate/correlations in time
double rate_e[T];//mean network firing rate in time
double rates[N][T];//Individual firing rates in time
double x[N];//spatial position of neurons in x 
double y[N];//spatial position of neurons in y 
double z[N];//spatial position of neurons in z 
double u_theo[T];//theoretical prediction of mean membrane voltage
double c[N][N];//conduction velocity matrix
double slope[N][N];//change in conduction velocity matrix
double w[N][N];//synaptic weight matrix
double l[N][N];//axonal lenghts matrix
int tau[N][N];//conduction delay matrix (int)
double delay[N][N];//condcution delay matrix (double)
int tau_o[N][N];//baseline conduction delays
double X_e[N][T];//noise array
double r_e[N][T];//instantaneous firing rates for display purposes
double Freq[T/2];//Array of output frequencies for PSD display
double PSD_network[T];//power spectrum network response of selected nodes
double PSD_signal[T];//power spectrum of stimulus signal
double MI_single_trial[Q];//Mutual information for a single trial
double MI[S][S];//mutual information as a function of correlation time and spatial correlation
double Var_MI[S][S];//variance of mutual information accross trials
double CORRELATION_TIME[S];//correlation times to be evaluated
double AMPLITUDE[S];//amplitudes to be evaluated
double SPATIAL_CORRELATION[S];//spatial correlations to be evaluated
double MEAN_MI_ACROSS_FREQUENCIES[S];//mean mutual information across frequencies
double 	MEAN_MI_ACROSS_FREQUENCIES_var[S];//variance mutual information across trials
double CV[S];//condcution velocities to be evaluated
double rate_stim[T];//mean firing rate of stimulated neurons
double rate_rec[T];//mean firing rate of receiving/recorded nodes
double Input[N][T];//Stimulus
double xi[T];//noise array for stimulus
double eta_c[T];//noise array for stimulus
double eta_i[N][T];//noise array for stimulus
double OU_process[N][T];//net resulting OU process stimulating the network

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
         	
         			
         		  			ufile>>c[i][j];
         		
			}
                         
      }  

					
		double sum=0;
		for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
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
						
         		
			}
                         
         }  
									
					
//conduction velocity to be evaluated		
CV[0]=0.1;
CV[1]=0.5;		
CV[2]=1;				
CV[3]=5;	
CV[4]=10;
CV[5]=50;		
CV[6]=100;				
//correlation times to be evaluated
CORRELATION_TIME[0]=  1;
CORRELATION_TIME[1] = 1;
CORRELATION_TIME[2] = 1;
CORRELATION_TIME[3] = 1;
CORRELATION_TIME[4] = 1;
CORRELATION_TIME[5] = 1;
CORRELATION_TIME[6] = 1;		
		
//spatial correlation to be evaluated		
SPATIAL_CORRELATION[0]=  0;
SPATIAL_CORRELATION[1] = 0.166;
SPATIAL_CORRELATION[2] = 0.333;
SPATIAL_CORRELATION[3] = 0.5;
SPATIAL_CORRELATION[4] = 0.666;
SPATIAL_CORRELATION[5] = 0.833;
SPATIAL_CORRELATION[6] = 1;			
					

			
				
			
for(int s1=0;s1<S;s1++)
{
			//CV[s1] =pow(10,s1-1);	;//c_min+100*s1/((double) S-1);
			
			MEAN_MI_ACROSS_FREQUENCIES[s1]=0;
			MEAN_MI_ACROSS_FREQUENCIES_var[s1]=0;
				
	for(int s2=0;s2<S;s2++)
	{
		
			AMPLITUDE[s2] = 0.1;//pow(10,s2-3);
			
		
			
				cout<<s1<<"	"<<s2<<endl;
				


				for (int q=0;q<Q;q++)
				{
					MI_single_trial[q]=0;
				}
				for(int q=0;q<Q;q++)
				{
					
					///////////////	reset connectivity	///////////////////////////////////////////////////////////////

						for (int i=0;i<N;i++)
						{
							for (int j=0;j<N;j++)
							{
								//c[i][j]=100;
								w[i][j]=0;
					    		
							}
						}
								double sum=0;
						for(int i=0;i<N;i++)
				         {
				         	for (int j=0;j<N;j++)
				         	{
				         			if(i==j)
											{
												w[i][j]=0;
											}
											else
											{
												if (fabs(i-j)<N)
												{
													long seed=(long)21+i*j*187198+56*j+s2+s1*q+12*1+18*q;
						                          	  
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
												//c[i][j]=c_min;
											c[i][j] = CV[s1];
											delay[i][j] = (double) (l[i][j]/(c[i][j]));
						         			tau[i][j] = (int) floor(l[i][j]/(c[i][j]));
								 			tau_o[i][j]=tau[i][j];
											//c[i][j]=c_min;
											
				         		
							}
				                         
				         }  
						for (int s=0;s<T-1;s++)
						{
						
								long d=rand()%105;      
			                    long noisyseed1=(long)21*s+s1*q*187198+56*d+12*q+14*s2;
			                    long noisyseed2=(long)69*s*s+121*s+q+45*d+2*1+q; 
	
								eta_c[s]= sqrt(-2*log(ran2(&noisyseed1)))*cos(2*3.14159*ran2(&noisyseed2));
								
								for (int k=0;k<N;k++)
								{
									 long noisyseed3=(long)88*k+s*q*67+56*d+12*k+d;
			                   		 long noisyseed4=(long)468*s*k+78*s+s1+322*d+q*221+k+s2; 
									 eta_i[k][s]=sqrt(-2*log(ran2(&noisyseed3)))*cos(2*3.14159*ran2(&noisyseed4));
								}
	
			            		for (int k=0;k<N;k++)
			            		{
		
				            		OU_process[k][s+1] = OU_process[k][s]+dt*(-(1/CORRELATION_TIME[s2])*OU_process[k][s])+sqrt(2*AMPLITUDE[s2]*dt)*(sqrt(SPATIAL_CORRELATION[s2])*eta_c[s]+sqrt(1-SPATIAL_CORRELATION[s2])*eta_i[k][s]);
									
								}
						
						}
						
			        	for (int t=0;t<T;t++)
						{
							
							double sum=0;
							for(int n=0;n<N;n++)
			        		{
			                            long d=rand()%105;      
			                            long noisyseed1=(long)21*t+n*s2*187198+34*s1+q+56*d+12*1+1;
			                            long noisyseed2=(long)69*t*n+11*s2+t+45*d+2*q+10*s1; 
			                          	long noisyseed3=(long)21*t+n*s2*666+34*s1+q+89*d+12*1+1;
			                            long noisyseed4=(long)34*t*n+235*s2+t+32*d+2*q+45*s1; 
			            				 
			            				u_e[n][t]=sqrt(-2*log(ran2(&noisyseed1)))*cos(2*3.14159*ran2(&noisyseed2));
			            				X_e[n][t]=0;
			            			
			            				if (t>t_stim_on&&t<t_stim_off)
										{
												if (n>node_1 && n<node_2)
												{
													Input[n][t] = 	OU_process[n][t]; 
												}
												else
												{
													Input[n][t]=0;
												}
									//I=I_o;
										}
							}
							mean_u_e[t]=0;
							rate_stim[t]=0;
							rate_rec[t]=0;
							PSD_signal[t]=0;
							PSD_network[t]=0;
						
			           }
///////////////////EULER SCHEME//////////////////////////////////////////////////////////////////////////				


					for (int t=0;t<T;t++)
					{
						
							
						double sum_ee; double sum_theo=0;double Gamma=0;
						for (int i=0;i<N;i++)
						{
						
							sum_ee=0;
						
							for (int j=0;j<N;j++)
							{								
									sum_ee=sum_ee+1/((double)N)*w[i][j]*X_e[j][t-tau[i][j]];
									
	 								
							}
							
						
							u_e[i][t+1]=u_e[i][t]+dt*(-1*u_e[i][t]+sum_ee+Input[i][t]);
					
										
										   // Poisson neurons
										 
										  long seed_q =rand()%101;  
						   				  long seed_e = (long) 12+i*s1+t+15+q*t*seed_q+56*s2+i+99*t+2*s1;              
			                              double p_e =ran2(&seed_e);
			                              double p_fire_e = ( 1-exp(-f(u_e[i][t+1],h)*dt)); 
			                              
			                              
			                              if (p_e<p_fire_e)//spike occurs
			                             
			                               {                   
			                                               X_e[i][t+1] =1/dt; 
			                                               
			        										
			                               }
			                               else//no spike
			                               {
			                                               X_e[i][t+1]=0; 
			                                               			                               
										   }
										   
										   
										   
              		  	
			                              
						}
								
					
			
				
						for (int i=0;i<N;i++)
						{
									mean_u_e[t+1] = mean_u_e[t+1]+1/((double)N)*u_e[i][t+1];
				
						}
					
					
									
						
					}// t loop
					

////////////Firing Rate////////////////////////////////////////////////
							int twindow=20;
								int range_stim= (int) node_2-node_1;
								int range_rec = (int) node_4-node_3 ;
								for (int t=twindow;t<T;t++)
								{
									double rate_stim_counter=0;
									double rate_rec_counter=0;
								
									double counts_stim=0;
									double counts_rec=0;
								
							
									for (int time=0;time<twindow;time++)
									{
										for (int index=node_1;index<node_2;index++)
										{
												if(X_e[index][t-(int)twindow/2+time]>0.1)
												{
													rate_stim_counter=rate_stim_counter+1;///((double) Ne);
													counts_stim = counts_stim+1;
												}
												else{}
										}
										for (int index=node_3;index<node_4;index++)
										{
												if(X_e[index][t-(int)twindow/2+time]>0.1)
												{
													rate_rec_counter=rate_rec_counter+1;///((double) Ne);
													counts_rec = counts_rec+1;
												}
												else{}
										}
									}
									
										rate_stim[t]=rate_stim_counter/((double) twindow*Dt)*1/((double)range_stim);
										rate_rec[t]=rate_rec_counter/((double) twindow*Dt)*1/((double)range_rec);
										
								}
									
									
			////////////////////////MUTUAL INFORMATION//////////////////////////////////////
						for (int k=0;k<(int)T;k++)
						{
												
												PSD_network[k] = 	rate_rec[k];
												PSD_signal[k]=		Input[node_1+1][k];													 		 
						}
									
						unsigned long nn=T/2;
						four1(PSD_network-1, nn,1);
																 
						for (int k=0;k<T/2;k++)
						{
											PSD_network[k] = 1/((double)T*T)*(fabs(PSD_network[k])*fabs(PSD_network[k])+fabs(PSD_network[(int)T-k])*fabs(PSD_network[(int)T-k]));
						}
						
						four1(PSD_signal-1, nn,1);
																 
						for (int k=0;k<T/2;k++)
						{
											PSD_signal[k] = 1/((double)T*T)*(fabs(PSD_signal[k])*fabs(PSD_signal[k])+fabs(PSD_signal[(int)T-k])*fabs(PSD_signal[(int)T-k]));
						}
						//normalize power spectra
						double norm1=0;double norm2=0;
						for (int k=0;k<T/2;k++)
						{
											norm1 = norm1 + PSD_network[k];
											norm2 = norm2 + PSD_signal[k];
			
						}
						for (int k=0;k<T/2;k++)
						{
										PSD_network[k] = PSD_network[k] /norm1;
										PSD_signal[k]=PSD_signal[k]/norm2;
			
						}
						double convo=0;
						double autoconv1 =0;
						double autoconv2 =0;
						for (int k=5;k<T/2;k++)
						{
											convo=convo+PSD_network[k]*PSD_signal[k];
											autoconv1 = autoconv1 + PSD_network[k]*PSD_network[k];
											autoconv2 = autoconv2 + PSD_signal[k]*PSD_signal[k];
										
						}
						
						//cout<<convo/sqrt(autoconv1*autoconv2)<<endl;;
				
						
						MI_single_trial[q]=0.5*log2(1+fabs(convo/sqrt(autoconv1*autoconv2))/(1-fabs(convo/sqrt(autoconv1*autoconv2))));
						
						}//q trial
						
						MI[s1][s2]=0;
						for (int q=0;q<Q;q++)
						{
							MI[s1][s2] = MI[s1][s2] +1/((double) Q)*MI_single_trial[q];
						}
						Var_MI[s1][s2]=0;
						for (int q=0;q<Q;q++)
						{
							Var_MI[s1][s2] = Var_MI[s1][s2] +1/((double) Q)*(MI[s1][s2]-MI_single_trial[q])*(MI[s1][s2]-MI_single_trial[q]);
						}
						
								MEAN_MI_ACROSS_FREQUENCIES[s1]=MEAN_MI_ACROSS_FREQUENCIES[s1]+1/((double)S)*MI[s1][s2];	
								MEAN_MI_ACROSS_FREQUENCIES_var[s1] = 	MEAN_MI_ACROSS_FREQUENCIES_var[s1] +1/((double)S)*Var_MI[s1][s2];
					}//s2 loop
					
					
}//s1 loop
							
		for(int k=0;k<(int) T/2;k++)
					
					{
								Freq[k] = k*1/((double) 2* (T)*Dt);
								//cout<<k*1/((double) 2* (T/2)*Dt)<<endl;
					}
		
/////////////OUTPUT//////////////////////////////////////////////////////////////////////////	
	
			
        ofstream outfile;
        
        
     	   outfile.open("ADM -  PSD.txt", ios::out);
         for(int k=5;k<(int) 50*(2*T)*Dt;k++)
         {
         	
         		  			outfile<<Freq[k]<<"	"<<PSD_network[k]<<"	"<<PSD_signal[k]<<endl;
         	
                         
         }  
        outfile.close(); 
     	 
	
        
          outfile.open("ADM  -  membrane potential dynamics.txt", ios::out);
         for(int t=0;t<T;t++)
         {
         	
         		  			outfile<<t<<"	"<<mean_u_e[t]<<"	"<<u_e[0][t]<<"	"<<u_e[2][t]<<"	"<<u_e[3][t]<<"	"<<u_e[N-1][t]<<"	"<<Input[node_1+1][t]<<endl;
         	
                         
         }  
        outfile.close(); 
        
       
        
        
     
          outfile.open("ADM  -  Mean activity and Firing rate.txt", ios::out);
         for(int t=0;t<T;t++)
         {
         	
         		  			outfile<<t<<"	"<<rate_stim[t]<<"	"<<rate_rec[t]<<"	"<<Input[node_1+1][t]<<endl;
         	
                         
         }  
        outfile.close(); 
        
                
        
    
    
        outfile.open("ADM  -  Connectivity.txt", ios::out);
         for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  			outfile<<i<<"	"<<j<<"	"<<c[i][j]<<"	"<<tau[i][j]<<"	"<<w[i][j]<<endl;
         		
			}
                         
         }  
        outfile.close(); 
        
  
        
          outfile.open("ADM  -  Mean mutual Information vs CV.txt", ios::out);
         for(int s1=0;s1<S;s1++)
         {
         	
			 
         	
         			
         		  			outfile<<CV[s1]<<"	"<<MEAN_MI_ACROSS_FREQUENCIES[s1]<<"	"<<MEAN_MI_ACROSS_FREQUENCIES_var[s1]<<endl;
         
			                         
         }  
          outfile.close(); 
        
          outfile.open("ADM  -  Mutual Information.txt", ios::out);
         for(int s1=0;s1<S;s1++)
         {
         	for(int s2=0;s2<S;s2++)
         	{
			 
         	
         			
         		  			outfile<<CV[s1]<<"	"<<CORRELATION_TIME[s2]<<"	"<<	SPATIAL_CORRELATION[s2]<<"	"<<MI[s1][s2]<<"	"<<Var_MI[s1][s2]<<endl;
         	}
			                         
         }  
        outfile.close(); 
        
         outfile.open("ADM  -  E Spikes.txt", ios::out);
         for(int t=0;t<T;t++)
         {
         	for (int n=0;n<N;n++)
         	{
         		if(X_e[n][t]>0.1)
         		{
         			
         		  			outfile<<t<<"	"<<n<<endl;
         		}
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

double f_prime(double u,double threshold)
{
	double output;
	output = 1/(1+exp(-beta*(u-threshold)))*1/(1+exp(-beta*(u-threshold)))*exp(-beta*(u-threshold))*beta;
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




