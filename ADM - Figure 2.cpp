
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
#define T_min 1000//minimal time post transient to initiate statistic calculations 
#define T 1000*100//Time of simulation in msec
#define T_window 1000 //Time window for firing rate calculations
#define Sample 10 //Number of nodes selected for plotting purposes
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

#define t_stim_on T//stimulus onset time
#define t_stim_off T//stimulus offset time

#define T_plastic_on 20000//plasticity onset time (post transient)
#define T_plastic_off T//plasticity offset time

#define t_delay_statistics_initial T_min+100//initial time in which delay statistics are computed
#define t_delay_statistics_final T-1000//final time in which delay statistics are computed
#define delay_graining 50//number of bins for delay statistics
#define speed_graining 50//number of bins for conduction velocity statistics
#define max_delay 200//maximal delay for delay statistics
#define max_speed 30//maximum conduction velocity for conduction velocity statistics

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
double CV_in_time[Sample][T];//selected conduction velocities for display purposes
double Delay_in_time[Sample][T];//selected conduction delays for display purposes
double  Delay_distribution_final[delay_graining];//delay distribution after learning
double  Delay_distribution_initial[delay_graining];//delay distribution before learning
double Speed_distribution_final[speed_graining];//conduction velocity distribution after learning
double Speed_distribution_initial[speed_graining];//conduction velocity distribution before learning

double temporal_mean_coherence[T];	//mean spike coherence in time
double temporal_ave_rate[T];//time averaged mean firing rate in time
double temporal_stability_metric[T];//stability metric in time

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
         
        
          for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  		c[i][j]=c_min;
         		
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
						sum=sum+1/((double)N*N)*w[i][j];	
					delay[i][j] = (double) (l[i][j]/(c[i][j]));
         			tau[i][j] = (int) floor(l[i][j]/(c[i][j]));
		 			tau_o[i][j]=tau[i][j];
         		
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
				
							
					
					//define history
		        	for (int t=0;t<T;t++)
					{
						
		
						for(int n=0;n<N;n++)
		        		{
		                            long d=rand()%105;      
		                            long noisyseed1=(long)21*t+n*n*187198+56*d+12*1+1;
		                            long noisyseed2=(long)69*t*n+11*n+t+45*d+2*1+1; 
		                          
		            				X_e[n][t]=0; 
		            				u_e[n][t]=0;
		            				
						}
					
		           }
		               
				
		
		
		 		
						
					

			
				
			
						
				
////////////////////EULER SCHEME//////////////////////////////////////////////////////////////////////////				
					
					//Define integration order
					
					for (int t=T_min;t<T;t++)
					{
						cout<<t<<endl;
							
						double sum_ee; double sum_theo=0;double Gamma=0;
						for (int i=0;i<N;i++)
						{
						
							sum_ee=0;
						
							for (int j=0;j<N;j++)
							{								
									sum_ee=sum_ee+1/((double)N)*w[i][j]*X_e[j][t-tau[i][j]];
						
	 								
							}
							double I=0;
							if (t>t_stim_on&&t<t_stim_off)
							{
								I=I_o;
							}

						
							u_e[i][t+1]=u_e[i][t]+dt*(-1*u_e[i][t]+sum_ee+I);
					
										
										   // Poisson neurons
										 
										  long seed_q =rand()%101;  
						   				  long seed_e = (long) 12+i*i+t+15+t*seed_q+56*i+i+99*t+2;              
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
										   for (int t_index=0;t_index<T_window;t_index++)
										   {
										   			r_e[i][t]=r_e[i][t]+1/((double)T_window*Dt)*X_e[i][t+1-t_index]*dt;
										   			
										   }
										   
										   
              		  	
			                              
						}
						
						
						if (t>T_plastic_on&&t<T_plastic_off)
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
									c[i][j]=c[i][j]+dt*((c_min-c[i][j])*alpha_retraction+0.3*X_e[j][t-1]*(l[i][j]/(c[i][j]))*alpha_formation);
									slope[i][j] = 	((c_min-c[i][j])*alpha_retraction+0.3*X_e[j][t-1]*(l[i][j]/(c[i][j]))*alpha_formation);
								}
								 if(c[i][j]<0)
								 {
								 	c[i][j]=c_min;
								 }
								 else
								 {
								 	
								 }
								 		delay[i][j] = (double) l[i][j]/(c[i][j]);
										tau[i][j] = (int) floor(l[i][j]/(c[i][j]));
										
							}
						}
						}
						for(int i=0;i<Sample;i++)
						{
							CV_in_time[i][t]=c[sampled_index_i[i]][sampled_index_j[i]];
							Delay_in_time[i][t]=delay[sampled_index_i[i]][sampled_index_j[i]];
						}
						
					
			
				
						for (int i=0;i<N;i++)
						{
									mean_u_e[t+1] = mean_u_e[t+1]+1/((double)N)*u_e[i][t+1];
				
						}
						
						for (int i=0;i<N;i++)
						{
							for (int j=0;j<N;j++)
							{
							
									stability_metric[t] = stability_metric[t]+1/((double)N*N)*(slope[i][j]);
							}
						}
						
										if (t==t_delay_statistics_initial)
						{
						 cout<<t<<endl;
													for (int s=0;s<delay_graining;s++)
													{
																			Delay_distribution_initial[s]=0;
													}
																					
													for (int tt=0;tt<delay_graining;tt++)
													{
																			double ttau = 0 + (max_delay)/((double) delay_graining)*tt;
																			double dttau = (max_delay)/((double) delay_graining);
																							
																			for (int i=0;i<N;i++)
																			{
																				for (int j=0;j<N;j++)
																				{
																					
																									double ttau = 0 + (max_delay)/((double) delay_graining)*tt;
																									double dttau = (max_delay)/((double) delay_graining);
																									double delay  = (double) tau[i][j];			
																									if (((double)delay>ttau-dttau*0.5)&&((double)delay<ttau+dttau*0.5))
																									{
																												Delay_distribution_initial[tt] = Delay_distribution_initial[tt]+w[i][j];
													
																									}	
																				}
																			}
													}
									
														
										for (int s=0;s<speed_graining;s++)
										{
													Speed_distribution_initial[s]=0;
										}
										for (int tt=0;tt<speed_graining;tt++)
										{
																double tc = 0 + (max_speed)/((double) speed_graining)*tt;
																double dcc = (max_speed)/((double) speed_graining);
																
																for (int i=0;i<N;i++)
																{
													
																		for (int j=0;j<N;j++)
																		{
																			double tc = 0 + (max_speed)/((double) speed_graining)*tt;
																			double dcc = (max_speed)/((double) speed_graining);
																			double speed = (double) c[i][j];			
																			if (((double)speed>tc-dcc*0.5)&&((double)speed<tc+dcc*0.5))
																			{
																						Speed_distribution_initial[tt] = Speed_distribution_initial[tt]+w[i][j];
																						
																						
																			}	
																		}
															
																}
																
										}
										
						}//initial distribution loop
								
						if (t==t_delay_statistics_final)
						{
								cout<<t<<endl;
													for (int s=0;s<delay_graining;s++)
													{
																			Delay_distribution_final[s]=0;
													}
																					
													for (int tt=0;tt<delay_graining;tt++)
													{
																			double ttau = 0 + (max_delay)/((double) delay_graining)*tt;
																			double dttau = (max_delay)/((double) delay_graining);
																							
																			for (int i=0;i<N;i++)
																			{
																				for (int j=0;j<N;j++)
																				{
																					
																									double ttau = 0 + (max_delay)/((double) delay_graining)*tt;
																									double dttau = (max_delay)/((double) delay_graining);
																									double delay  = (double) tau[i][j];			
																									if (((double)delay>ttau-dttau*0.5)&&((double)delay<ttau+dttau*0.5))
																									{
																												Delay_distribution_final[tt] = Delay_distribution_final[tt]+w[i][j];
													
																									}	
																				}
																			}
													}
									
														
										for (int s=0;s<speed_graining;s++)
										{
													Speed_distribution_final[s]=0;
										}
										for (int tt=0;tt<speed_graining;tt++)
										{
																double tc = 0 + (max_speed)/((double) speed_graining)*tt;
																double dcc = (max_speed)/((double) speed_graining);
																
																for (int i=0;i<N;i++)
																{
													
																		for (int j=0;j<N;j++)
																		{
																			double tc = 0 + (max_speed)/((double) speed_graining)*tt;
																			double dcc = (max_speed)/((double) speed_graining);
																			double speed = (double) c[i][j];			
																			if (((double)speed>tc-dcc*0.5)&&((double)speed<tc+dcc*0.5))
																			{
																						Speed_distribution_final[tt] = Speed_distribution_final[tt]+w[i][j];
																						
																						
																			}	
																		}
															
																}
																
										}
										
								}// final distribution loop
						
						
						
					
									
						
					}// t loop
					
					
					
////////////Firing Rate////////////////////////////////////////////////
							int twindow=1000;
								int range_e=N;
							
								for (int t=twindow;t<T-twindow;t++)
								{
									double rate_e_counter=0;
								
									double countse=0;
								
							
									for (int time=0;time<twindow;time++)
									{
										for (int index=0;index<N;index++)
										{
										if(X_e[index][t-(int)twindow/2+time]>0.1)
										{
											rate_e_counter=rate_e_counter+1;///((double) Ne);
											countse = countse+1;
										}
										else{}
										}
									}
									
										rate_e[t]=rate_e_counter/((double) twindow*Dt)*1/((double)range_e);
									}
									
////////////mean firing rate correlation/////////////////////////////////////////////


int t_window=100;//time window for correlation calculations


for (int i=0;i<N;i++)
{
	
	for (int t=t_window;t<T-t_window;t++)
	{
		int rate_counter=0;
		for (int s=0;s<t_window;s++)
		{
			if(X_e[i][t-(int)twindow/2+s]>0.1)
			{
											rate_counter=rate_counter+1;
											
			}
			else{}
			
		}
		rates[i][t] = rate_counter/((double) twindow*Dt);
		
		
	}
	
	
	
}


int no_of_intervals  = 10;//sumber of intervals in which correlation and spatibility metrics will be computed

int t_interval =  (int) (T)/((double)no_of_intervals);   
double mean_correlation[no_of_intervals];
double mean_rate[no_of_intervals];
double mean_stability_metric[no_of_intervals];
double mean_conditional_firing_rate[no_of_intervals];
double times[no_of_intervals];                                          
int Pairings = N*N;

for (int k=0;k<no_of_intervals;k++)
{
							mean_rate[k]=0;
							mean_correlation[k]=0;
							mean_conditional_firing_rate[k]=0;
							for (int p=0;p<Pairings;p++)
							{
									long seed1 = 2342*p+56*5+3;
									long seed2 = 45*p+556*1+45;	
									int x1 = (int) (ran2(&seed1)*N);
									int x2 = (int) (ran2(&seed2)*N);
									
									
									
									if (x1==x2)
									{
									}
									else
									{
									
									double spike_cov12=0; 
									double mean_1=0;double mean_2=0;double cov_12=0;double var_1=0;double var_2=0;
									for (int s=0;s<t_interval;s++)
									{
											spike_cov12 = spike_cov12+1/((double)t_interval)*(X_e[x1][(int) 2*0+k*t_interval+s])*(X_e[x2][(int) 2*0+ k*t_interval+s]);
									}
									for (int s=0;s<t_interval;s++)
									{
										mean_1 = mean_1 +1/((double)t_interval)*rates[x1][2*0+(int) k*t_interval+s];
										mean_2 = mean_2 +1/((double)t_interval)*rates[x2][2*0+(int) k*t_interval+s];
										
									}
									for (int s=0;s<t_interval;s++)
									{
										var_1 = var_1 +1/((double)t_interval)*(rates[x1][(int) 2*0+k*t_interval+s]-mean_1)*(rates[x1][(int) 2*0+k*t_interval+s]-mean_1);
										var_2 = var_2 +1/((double)t_interval)*(rates[x2][(int) 2*0+k*t_interval+s]-mean_2)*(rates[x2][(int) 2*0+k*t_interval+s]-mean_2);
										cov_12 = cov_12 +1/((double)t_interval)*(rates[x1][(int) 2*0+k*t_interval+s]-mean_1)*(rates[x2][(int) 2*0+ k*t_interval+s]-mean_2);
									
									}
									mean_correlation[k] = mean_correlation[k] +1/((double) Pairings)*cov_12/(sqrt(var_1*var_2));
									mean_rate[k] = mean_rate[k] +10/((double) Pairings)*(0.5*mean_1+0.5*mean_2);
									mean_conditional_firing_rate[k] = mean_conditional_firing_rate[k]+1/((double) Pairings)*spike_cov12/sqrt(mean_1*mean_2);
								
									}
							}
							mean_stability_metric[k] = 0;
							for (int s=0;s<t_interval;s++)
							{
											mean_stability_metric[k] = mean_stability_metric[k]+1/((double) t_interval)*stability_metric[2*0+(int)k*t_interval+s] ;	
							}
							
											
							
							times[k] = (double)2*0+ k*t_interval;
							
}
								




				
				
//computing firing rate variances before and after plasticity		
						
		double sumr_baseline=0;double var_baseline=0;
		for (int t=T_min;t<T_min+10000;t++)
		{
			sumr_baseline=sumr_baseline+1/((double)10000)*rate_e[t];
		}
		for (int t=T_min;t<T_min+10000;t++)
		{
			var_baseline=var_baseline+1/((double)10000)*(rate_e[t]-sumr_baseline)*(rate_e[t]-sumr_baseline);
		}
		cout<<"	"<<sumr_baseline<<"	"<<sqrt(var_baseline)<<endl;
		
		double sum_stabilized=0;double var_stabilized=0;
		for (int t=T-10000;t<T;t++)
		{
			sum_stabilized=sum_stabilized+1/((double)10000)*rate_e[t];
		}
		for (int t=T-10000;t<T;t++)
		{
			var_stabilized=var_stabilized+1/((double)10000)*(rate_e[t]-sum_stabilized)*(rate_e[t]-sum_stabilized);
		}
		cout<<"	"<<sum_stabilized<<"	"<<sqrt(var_stabilized)<<endl;
					
		
/////////////OUTPUT//////////////////////////////////////////////////////////////////////////	
	
			
        ofstream outfile;
        
	
        
          outfile.open("ADM  -  sample membrane potential dynamics.txt", ios::out);
         for(int t=0;t<T;t++)
         {
         	
         		  			outfile<<t<<"	"<<mean_u_e[t]<<"	"<<u_e[0][t]<<"	"<<u_e[2][t]<<"	"<<u_e[3][t]<<"	"<<u_e[N-1][t]<<endl;
         	
                         
         }  
        outfile.close(); 
        
        outfile.open("ADM  -  Conduction velocities in time.txt", ios::out);
         for(int t=0;t<T;t++)
         {
         	
         		  			outfile<<t<<"	"<<CV_in_time[0][t]<<"	"<<CV_in_time[1][t]<<"	"<<CV_in_time[2][t]<<"	"<<CV_in_time[3][t]<<"	"<<CV_in_time[4][t]<<endl;
         	
                         
         }  
        outfile.close(); 
        
         outfile.open("ADM  -  Conduction delays in time.txt", ios::out);
         for(int t=0;t<T;t++)
         {
         	
         		  			outfile<<t<<"	"<<Delay_in_time[0][t]<<"	"<<Delay_in_time[1][t]<<"	"<<Delay_in_time[2][t]<<"	"<<Delay_in_time[3][t]<<"	"<<Delay_in_time[4][t]<<endl;
         	
                         
         }  
        outfile.close(); 
        
        
        outfile.open("ADM  -  Single connection activity and change in CV.txt", ios::out);
         for(int t=0;t<T;t++)
         {
         	
         		  			outfile<<t<<"	"<<u_e[sampled_index_j[1]][t]<<"	"<<X_e[sampled_index_j[1]][t]<<"	"<<r_e[sampled_index_j[1]][t]<<"	"<<CV_in_time[1][t]<<"	"<<Delay_in_time[1][t]<<endl;
         	
                         
         }  
        outfile.close(); 
        
          outfile.open("ADM  -  Mean corelations, conditional firing rate and stability metric in time.txt", ios::out);
         for(int k=0;k<no_of_intervals;k++)
         {
         	
         		  			outfile<<times[k]<<"	"<<mean_rate[k]<<"	"<<mean_correlation[k]<<"	"<<mean_stability_metric[k]<<"	"<<mean_conditional_firing_rate[k]<<endl;
         	
                         
         }  
        outfile.close(); 
        
         outfile.open("ADM  -  Mean activity and Firing rate.txt", ios::out);
         for(int t=0;t<T;t++)
         {
         	
         		  			outfile<<t<<"	"<<mean_u_e[t]<<"	"<<rate_e[t]<<"	"<<u_e[0][t]<<"	"<<f(mean_u_e[t],h)<<endl;
         	
                         
         }  
        outfile.close(); 
        
        
        
          outfile.open("ADM  -  Spikes.txt", ios::out);
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
        
    
        outfile.open("ADM  -  Connectivity.txt", ios::out);
         for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  			outfile<<i<<"	"<<j<<"	"<<c[i][j]<<"	"<<tau[i][j]<<"	"<<w[i][j]<<endl;
         		
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
   
double coherence(double X1[], double X2[], int window, int total_time)
{
		int Bins = total_time/((int) window);
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




