
//#1 DECLARATION SEQUENCE------------------------------------------------------------------------------------------------------------------------------------------

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

#define T 1000*140//Time of simulation in msec
#define T_window 1000
#define T_min 1000
#define Sample 10
#define N 100
#define dt 0.1 //time step for integration
#define Dt 0.001// Value in sec of 1 dt e.g. Dt=0.001 means that 1 dt equals 1 msec

#define I_o -0.0
#define temporal_window 5 //Temporal window for spike coherence measures
#define beta 25 //300//activation function non-linear gain
#define rho 0.15
#define Omega 1.0 //(mm)
#define h 0.16//activation function threshold (inflexion point) 
#define D 0.00000//0.015// ADD NOISE TO SUPPRESS SEIZURE Noise variance
#define c_min 0.1
#define c_max 100

#define alpha_retraction 0.0001
#define alpha_formation 0.001

#define t_stim_on T
#define t_stim_off T

#define T_plastic_on 40000
#define T_plastic_off T

#define t_delay_statistics_initial 1000
#define t_delay_statistics_final T-100
#define delay_graining 50
#define speed_graining 50
#define max_delay 200
#define max_speed 30

double f(double u, double threshold);
double f_prime(double u, double threshold);
double Heaviside(double input);
float ran2(long *idum);//random number generator - initial cond.
void shuffle_indices(int Array[], int size);//shuffle an array with int  entries
void shuffle(double Array[], int size);//shuffle an array with double  entries
void sort(double Array[], int size);//sort array 
void four1(double data[], unsigned long nn, int isign);
double coherence(double X1[], double X2[], int window);



double u_e[N][T];//Dynamics of E units
double mean_u_e[T];//mean network activity E units

double mean_coherence[T];
double stability_metric[T];
double mean_correlation[T]; 
double rate_e[T];
double rates[N][T];

double x[N];
double y[N];
double z[N];
double u_theo[T];
double c[N][N];
double slope[N][N];
double w[N][N];
double l[N][N];

int tau[N][N];
double delay[N][N];
int tau_o[N][N];
double X_e[N][T];
double r_e[N][T];
double CV_in_time[Sample][T];
double Delay_in_time[Sample][T];
double  Delay_distribution_final[delay_graining];
double  Delay_distribution_initial[delay_graining];
double Speed_distribution_final[speed_graining];
double Speed_distribution_initial[speed_graining];
using namespace std;


int main()
{
	
srand (time(NULL));




		cout<<"Importing data..."<<endl;
		ifstream sfile("LENGTHS1.txt", ios::out);
         for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  			sfile>>l[i][j];
         		
			}
                         
         }  
        
//        
//        ifstream ufile("CVS10.txt", ios::out);
//         for(int i=0;i<N;i++)
//         {
//         	for (int j=0;j<N;j++)
//         	{
//         	
//         			
//         		  			ufile>>c[i][j];
//         		
//			}
//                         
//         }  
         
          for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  		c[i][j]=c_min;
         		
			}
                         
         }  
         
         
     /*
	

					for (int i=0;i<N;i++)
					{
						long seed1=(long)21+i*187198+56*i+12*1+1;
						long seed2=(long)21+i*56+56*i+11*1+1;
						long seed3=(long)21+i*789+56*i+8*1+1;
		                          	  					
						x[i] = ran2(&seed1)*Omega;
						y[i] = ran2(&seed2)*Omega;
						z[i] = ran2(&seed3)*Omega;

					}
	
					double mean_w=0;
					for(int i=0;i<N;i++)
					{
						for (int j=0;j<N;j++)
						{
							l[i][j]= fabs(sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j])));
							c[i][j]=c_min;//c_o;
							
						}
						
					}
				*/	
					
		double Me[N][N];
		double Mi[N][N];
		double mu_e=7;
		double var_w=0.8;//
		double evsi=0.8;
		double List[N];
		for (int i=0;i<N;i++)
		{
			for (int j=0;j<N;j++)
			{
				long d=rand()%105;      
				long noisyseed1=(long)21*j+j*j*187198+56*d+12*1+1;
				long noisyseed2=(long)69*i*j+11*i+j+45*d+2*1+1; 
				long noisyseed3=(long)99*j+j*j*677+56*d+12*1+1;
				long noisyseed4=(long)77*i*j+55*i+j+6678*d+34*1+1; 
				List[j]=0;
				Me[i][j]= mu_e + sqrt(var_w)*sqrt(-2*log(ran2(&noisyseed1)))*cos(2*3.14159*ran2(&noisyseed2));
				Mi[i][j]= -(evsi)/(1-evsi)*mu_e + sqrt(var_w)*sqrt(-2*log(ran2(&noisyseed3)))*cos(2*3.14159*ran2(&noisyseed4));
			}
		}
		
		
		
			double sum=0;
		for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         			long seed=(long)21+i*j*187198+56*j+12*1+1;
		            if(ran2(&seed)<rho)	
					{
						if(j<(int)(evsi*N))
						{
							List[j] = Me[i][j];
						}
						else
						{
							List[j] = Mi[i][j];
						}
						
					}
					else
					{
						List[j]=0	;
					}	
					
							
							
						
					delay[i][j] = (double) (l[i][j]/(c[i][j]));
         			tau[i][j] = (int) floor(l[i][j]/(c[i][j]));
		 			tau_o[i][j]=tau[i][j];
		 	}
		 	shuffle(List,N);
		 	for (int j=0;j<N;j++)
		 	{
		 		
		 		if (i==j)
		 		{
		 			w[i][j]=0;
				 }
				 else
				 {
				 	w[i][j] = List[j];
				 }
				  
			 }
			
         		
		}
		double sum_row[N];
		for (int i=0;i<N;i++)
		{
			sum_row[i]=0;
			for (int j=0;j<N;j++)
			{
				sum_row[i] = sum_row[i]+w[i][j];
			}
			int number_of_nonzero_entries=0;
			for (int j=0;j<N;j++)
			{
				if (fabs(w[i][j])>0)
				{
					number_of_nonzero_entries=number_of_nonzero_entries+1;
				}
			}
			for (int j=0;j<N;j++)
			{
					if (fabs(w[i][j])>0)
					{
						w[i][j]= w[i][j] - sum_row[i]/((double)number_of_nonzero_entries);
					}
					sum=sum+w[i][j];
			}
			
			
		}
//		cout<<f(0,h)<<endl;
//		cout<<f_prime(0,h)<<endl;
//		cout<<sqrt(N*rho*(evsi*var_w+(1-evsi)*var_w+evsi/(1-evsi)*mu_e*mu_e)*f_prime(0,h)*f_prime(0,h))<<endl;
//		cout<<sum<<endl;
//		getch();
	
                         
          
          									
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
				
							
					
					// define noise arrays
		        	for (int t=0;t<T;t++)
					{
						
						double sum=0;
						for(int n=0;n<N;n++)
		        		{
		                            long d=rand()%105;      
		                            long noisyseed1=(long)21*t+n*n*187198+56*d+12*1+1;
		                            long noisyseed2=(long)69*t*n+11*n+t+45*d+2*1+1; 
		                          
		            				 
		            				u_e[n][t]=sqrt(-2*log(ran2(&noisyseed1)))*cos(2*3.14159*ran2(&noisyseed2));
		            				
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
									sum_theo  =sum_theo + 1/((double)N*N)*f(u_theo[t-tau[i][j]],h);
									Gamma = Gamma +1/((double)N*N)*exp(-delay[i][j]*dt);
	 								
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
//						 long noisyseed1=(long)21*t+t*t*187198+56*t+12*1+1;
//		                 long noisyseed2=(long)69*t*t+11*t+t+45*t+2*1+1; 
//						double xi = sqrt(-2*log(ran2(&noisyseed1)))*cos(2*3.14159*ran2(&noisyseed2));
////						u_theo[t+1] = u_theo[t]+dt*(-1*u_theo[t]+w_o*sum_theo)+sqrt(2*w_o*w_o*f(0.2,h)/2.0*1/((double)N)*dt)*xi;;
					
//						v[t+1] = v[t] + dt*(-v[t])+ sqrt(2*w_o*w_o*f(fixed_point,h)/2.0*1/(1-f_prime(fixed_point,h)*w_o*Gamma)*1/((double)N)*dt)*xi;
//						//v2[t+1] = v2[t]+dt*(-v2[t]+ w_o*X_e[0][t]-w_o*f(u_e[0][t],h));	
						
						
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
									double c_now = c[i][j];
									//c[i][j]=c[i][j]+dt*alpha_c*Heaviside(w[i][j])*(0.5*exp(0.03*r_e[j][t-1])-c[i][j]);
									//c[i][j]=c[i][j]+dt*(-(c[i][j])*alpha_retraction+0.1*X_e[j][t-1]*(1-c[i][j]/c_max)*alpha_formation);//kappa*(r-1));
									c[i][j]=c[i][j]+dt*((c_min-c[i][j])*alpha_retraction+2*X_e[j][t-1]*(l[i][j]/(c[i][j]))*alpha_formation);//kappa*(r-1));
									slope[i][j] = 	((c_min-c[i][j])*alpha_retraction+2*X_e[j][t-1]*(l[i][j]/(c[i][j]))*alpha_formation);//kappa*(r-1));
										
								}
								 if(c[i][j]<c_min)
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
						//coincidence_theo[t]=1/dt*1/dt*( 1-exp(-f(mean_u_e[t+1],h)*dt))*( 1-exp(-f(mean_u_e[t+1],h)*dt));//
						
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
							int twindow=20;
								int range_e=N;
							
								for (int t=twindow;t<T-twindow;t++)
								{
									double rate_e_counter=0;
								
								
							
									for (int time=0;time<twindow;time++)
									{
										for (int index=0;index<N;index++)
										{
										if(X_e[index][t-(int)twindow/2+time]>0.1)
										{
											rate_e_counter=rate_e_counter+1;///((double) Ne);
										
										}
										else{}
										}
									}
									
										rate_e[t]=rate_e_counter/((double) twindow*Dt)*1/((double)range_e);
									}
									
///////////////////////CORREALTIONS///////////////////////////////////////

int t_window=100;


for (int i=0;i<N;i++)
{
	
	for (int t=t_window;t<T-t_window;t++)
	{
		int rate_counter=0;
		for (int s=0;s<t_window;s++)
		{
			if(X_e[i][t-(int)twindow/2+s]>0.1)
			{
											rate_counter=rate_counter+1;///((double) Ne);
											
			}
			else{}
			
		}
		rates[i][t] = 0.1*rate_counter/((double) twindow*Dt);
		
		
	}
	
	
	
}

int no_of_intervals  = 14;
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
									mean_rate[k] = mean_rate[k] +1/((double) Pairings)*(0.5*mean_1+0.5*mean_2)*0.1;
									mean_conditional_firing_rate[k] = mean_conditional_firing_rate[k]+1/((double) Pairings)*spike_cov12/sqrt(mean_1*mean_2);
									//mean_conditional_firing_rate_theo[k] = mean_conditional_firing_rate_theo[k]+1/((double) Pairings)*spike_cov12/sqrt(mean_1*mean_2);

									}
							}
							mean_stability_metric[k] = 0;
							for (int s=0;s<t_interval;s++)
							{
											mean_stability_metric[k] = mean_stability_metric[k]+1/((double) t_interval)*stability_metric[2*0+(int)k*t_interval+s] ;	
							}
							
											
							
							times[k] = (double)2*0+ k*t_interval;//+ (double)t_interval/2;
							
}									
								
								
							
								
////////////Spike coherence/////////////////////////////////////////////
//							int Pairings = N;
//							double X1[T];
//							double X2[T];
//							double mean_corr =0;
//							for (int p=0;p<Pairings;p++)
//							{
//									long seed1 = 2342*p+56*5+3;
//									long seed2 = 45*p+556*1+45;	
//									int x1 = (ran2(&seed1)*N);
//									int x2 = (ran2(&seed2)*N);
//									for (int t=0;t<T;t++)
//									{
//										X1[t]=X_e[x1][t];
//										X2[t]=X_e[x2][t];
//										
//									}							
//									
//									mean_corr = mean_corr+1/((double)Pairings)*coherence(X1,X2,temporal_window);
//							}
//						
		double sumr_baseline=0;double var_baseline=0;
		int T_averaging = 20000;
		for (int t=T_min;t<T_min+T_averaging;t++)
		{
			sumr_baseline=sumr_baseline+1/((double)T_averaging)*rate_e[t];
		}
		for (int t=T_min;t<T_min+T_averaging;t++)
		{
			var_baseline=var_baseline+1/((double)T_averaging)*(rate_e[t]-sumr_baseline)*(rate_e[t]-sumr_baseline);
		}
		cout<<"	"<<sumr_baseline<<"	"<<(var_baseline)<<endl;
		
		double sum_stabilized=0;double var_stabilized=0;
		for (int t=T-T_averaging;t<T;t++)
		{
			sum_stabilized=sum_stabilized+1/((double)T_averaging)*rate_e[t];
		}
		for (int t=T-T_averaging;t<T;t++)
		{
			var_stabilized=var_stabilized+1/((double)T_averaging)*(rate_e[t]-sum_stabilized)*(rate_e[t]-sum_stabilized);
		}
		cout<<"	"<<sum_stabilized<<"	"<<(var_stabilized)<<endl;
					
					
		
/////////////OUTPUT//////////////////////////////////////////////////////////////////////////	
	
			
        ofstream outfile;
        
        
     	    outfile.open("ADM  -  Mean corelations and stability metric.txt", ios::out);
         for(int k=0;k<no_of_intervals;k++)
         {
         	
         		  			outfile<<times[k]<<"	"<<mean_rate[k]<<"	"<<mean_correlation[k]<<"	"<<mean_stability_metric[k]<<"	"<<mean_conditional_firing_rate[k]<<endl;
         	
                         
         }  
        outfile.close(); 
     	 
	
        
          outfile.open("ADM  -  membrane potential dynamica.txt", ios::out);
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
        
        
        outfile.open("ADM  -  Single connection.txt", ios::out);
         for(int t=0;t<T;t++)
         {
         	
         		  			outfile<<t<<"	"<<u_e[sampled_index_j[1]][t]<<"	"<<X_e[sampled_index_j[1]][t]<<"	"<<r_e[sampled_index_j[1]][t]<<"	"<<CV_in_time[1][t]<<"	"<<Delay_in_time[1][t]<<endl;
         	
                         
         }  
        outfile.close(); 
        
          outfile.open("ADM  -  Mean activity and Firing rate.txt", ios::out);
         for(int t=0;t<T;t++)
         {
         	
         		  			outfile<<t<<"	"<<mean_u_e[t]<<"	"<<rate_e[t]<<"	"<<u_e[0][t]<<"	"<<f(mean_u_e[t],h)<<endl;
         	
                         
         }  
        outfile.close(); 
        
         outfile.open("ADM  -  Individual Firing rates.txt", ios::out);
         for(int t=0;t<T;t++)
         {
         	
         		  			outfile<<t<<"	"<<rates[1][t]<<"	"<<rates[10][t]<<"	"<<rates[22][t]<<"	"<<rate_e[t]<<endl;
         	
                         
         }  
        outfile.close(); 
        
        outfile.open("ADM  -  Spike response .txt", ios::out);
         for(int t=0;t<T;t++)
         {
         	
         		  			outfile<<t<<"	"<<u_theo[t]<<"	"<<u_e[1][t]<<"	"<<endl;
         	
                         
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
        
    
        outfile.open("ADM  -  Connectivity.txt", ios::out);
         for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  			outfile<<i<<"	"<<j<<"	"<<c[i][j]<<"	"<<tau[i][j]<<"	"<<w[i][j]<<endl;
         		
			}
                         
         }  
        outfile.close(); 
        
        
        outfile.open("ADM  -  theoretical predictions.txt", ios::out);
         for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         					double theo = (alpha_retraction*c_min+sqrt(alpha_retraction*c_min*alpha_retraction*c_min+4*0.3*alpha_formation*alpha_retraction*l[i][j]*0.5*(f(0.03,h)+f(0.2,h))))/(2*alpha_retraction);
         		  			outfile<<i<<"	"<<j<<"	"<<c[i][j]<<"	"<<theo<<"	"<<delay[i][j]<<"	"<<l[i][j]/theo<<endl;
         		
			}
                         
         }  
        outfile.close(); 
        
        
               
    
       	outfile.open("ADM  - Initial and Final Delay distribution.txt", ios::out);
		
		
		
				for (int k=0;k<delay_graining;k++)
				{
				 
				 	
		         					outfile<<(max_delay)/((double) delay_graining)*k<<"	"<<Delay_distribution_initial[k]<<"	"<<Delay_distribution_final[k]<<endl;
		         		 
		     	}
		     	
   
        outfile.close(); 
        
        	outfile.open("ADM  - Initial and Final Conduction Velocity distribution.txt", ios::out);
		
		
		
				for (int k=0;k<speed_graining;k++)
				{
				 
				 	
		         					outfile<<(max_speed)/((double) speed_graining)*k<<"	"<<Speed_distribution_initial[k]<<"	"<<Speed_distribution_final[k]<<endl;
		         		 
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




