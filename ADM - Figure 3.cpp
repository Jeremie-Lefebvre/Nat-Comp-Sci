
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
#define twindow 200//Time window for firing rate calculations
#define T_min 1000 ///minimal time post transient to initiate statistic calculations 
#define Q 10 //number of trials
#define Sample 10 ////Number of nodes selected for plotting purposes
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

#define T_plastic_on T//plasticity onset time (post transient)
#define T_plastic_off T//plasticity offset time

#define t_delay_statistics_initial T_min+100//initial time in which delay statistics are computed
#define t_delay_statistics_final T-1000//final time in which delay statistics are computed
#define delay_graining 50//number of bins for delay statistics
#define speed_graining 50//number of bins for conduction velocity statistics
#define max_delay 200//maximal delay for delay statistics
#define max_speed 30//maximum conduction velocity for conduction velocity statistics

double f(double u, double threshold);
double f_prime(double u, double threshold);
double Heaviside(double input);
float ran2(long *idum);//random number generator - initial cond.
void shuffle_indices(int Array[], int size);//shuffle an array with int  entries
void shuffle(double Array[], int size);//shuffle an array with double  entries
void sort(double Array[], int size);//sort array 
void four1(double data[], unsigned long nn, int isign);
double coherence(double X1[], double X2[], int window);



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


double MEAN[Q] ;//mean firing network rate 
double VARIANCE[Q]; //mean network firing rate variance in time
double CORR[Q]; //mean network firing rate correlation in time
double LOCAL_VARIANCES[Q];//local single neuron variance in firing rate in time


using namespace std;


int main()
{
	
srand (time(NULL));

double U_max=1;
double U_min=-1;
int K =10000;
double fixed_point=0;
for (int k=0;k<K;k++)
{
	double u1=U_max-k*fabs(U_max-U_min)/((double)K);
	double u2=U_max-(k+1)*fabs(U_max-U_min)/((double)K);
	double G1 = u1 - w_o*f(u1,h);
	double G2 = u2 - w_o*f(u2,h);

	if (G1*G2<0)
	{
		fixed_point=0.5*u1+0.5*u2;
	}
	
	
}


		cout<<"Importing data..."<<endl;
		ifstream sfile("LENGTHS100.txt", ios::out);
         for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  			sfile>>l[i][j];
         		
			}
                         
         }  
        
        
      ifstream ufile("CVS1000.txt", ios::out);
     for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  			ufile>>c[i][j];
         		//c[i][j]=c_max;
			}
                         
         }  
        
         
//          for(int i=0;i<N;i++)
//         {
//         	for (int j=0;j<N;j++)
//         	{
//         	
//         			
//         		  		c[i][j]=c_min;
//         		
//			}
//                         
//         }  
         
         
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
				
							
					
					// define noise arrays
		        	for (int t=0;t<T;t++)
					{
						
						double sum=0;
						for(int n=0;n<N;n++)
		        		{
		                            long d=rand()%105;      
		                            long noisyseed1=(long)21*t+n*n*187198+56*d+12*1+1;
		                            long noisyseed2=(long)69*t*n+11*n+t+45*d+2*1+1; 
		                          
		            				 
		            				u_e[n][t]=-0.1;//sqrt(-2*log(ran2(&noisyseed1)))*cos(2*3.14159*ran2(&noisyseed2));
		            				
						}
					
		           }
		               
				
		
		
		 		
						
					

			
				
			
						
				
////////////////////EULER SCHEME//////////////////////////////////////////////////////////////////////////				
	for (int q=0;q<Q;q++)
	{
					MEAN[q]=0;
					VARIANCE[q]=0;
					CORR[q] = 	0;
					LOCAL_VARIANCES[q] = 0;
			cout<<q<<endl;
								// define noise arrays
		        	for (int t=0;t<T;t++)
					{
						
						double sum=0;
						for(int n=0;n<N;n++)
		        		{
		                            long d=rand()%105;      
		                            long noisyseed1=(long)21*t+n*n*187198+76*q+56*d+12*1+1;
		                            long noisyseed2=(long)69*t*n+11*n*q+t+45*d+2*1+1; 
		                          
		            				 
		            				u_e[n][t]=sqrt(-2*log(ran2(&noisyseed1)))*cos(2*3.14159*ran2(&noisyseed2));
		            				X_e[n][t]=0;
		            				r_e[n][t]=0;
		            				
						}
						rate_e[t];
						
		           }
		               
					for (int t=T_min;t<T;t++)
					{
					
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
						   				  long seed_e = (long) 12+i*i+t+q*10+15+t*seed_q+56*i+i+99*t+2;              
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
										   		
										   for (int t_index=0;t_index<twindow;t_index++)
										   {
										   			r_e[i][t]=r_e[i][t]+1/((double)twindow*Dt)*X_e[i][t+1-t_index]*dt;
										   			
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
									double c_now = c[i][j];
									//c[i][j]=c[i][j]+dt*alpha_c*Heaviside(w[i][j])*(0.5*exp(0.03*r_e[j][t-1])-c[i][j]);
									//c[i][j]=c[i][j]+dt*(-(c[i][j])*alpha_retraction+0.1*X_e[j][t-1]*(1-c[i][j]/c_max)*alpha_formation);//kappa*(r-1));
									c[i][j]=c[i][j]+dt*((c_min-c[i][j])*alpha_retraction+0.3*X_e[j][t-1]*(l[i][j]/(c[i][j]))*alpha_formation);//kappa*(r-1));
										
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
						//coincidence_theo[t]=1/dt*1/dt*( 1-exp(-f(mean_u_e[t+1],h)*dt))*( 1-exp(-f(mean_u_e[t+1],h)*dt));//
//						
//						if (t==t_delay_statistics_initial)
//						{
//						 cout<<t<<endl;
//													for (int s=0;s<delay_graining;s++)
//													{
//																			Delay_distribution_initial[s]=0;
//													}
//																					
//													for (int tt=0;tt<delay_graining;tt++)
//													{
//																			double ttau = 0 + (max_delay)/((double) delay_graining)*tt;
//																			double dttau = (max_delay)/((double) delay_graining);
//																							
//																			for (int i=0;i<N;i++)
//																			{
//																				for (int j=0;j<N;j++)
//																				{
//																					
//																									double ttau = 0 + (max_delay)/((double) delay_graining)*tt;
//																									double dttau = (max_delay)/((double) delay_graining);
//																									double delay  = (double) tau[i][j];			
//																									if (((double)delay>ttau-dttau*0.5)&&((double)delay<ttau+dttau*0.5))
//																									{
//																												Delay_distribution_initial[tt] = Delay_distribution_initial[tt]+w[i][j];
//													
//																									}	
//																				}
//																			}
//													}
//									
//														
//										for (int s=0;s<speed_graining;s++)
//										{
//													Speed_distribution_initial[s]=0;
//										}
//										for (int tt=0;tt<speed_graining;tt++)
//										{
//																double tc = 0 + (max_speed)/((double) speed_graining)*tt;
//																double dcc = (max_speed)/((double) speed_graining);
//																
//																for (int i=0;i<N;i++)
//																{
//													
//																		for (int j=0;j<N;j++)
//																		{
//																			double tc = 0 + (max_speed)/((double) speed_graining)*tt;
//																			double dcc = (max_speed)/((double) speed_graining);
//																			double speed = (double) c[i][j];			
//																			if (((double)speed>tc-dcc*0.5)&&((double)speed<tc+dcc*0.5))
//																			{
//																						Speed_distribution_initial[tt] = Speed_distribution_initial[tt]+w[i][j];
//																						
//																						
//																			}	
//																		}
//															
//																}
//																
//										}
//										
//						}//initial distribution loop
								
//						if (t==t_delay_statistics_final)
//						{
//								cout<<t<<endl;
//													for (int s=0;s<delay_graining;s++)
//													{
//																			Delay_distribution_final[s]=0;
//													}
//																					
//													for (int tt=0;tt<delay_graining;tt++)
//													{
//																			double ttau = 0 + (max_delay)/((double) delay_graining)*tt;
//																			double dttau = (max_delay)/((double) delay_graining);
//																							
//																			for (int i=0;i<N;i++)
//																			{
//																				for (int j=0;j<N;j++)
//																				{
//																					
//																									double ttau = 0 + (max_delay)/((double) delay_graining)*tt;
//																									double dttau = (max_delay)/((double) delay_graining);
//																									double delay  = (double) tau[i][j];			
//																									if (((double)delay>ttau-dttau*0.5)&&((double)delay<ttau+dttau*0.5))
//																									{
//																												Delay_distribution_final[tt] = Delay_distribution_final[tt]+w[i][j];
//													
//																									}	
//																				}
//																			}
//													}
//									
//														
//										for (int s=0;s<speed_graining;s++)
//										{
//													Speed_distribution_final[s]=0;
//										}
//										for (int tt=0;tt<speed_graining;tt++)
//										{
//																double tc = 0 + (max_speed)/((double) speed_graining)*tt;
//																double dcc = (max_speed)/((double) speed_graining);
//																
//																for (int i=0;i<N;i++)
//																{
//													
//																		for (int j=0;j<N;j++)
//																		{
//																			double tc = 0 + (max_speed)/((double) speed_graining)*tt;
//																			double dcc = (max_speed)/((double) speed_graining);
//																			double speed = (double) c[i][j];			
//																			if (((double)speed>tc-dcc*0.5)&&((double)speed<tc+dcc*0.5))
//																			{
//																						Speed_distribution_final[tt] = Speed_distribution_final[tt]+w[i][j];
//																						
//																						
//																			}	
//																		}
//															
//																}
//																
//										}
//										
//								}// final distribution loop
									
									
						
					}// t loop
////////////Firing Rate////////////////////////////////////////////////
					
								int range_e=N;
							
								for (int t=twindow;t<T;t++)
								{
									double rate_e_counter=0;
									double rate_i_counter=0;
									double rate_fdbck_counter=0;
									double countse=0;
									double countsi=0;
									double countsfdbck=0;
							
									for (int time=0;time<twindow;time++)
									{
										for (int index=0;index<N;index++)
										{
										if(X_e[index][t+1-time]>0.1)
										{
											rate_e_counter=rate_e_counter+1;///((double) Ne);
											countse = countse+1;
										}
										else{}
										}
									}
									
										rate_e[t]=rate_e_counter/((double) twindow*Dt)*1/((double)range_e);
									}
									
									
								
								
							
								
int Pairings = (int) N*N-N;
double mean_correlation=0;
double mean_correlation_without_mean_field=0;
double mean_local_variances=0;
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
									
										
										double mean_c1=0;double mean_c2=0;double cov_c12=0;double var_c1=0;double var_c2=0;
								        double cov_cc12=0;double var_cc1=0;double var_cc2=0;
										for (int s=0;s<T;s++)
										{
											mean_c1 = mean_c1 +1/((double)T)*r_e[x1][s];
											mean_c2 = mean_c2 +1/((double)T)*r_e[x2][s];
											
										}
										for (int s=0;s<T;s++)
										{
											var_c1 = var_c1 +1/((double)T)*(r_e[x1][s]-rate_e[s])*(r_e[x1][s]-rate_e[s]);
											var_c2 = var_c2 +1/((double)T)*(r_e[x2][s]-rate_e[s])*(r_e[x2][s]-rate_e[s]);
											cov_c12 = cov_c12 +1/((double)T)*(r_e[x1][s]-rate_e[s])*(r_e[x2][s]-rate_e[s]);
											var_cc1 = var_cc1 +1/((double)T)*(r_e[x1][s]-mean_c1)*(r_e[x1][s]-mean_c1);
											var_cc2 = var_cc2 +1/((double)T)*(r_e[x2][s]-mean_c2)*(r_e[x2][s]-mean_c2);
											cov_cc12 = cov_cc12 +1/((double)T)*(r_e[x1][s]-mean_c1)*(r_e[x2][s]-mean_c2);
										
										}
										mean_correlation_without_mean_field = mean_correlation_without_mean_field +1/((double) Pairings)*cov_c12/(sqrt(var_c1*var_c2));
										mean_local_variances = mean_local_variances+1/((double) Pairings)*(0.5*var_c1+0.5*var_c1);
										mean_correlation=mean_correlation +1/((double) Pairings)*cov_cc12/(sqrt(var_cc1*var_cc2));
									
									}
		}

        double mean_ui=0;double var_ui=0;
        for (int t=(int)T/2;t<T;t++)
        {
        	mean_ui=mean_ui+1/((double)T/2)*u_e[3][t];
		}
		for (int t=(int)T/2;t<T;t++)
        {
        	var_ui=var_ui+1/((double)T/2)*(u_e[3][t]-mean_ui)*((u_e[3][t]-mean_ui));
		}
	
		double mean_u[T];
		for (int t=0;t<T;t++)
		{
			mean_u[t]=0;
		}				
		for (int t=T/2;t<T;t++)
		{
			for (int i=0;i<N;i++)
			{
				mean_u[t] =mean_u[t]+1/((double)N)*u_e[i][t];
			}
		}
		double m_u=0;
		double v_u=0;
		for (int t=T/2;t<T;t++)
		{
			m_u = m_u +1/((double)T/2)*mean_u[t];
		}
		for (int t=T/2;t<T;t++)
		{
			v_u = v_u +1/((double)T/2)*(mean_u[t]-m_u)*(mean_u[t]-m_u);
		}
		//cout<<m_u<<"	"<<v_u<<endl;
	
	
	
	
	
		double meanr=0;double var=0;
		for (int t=(int)T/2;t<T;t++)
		{
			meanr=meanr+1/((double)T/2)*rate_e[t];
		}
		for (int t=(int)T/2;t<T;t++)
		{
			var=var+1/((double)T/2)*(rate_e[t]-meanr)*(rate_e[t]-meanr);
		}
		
		MEAN[q] = 	meanr;
		VARIANCE[q] = 	var;
		CORR[q] = 	mean_correlation;
		LOCAL_VARIANCES[q] = mean_local_variances;
		cout<<mean_correlation_without_mean_field<<endl;
		cout<<q<<"	"<<MEAN[q]<<"	"<<VARIANCE[q] <<"	"<<CORR[q]<<"	"<<LOCAL_VARIANCES[q]<<endl;
		
}//q loop

double SD_mean_across_trials =0;
double SD_variance_across_trials =0;
double SD_corr_across_trials =0;
double mean_1=0;double mean_2=0;double mean_3=0;double mean_4=0;
for (int q=0;q<Q;q++)
{
	mean_1=mean_1+1/((double)Q)*MEAN[q];
	mean_2=mean_2+1/((double)Q)*VARIANCE[q];
	mean_3=mean_3+1/((double)Q)*CORR[q];
	mean_4=mean_4+1/((double)Q)*LOCAL_VARIANCES[q];
	
}
double var_1=0;double var_2=0;double var_3=0;double var_4=0;
for (int q=0;q<Q;q++)
{
	var_1=var_1+1/((double)Q)*(MEAN[q]-mean_1)*(MEAN[q]-mean_1);
	var_2=var_2+1/((double)Q)*(VARIANCE[q]-mean_2)*(VARIANCE[q]-mean_2);
	var_3=var_3+1/((double)Q)*(CORR[q]-mean_3)*(CORR[q]-mean_3);
	var_4 = var_4+1/((double)Q)*(LOCAL_VARIANCES[q]-mean_4)*(LOCAL_VARIANCES[q]-mean_4);
}

double Gamma=0;
double s=0.1;
			
for (int i=0;i<N;i++)
{
        				for (int j=0;j<N;j++)
        				{
        					Gamma = Gamma+1/((double)N*N)*exp(-tau[i][j]*dt*s);
						}
}

double fixedPoint_1 = 0.0335;
double S=w_o*w_o/((double)N*(1-w_o*f_prime(fixedPoint_1,h)*Gamma));   

					
double mean_firing_rate_theo =   0.2537354531e-1 * sqrt(-(0.2383736622e1 * S + 0.6684252738e0) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.3805155286e-1 * S - 0.1923133216e-2) / S) / (0.1e1 / 0.3141592654e1 * sqrt(-(0.2383736622e1 * S + 0.6684252738e0) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.3805155286e-1 * S - 0.1923133216e-2) / S) / 0.2e1 + 0.1e1 / 0.3141592654e1 * sqrt(-(0.2101433283e2 * S + 0.1162972721e2) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.7165986182e-1 * S + 0.1904648304e-2) / S) / 0.2e1) + 0.1439796044e0 * sqrt(-(0.2101433283e2 * S + 0.1162972721e2) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.7165986182e-1 * S + 0.1904648304e-2) / S) / (0.1e1 / 0.3141592654e1 * sqrt(-(0.2383736622e1 * S + 0.6684252738e0) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.3805155286e-1 * S - 0.1923133216e-2) / S) / 0.2e1 + 0.1e1 / 0.3141592654e1 * sqrt(-(0.2101433283e2 * S + 0.1162972721e2) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.7165986182e-1 * S + 0.1904648304e-2) / S) / 0.2e1);

					
double mean_variance_theo =  0.1406739725e-1 * sqrt(-(0.2383736622e1 * S + 0.6684252738e0) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.3805155286e-1 * S - 0.1923133216e-2) / S) * pow(0.1e1 / 0.3141592654e1 * sqrt(-(0.2383736622e1 * S + 0.6684252738e0) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.3805155286e-1 * S - 0.1923133216e-2) / S) / 0.2e1 + 0.1e1 / 0.3141592654e1 * sqrt(-(0.2101433283e2 * S + 0.1162972721e2) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.7165986182e-1 * S + 0.1904648304e-2) / S) / 0.2e1, -0.2e1) * sqrt(-(0.2101433283e2 * S + 0.1162972721e2) * (0.1436106292e2 * S - 0.1348106324e1)) * exp(-(-0.7165986182e-1 * S + 0.1904648304e-2) / S);

					
													
double MEAN_FIRING_RATE_THEO = 100*mean_firing_rate_theo;
double MEAN_VARIANCE_THEO =100*100* mean_variance_theo;
double MEAN_CORRELATION_THEO =mean_variance_theo/(mean_variance_theo+0.5*f_prime(fixedPoint_1,h)*f_prime(fixedPoint_1,h)*w_o*w_o*f(fixedPoint_1,h));

cout<<"mean firing rate:"<<mean_1<<endl;
cout<<"mean firing rate SD:"<<sqrt(var_1)<<endl;

cout<<"mean variance:"<<mean_2<<endl;
cout<<"mean variance SD:"<<sqrt(var_2)<<endl;

cout<<"mean correlation:"<<mean_3<<endl;
cout<<"mean correlation SD:"<<sqrt(var_3)<<endl;

cout<<"mean local variances:"<<mean_4<<endl;
cout<<"mean local variances SD:"<<sqrt(var_4)<<endl;

cout<<"mean firing rate theo:"<<MEAN_FIRING_RATE_THEO<<endl;
cout<<"mean correlation theo:"<<MEAN_CORRELATION_THEO<<endl;

/////////////OUTPUT//////////////////////////////////////////////////////////////////////////	
	
			
        ofstream outfile;
        
        
     	       
          outfile.open("ADM  -  Mean variances correlations across trials.txt", ios::out);
         for(int q=0;q<Q;q++)
         {
         	
         		  			outfile<<q<<"	"<<MEAN[q]<<"	"<<VARIANCE[q]<<"	"<<CORR[q]<<""<<LOCAL_VARIANCES[q]<<endl;
         	
                         
         }  
        outfile.close(); 
        
          outfile.open("ADM  -  Individual rates in time.txt", ios::out);
         for(int t=0;t<T;t++)
         {
         	
         		  			outfile<<t<<"	"<<rate_e[t]<<"	"<<r_e[2][t]<<"	"<<r_e[20][t]<<"	"<<r_e[12][t]<<"	"<<r_e[6][t]<<endl;
         	
                         
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
 getch();
      
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




