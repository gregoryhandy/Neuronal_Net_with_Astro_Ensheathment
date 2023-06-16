/*
** Mex code for the spiking simulations of the spatial model found in Handy and Borisyuk (2023)
** It interfaces with Matlab (see that code for additional details)
*/



#include "mex.h"
#include "math.h"
#include "time.h"
#include "matrix.h"


/* A fast approximation of the exp function */
/* If these lines cause compilation errors, just
   delete them and replace instances of EXP with exp
   below.  It will be slightly slower, but no biggie. */
static union 
{
	  double d;
	    struct {
#ifdef LITTLE_ENDIAN
	    int j,i;
#else 
		    int i,j;
#endif
	  } n;
} _eco;
#define EXP_A (1048576/0.69314718055994530942)
#define EXP_C 60801
#define EXP(y) (_eco.n.i = EXP_A*(y) + (1072693248 - EXP_C), _eco.d)



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	int N,Ne,Ni,Ne1,Ni1,Nt, maxns, Nrecord,Ntref[2],ns,total_num_synapses,total_types_synapses;
	int i,j,k,jj,m1,m2,flag;
	double dt,*s,*v,*v0, *etare,*etari;
	double *eta_adjusted_e0,*eta_adjusted_e1,*eta_adjusted_i0,*eta_adjusted_i1;
	double *eta_vec_next_e0,*eta_vec_next_e1,*eta_vec_next_i0,*eta_vec_next_i1;
	double *tausyn,*vr;
	double *Ix1e,*Ix2e,*Ix1i,*Ix2i;
	double *J,T,*C,*Vleak,*DeltaT,*VT,*tref,*gl,*Vth,*Vre,*Vlb;
	int *refstate;

	int *Irecord,*synapse_num_default;
	int *id_start_firing, *firing_neuron, id_current, *target_neuron,*ensh_realizations;

	// new for spatial
	double *JnextX,xloc,yloc;
	int printfperiod, Nsx,Nx1,Nx, j1, j2;
	double *alphax,*alphaxr,*sx,Jex,Jix,tausynx;
	int iXspike,jspike,postcell;
	int *Wex1,*Wex2,*Wix1,*Wix2;
	int Kex, Kix;
	
	mxArray *temp0;

/******
 * Import variables from matlab
 * This is messy looking and is specific to mex.
 * Ignore if you're implementing this outside of mex.
 *******/
	
// sx,Nx1,Ne1,Ni1,Jex,Jix  0-5
sx = mxGetPr(prhs[0]);
m1 = mxGetM(prhs[0]);
Nsx = mxGetN(prhs[0]);
if(m1!=3){
    mexErrMsgTxt("sx should be Nsxx3");
}

/* Number of neurons in the ffwd layer in each direction. */ 
Nx1=(int)mxGetScalar(prhs[1]);
/* Number of exc neurons in each direction. */ 
Ne1=(int)mxGetScalar(prhs[2]);
/* Number of inh neurons in each direction. */ 
Ni1=(int)mxGetScalar(prhs[3]);

// strength of feedforward connections
Jex=mxGetScalar(prhs[4]);
Jix=mxGetScalar(prhs[5]);

// J,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb, 6-15
J = mxGetPr(prhs[6]);  // Connection weights
total_num_synapses = mxGetM(prhs[6]);
m2 = mxGetN(prhs[6]);
if(total_num_synapses==1 && m2!=1)
    total_num_synapses=m2;
C=mxGetPr(prhs[7]);
m1 = mxGetN(prhs[7]);
m2 = mxGetM(prhs[7]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
gl=mxGetPr(prhs[8]);
m1 = mxGetN(prhs[8]);
m2 = mxGetM(prhs[8]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
Vleak=mxGetPr(prhs[9]);
m1 = mxGetN(prhs[9]);
m2 = mxGetM(prhs[9]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
DeltaT=mxGetPr(prhs[10]);
m1 = mxGetN(prhs[10]);
m2 = mxGetM(prhs[10]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
VT=mxGetPr(prhs[11]);
m1 = mxGetN(prhs[11]);
m2 = mxGetM(prhs[11]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
tref=mxGetPr(prhs[12]);
m1 = mxGetN(prhs[12]);
m2 = mxGetM(prhs[12]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
Vth=mxGetPr(prhs[13]);
m1 = mxGetN(prhs[13]);
m2 = mxGetM(prhs[13]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
Vre=mxGetPr(prhs[14]);
m1 = mxGetN(prhs[14]);
m2 = mxGetM(prhs[14]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
Vlb=mxGetPr(prhs[15]);
m1 = mxGetN(prhs[15]);
m2 = mxGetM(prhs[15]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");


// tausyn, s_en, tausynx,V0,T,dt,maxns,Irecord 16-23
tausyn=mxGetPr(prhs[16]); // synaptic decay parameter (varies between synapses)
total_types_synapses = mxGetM(prhs[16]);
m2 = mxGetN(prhs[16]);
if(total_types_synapses==1 && m2!=1)
    total_types_synapses=m2;

double s_en;
s_en=mxGetScalar(prhs[17]); 

tausynx=mxGetScalar(prhs[18]);

v0 = mxGetPr(prhs[19]);
N = mxGetM(prhs[19]);
m2 = mxGetN(prhs[19]);
if(N==1 && m2!=1)
    N=m2;

T = mxGetScalar(prhs[20]);
dt = mxGetScalar(prhs[21]);

maxns = ((int)mxGetScalar(prhs[22]));


// These are vectors of integers, which is important
// since they will be used to access array locations
// As a result, the mex call has been motified
mxArray *temp[1];
mxArray *lhs[1]; 

temp[0] = (mxArray *) prhs[23];
mexCallMATLAB(1,lhs,1,&temp[0],"int32");
Irecord=(int*)mxGetData(lhs[0]);
Nrecord = mxGetN(prhs[23]);
m1 = mxGetM(prhs[23]);
if(Nrecord==1 || m1!=2)
    mexErrMsgTxt("Irecord should be Nrecordx2 or 1xNrecord");

// Irecord=mxGetPr(prhs[39]);
// Nrecord = mxGetN(prhs[39]);
// m2 = mxGetM(prhs[39]);
// if(m2!=2)
//     mexErrMsgTxt("Irecord should be Nx2.");


// id_start_firing, firing_neuron, synapse_num_default 24-26
temp[0] = (mxArray *) prhs[24];
mexCallMATLAB(1,lhs,1,&temp[0],"int32");
id_start_firing=(int*)mxGetData(lhs[0]); // synaptic decay parameter (varies between synapses)
m1 = mxGetN(prhs[24]);
m2 = mxGetM(prhs[24]);
if(m1!=N && m2!=N)
    mexErrMsgTxt("id_start_rec is not the right size");
temp[0] = (mxArray *) prhs[25];
mexCallMATLAB(1,lhs,1,&temp[0],"int32");
firing_neuron=(int*)mxGetData(lhs[0]);   // synaptic decay parameter (varies between synapses)
m1 = mxGetN(prhs[25]);
m2 = mxGetM(prhs[25]);
if(m1!=(total_num_synapses+1) && m2!=(total_num_synapses+1))
    mexErrMsgTxt("firing_neuron is not the right size");
temp[0] = (mxArray *) prhs[26];
mexCallMATLAB(1,lhs,1,&temp[0],"int32");
synapse_num_default = (int*)mxGetData(lhs[0]);
m1 = mxGetN(prhs[26]);
m2 = mxGetM(prhs[26]);
if(m1!=total_num_synapses && m2!=total_num_synapses)
    mexErrMsgTxt("synapse_num_default is not the right size");

// target_neuron, ensh_realizations 27-28
temp[0] = (mxArray *) prhs[27];
mexCallMATLAB(1,lhs,1,&temp[0],"int32");
target_neuron = (int*)mxGetData(lhs[0]);
m1 = mxGetN(prhs[27]);
m2 = mxGetM(prhs[27]);
if(m1!=total_num_synapses && m2!=total_num_synapses)
    mexErrMsgTxt("target_neuron is not the right size");


temp[0] = (mxArray *) prhs[28];
mexCallMATLAB(1,lhs,1,&temp[0],"int32");
ensh_realizations = (int*)mxGetData(lhs[0]);
m1 = mxGetN(prhs[28]);
m2 = mxGetM(prhs[28]);
if(m1!=total_num_synapses && m2!=total_num_synapses)
    mexErrMsgTxt("ensh_realizations is not the right size");

// Kex, Kix, Wex1, Wex2, Wix1, Wix2 29-32

/* Number of connections for the ffwd layer to excitatory cells */ 
Kex=(int)mxGetScalar(prhs[29]);
/* Number of connections for the ffwd layer to inhibitory cells  */ 
Kix=(int)mxGetScalar(prhs[30]);

/* Total number of ffwd neurons */
Nx = Nx1*Nx1;

temp[0] = (mxArray *) prhs[31];
mexCallMATLAB(1,lhs,1,&temp[0],"int32");
Wex1 = (int*)mxGetData(lhs[0]);
m1 = mxGetN(prhs[31]);
m2 = mxGetM(prhs[31]);
if(m1!=Nx*Kex && m2!=Nx*Kex)
    mexErrMsgTxt("Wex1 is not the right size");

temp[0] = (mxArray *) prhs[32];
mexCallMATLAB(1,lhs,1,&temp[0],"int32");
Wex2 = (int*)mxGetData(lhs[0]);
m1 = mxGetN(prhs[32]);
m2 = mxGetM(prhs[32]);
if(m1!=Nx*Kex && m2!=Nx*Kex)
    mexErrMsgTxt("Wex2 is not the right size");

temp[0] = (mxArray *) prhs[33];
mexCallMATLAB(1,lhs,1,&temp[0],"int32");
Wix1 = (int*)mxGetData(lhs[0]);
m1 = mxGetN(prhs[33]);
m2 = mxGetM(prhs[33]);
if(m1!=Nx*Kix && m2!=Nx*Kix)
    mexErrMsgTxt("Wix1 is not the right size");

temp[0] = (mxArray *) prhs[34];
mexCallMATLAB(1,lhs,1,&temp[0],"int32");
Wix2 = (int*)mxGetData(lhs[0]);
m1 = mxGetN(prhs[34]);
m2 = mxGetM(prhs[34]);
if(m1!=Nx*Kix && m2!=Nx*Kix)
    mexErrMsgTxt("Wix2 is not the right size");

/******
 * Finished importing variables.
 *******/

/* Total number of each type of neuron */
Ne=Ne1*Ne1;
Ni=Ni1*Ni1;

/* Check for consistency with total number of neurons */
if(N!=Ne+Ni)
    mexErrMsgTxt("Ne1 and/or Ni1 not consistent with size of V0");

/* Numebr of time bins */
Nt=(int)(T/dt);

/******
 * Now allocate new variables.
 * This is also mex specific.  Use malloc in C, etc.
 *****/

/* Allocate output vector */
plhs[0] = mxCreateDoubleMatrix(3, maxns, mxREAL);
s=mxGetPr(plhs[0]);

/* Allocate output vector */
plhs[1] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
alphaxr=mxGetPr(plhs[1]);


/* Allocate output vector */
plhs[2] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
etare=mxGetPr(plhs[2]);
/* Allocate output vector */
plhs[3] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
etari=mxGetPr(plhs[3]);

/* Allocate output vector */
plhs[4] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
vr=mxGetPr(plhs[4]);

/* Allocate membrane potential */
temp0=mxCreateDoubleMatrix(N, 1, mxREAL);
v = mxGetPr(temp0);
// Another possibility: v=mxMalloc(N*sizeof(double));

// four types of inputs
eta_vec_next_e0=mxMalloc(N*sizeof(double));
eta_vec_next_e1=mxMalloc(N*sizeof(double));
eta_vec_next_i0=mxMalloc(N*sizeof(double));
eta_vec_next_i1=mxMalloc(N*sizeof(double));

eta_adjusted_e0=mxMalloc(N*sizeof(double));
eta_adjusted_e1=mxMalloc(N*sizeof(double));
eta_adjusted_i0=mxMalloc(N*sizeof(double));
eta_adjusted_i1=mxMalloc(N*sizeof(double));

// Feedforward inputs
JnextX=mxMalloc(N*sizeof(double));
alphax=mxMalloc(N*sizeof(double));

refstate=mxMalloc(N*sizeof(int));

/*****
 * Finished allocating variables
 ****/

/* Inititalize variables */
for(j=0;j<N;j++){
    v[j]=v0[j]; 
    refstate[j]=0;
	
	eta_adjusted_e0[j]=0;
	eta_adjusted_e1[j]=0;
	eta_adjusted_i0[j]=0;
	eta_adjusted_i1[j]=0;
	
	eta_vec_next_e0[j]=0;
	eta_vec_next_e1[j]=0;
	eta_vec_next_i0[j]=0;
	eta_vec_next_i1[j]=0;
	
    alphax[j]=0;
	JnextX[j]=0;        
}


/* Record first time bin */
for(jj=0;jj<Nrecord;jj++){
      /* Find index into local variables */
    
        /* Find index into local variables */
        j1=(int)round(Irecord[2*jj]-1);
        j2=(int)round(Irecord[2*jj+1]-1);

        if(j1<Ne1 && j2<Ne1){                 
            j=j1+Ne1*j2;
        }
        else
          if(j1>=Ne1 && j2>=Ne1){
            j=(j1-Ne1)+(j2-Ne1)*Ni1+Ne;
          }
          else
            mexErrMsgTxt("Indices in Irecord must have both terms <Ne1 or both terms >Ne1");


                
      if(j>=N || j<0){         
         mexErrMsgTxt("Irecord contains out of bounds indices.");
      }

  	  etare[jj+Nrecord*0]=eta_adjusted_e0[j]+eta_adjusted_e1[j];
  	  etari[jj+Nrecord*0]=eta_adjusted_i0[j]+eta_adjusted_i1[j];
      alphaxr[jj+Nrecord*0]=alphax[j];
      vr[jj+Nrecord*0]=v[j];
}
    
/* Initialize connections: now occurs in MATLAB */


/* Refractory states */
Ntref[0]=(int)round(tref[0]/dt);
Ntref[1]=(int)round(tref[1]/dt);

printfperiod=(int)(round(Nt/20.0));

/* Initialize number of spikes */
ns=0;
flag=0;
iXspike=0;


/* Main loop */
/* Exit loop and issue a warning if max number of spikes is exceeded */
for(i=1;i<Nt && ns<maxns;i++){
      
    
	/* Find all spikes in feedforwar layer at this time bin */
	/* Add to corresponding elements of JnextX */
	while(sx[iXspike*3+0]<=i*dt && iXspike<Nsx){
		jspike=(int)round(sx[iXspike*3+1]-1)*Nx1+(int)round(sx[iXspike*3+2]-1);
		if(jspike<0 || jspike>=Nx){
			mexPrintf("\n %d %d %d %d %d %d\n",(int)round(sx[iXspike*3+0]/dt),iXspike,i,jspike,(int)round(sx[iXspike*3+1]-1),(int)round(sx[iXspike*3+2]-1));
			mexErrMsgTxt("Out of bounds index in sx.");
		}
		
		for(k=0;k<Kex;k++){ 
			if(jspike*Kex+k>=Nx*Kex || jspike*Kex+k<0){
				mexErrMsgTxt("Out of bounds jspike E.");
			}
			
			postcell=Wex1[jspike*Kex+k]*Ne1+Wex2[jspike*Kex+k];    
			if(postcell<0 || postcell>=N){
				mexErrMsgTxt("Out of bounds index ino JnextX");
			}	
			
			JnextX[postcell]+=Jex;
		}
		
		for(k=0;k<Kix;k++){
			if(jspike*Kix+k>=Nx*Kix || jspike*Kix+k<0){
				mexErrMsgTxt("Out of bounds jspike I.");
			}
				
			postcell=Ne+Wix1[jspike*Kix+k]*Ni1+Wix2[jspike*Kix+k];
			if(postcell<0 || postcell>=N){
				mexErrMsgTxt("Out of bounds index ino JnextX i");
			}
				
			JnextX[postcell]+=Jix;
		}

		iXspike++;
	}

	// loop through each neuron
	for(j=0;j<N;j++){    
    
        
		/* Update synaptic variables */
		//alphae[j]-=alphae[j]*(dt/tausyne);
		//alphai[j]-=alphai[j]*(dt/tausyni);  
		alphax[j]-=alphax[j]*(dt/tausynx);  
		 
		// Uses Euler's method to approximate the exponential decay     		 		
		eta_adjusted_e0[j] = eta_adjusted_e0[j]-eta_adjusted_e0[j]*(dt/tausyn[0]);
		eta_adjusted_e1[j] = eta_adjusted_e1[j]-eta_adjusted_e1[j]*(dt/tausyn[1]);
		eta_adjusted_i0[j] = eta_adjusted_i0[j]-eta_adjusted_i0[j]*(dt/tausyn[2]);
		eta_adjusted_i1[j] = eta_adjusted_i1[j]-eta_adjusted_i1[j]*(dt/tausyn[3]);
      
		if(j<Ne){
             
			// if the cell is not in the refractory state
			if(refstate[j]<=0){
				//v[j]+=fmax((alphae[j]+alphai[j]+alphax[j]-gl[0]*(v[j]-Vleak[0])+gl[0]*DeltaT[0]*EXP((v[j]-VT[0])/DeltaT[0]))*dt/C[0],Vlb[0]-v[j]);
				
				v[j]=v[j]+(eta_adjusted_e0[j]+eta_adjusted_e1[j]+eta_adjusted_i0[j]+eta_adjusted_i1[j]+alphax[j]
					-gl[0]*(v[j]-Vleak[0])
						+gl[0]*DeltaT[0]*exp((v[j]-VT[0])/DeltaT[0]))*dt/C[0];
			 
				v[j]=fmax(v[j],Vlb[0]);
			}
			else{       
				// Hold the voltage at threshold for the refractory time
				// and drop it to Vre at the last time point          
				if(refstate[j]>1){
					v[j]=Vth[0];
				}
				else{
					v[j]=Vre[0];
				}
				 
				refstate[j]--;
			}
             
             
			/* If a spike occurs */
			if(v[j]>=Vth[0] && refstate[j]<=0 && ns<maxns){
				refstate[j]=Ntref[0];
				v[j]=Vth[0];       /* reset membrane potential */
				s[0+3*ns]=i*dt; /* spike time */
				s[2+3*ns]=j/Ne1+1;     /* neuron index 1 */
				s[1+3*ns]=j%Ne1+1;     /* neuron index 2 */
				ns++;           /* update total number of spikes */


				// For each postsynaptic target, propagate spike into JnextE   
				id_current=id_start_firing[j];
				while(firing_neuron[id_current] == j){
					if(ensh_realizations[synapse_num_default[id_current]]==0){
						eta_vec_next_e0[target_neuron[id_current]] = eta_vec_next_e0[target_neuron[id_current]]
							+J[synapse_num_default[id_current]];
						
					}else{
						eta_vec_next_e1[target_neuron[id_current]] = eta_vec_next_e1[target_neuron[id_current]]
							+J[synapse_num_default[id_current]]*(1-s_en);
					}

					id_current = id_current + 1;
				}
			}
             
                 
		}
		else{ /* If cell is inhibitory */
            
			// if the cell is not refractory
			if(refstate[j]<=0){
				//v[j]+=fmax((alphae[j]+alphai[j]+alphax[j]-gl[1]*(v[j]-Vleak[1])+
				//     gl[1]*DeltaT[1]*EXP((v[j]-VT[1])/DeltaT[1]))*dt/C[1],Vlb[1]-v[j]);
				
				v[j]=v[j]+(eta_adjusted_e0[j]+eta_adjusted_e1[j]+eta_adjusted_i0[j]+eta_adjusted_i1[j]+alphax[j]
					-gl[1]*(v[j]-Vleak[1])
						+gl[1]*DeltaT[1]*exp((v[j]-VT[1])/DeltaT[1]))*dt/C[1];
				
				v[j]=fmax(v[j],Vlb[1]);
			}
			else{                 
				if(refstate[j]>1){
					v[j]=Vth[1];
				}
				else{
					v[j]=Vre[1];
				}
			   
				refstate[j]--;
			}
             
             
			/* If a spike occurs */
			if(v[j]>=Vth[1] && refstate[j]<=0 && ns<maxns){                                                            
                  
                  
				refstate[j]=Ntref[1];
				v[j]=Vth[1];        /* reset membrane potential */
				s[0+3*ns]=i*dt;     /* spike time */
				s[2+3*ns]=-((j-Ne)/Ni1)-1;     /* neuron index 1 */
				s[1+3*ns]=-((j-Ne)%Ni1)-1;     /* neuron index 2 */
                 
				ns++;           /* update total number of spikes */


				// For each postsynaptic target             
				id_current=id_start_firing[j];
				while(firing_neuron[id_current] ==j){
					
					
					if(ensh_realizations[synapse_num_default[id_current]]==0){
						eta_vec_next_i0[target_neuron[id_current]] = eta_vec_next_i0[target_neuron[id_current]]
							+J[synapse_num_default[id_current]];
					}else{
						eta_vec_next_i1[target_neuron[id_current]] = eta_vec_next_i1[target_neuron[id_current]]
							+J[synapse_num_default[id_current]]*(1-s_en);
					}           
					id_current = id_current + 1;  
				}
			}
		}
	}
  

	/* Store recorded variables */
	for(jj=0;jj<Nrecord;jj++){
            
		/* Find index into local variables */
		j1=(int)round(Irecord[2*jj+0]-1);
		j2=(int)round(Irecord[2*jj+1]-1);

		if(j1<Ne1 && j2<Ne1){                 
			j=j1+Ne1*j2;
		}
		else
		if(j1>=Ne1 && j2>=Ne1){
			j=(j1-Ne1)+(j2-Ne1)*Ni1+Ne;
		}
		else
			mexErrMsgTxt("Indices in Irecord must have both terms <Ne1 or both terms >Ne1");

		etare[jj+Nrecord*i]=eta_adjusted_e0[j]+eta_adjusted_e1[j];
		etari[jj+Nrecord*i]=eta_adjusted_i0[j]+eta_adjusted_i1[j];
		alphaxr[jj+Nrecord*i]=alphax[j];
		vr[jj+Nrecord*i]=v[j]; 
			  
		// alphaer[jj+Nrecord*i]=alphae[j];
		// alphair[jj+Nrecord*i]=alphai[j];
		// alphaxr[jj+Nrecord*i]=alphax[j];
		// vr[jj+Nrecord*i]=v[j];
	}    

	/* Propagate spikes */
	for(j=0;j<N;j++){
		eta_adjusted_e0[j] = eta_adjusted_e0[j]+eta_vec_next_e0[j]/tausyn[0];
		eta_adjusted_e1[j] = eta_adjusted_e1[j]+eta_vec_next_e1[j]/tausyn[1];
		eta_adjusted_i0[j] = eta_adjusted_i0[j]+eta_vec_next_i0[j]/tausyn[2];
		eta_adjusted_i1[j] = eta_adjusted_i1[j]+eta_vec_next_i1[j]/tausyn[3];
		
		alphax[j]+=JnextX[j]/tausynx;
		
		eta_vec_next_e0[j] = 0;
		eta_vec_next_e1[j] = 0;
		eta_vec_next_i0[j] = 0;
		eta_vec_next_i1[j] = 0;
			
		JnextX[j]=0;
	}


	// if(i%printfperiod==0){
	// 	mexPrintf("%d percent complete  rate = %2.2fHz\n",i*100/Nt,1000*((double)(ns))/(((double)(N))*((double)(i))*dt));
	// 	mexEvalString("drawnow;");
	// }

}


/* Issue a warning if max number of spikes reached */
if(ns>=maxns)
   mexWarnMsgTxt("Maximum number of spikes reached, simulation terminated.");

/* Free allocated memory */
mxDestroyArray(temp0);
// mxFree(v);
mxFree(JnextX);
mxFree(alphax);
mxFree(refstate);

mxFree(eta_adjusted_e0);
mxFree(eta_adjusted_e1);
mxFree(eta_adjusted_i0);
mxFree(eta_adjusted_i1);

mxFree(eta_vec_next_i0);
mxFree(eta_vec_next_i1);
mxFree(eta_vec_next_e0);
mxFree(eta_vec_next_e1);

}










