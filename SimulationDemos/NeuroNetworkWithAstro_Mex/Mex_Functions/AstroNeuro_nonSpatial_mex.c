/*
** Mex code for the spiking simulations of the non-spatial model found in Handy and Borisyuk (2023)
** It interfaces with Matlab (see that code for additional details)
*/

#include "mex.h"
#include "math.h"
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

double s_en;

int *Irecord,*synapse_num_default;
int *id_start_firing, *firing_neuron, id_current, *target_neuron,*ensh_realizations;

double frac_complete;
 
/******
 * Import variables from matlab
 * This is messy looking and is specific to mex.
 * Ignore if you're implementing this outside of mex.
 *******/

// Ix1e,Ix2e,Ix1i,Ix2i,Ne,Ni,Ne1,Ni1,  0-7
Ix1e = mxGetPr(prhs[0]);
Nt = mxGetM(prhs[0]);
m1 = mxGetN(prhs[0]);
if(Nt==1 && m1!=1){
    Nt=m1;
    m1=1;
}
if(Nt==1 || m1!=1)
    mexErrMsgTxt("Ix1e should be Ntx1 or 1xNt");

Ix2e = mxGetPr(prhs[1]);
m1 = mxGetM(prhs[1]);
m2 = mxGetN(prhs[1]);
if(!((m1==Nt && m2==1)||(m1==1 && m2==Nt))){
   mexPrintf("\n%d %d %d\n",m1,m2,Nt);    
   mexErrMsgTxt("Ix2e should be Ntx1 or 1xNt");
}

Ix1i = mxGetPr(prhs[2]);
m1 = mxGetM(prhs[2]);
m2 = mxGetN(prhs[2]);
if(!((m1==Nt && m2==1)||(m1==1 && m2==Nt)))
   mexErrMsgTxt("Ix1i should be Ntx1 or 1xNt");

Ix2i = mxGetPr(prhs[3]);
m1 = mxGetM(prhs[3]);
m2 = mxGetN(prhs[3]);
if(!((m1==Nt && m2==1)||(m1==1 && m2==Nt)))
   mexErrMsgTxt("Ix2i should be Ntx1 or 1xNt");


Ne=(int)mxGetScalar(prhs[4]);
Ni=(int)mxGetScalar(prhs[5]);

Ne1=(int)mxGetScalar(prhs[6]);
Ni1=(int)mxGetScalar(prhs[7]);


// J,Cm,gl,vl,DeltaT,vT,tref,vth,vre,vlb, 8-17
J = mxGetPr(prhs[8]);  // Connection weights
total_num_synapses = mxGetM(prhs[8]);
m2 = mxGetN(prhs[8]);
if(total_num_synapses==1 && m2!=1)
    total_num_synapses=m2;
C=mxGetPr(prhs[9]);
m1 = mxGetN(prhs[9]);
m2 = mxGetM(prhs[9]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
gl=mxGetPr(prhs[10]);
m1 = mxGetN(prhs[10]);
m2 = mxGetM(prhs[10]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
Vleak=mxGetPr(prhs[11]);
m1 = mxGetN(prhs[11]);
m2 = mxGetM(prhs[11]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
DeltaT=mxGetPr(prhs[12]);
m1 = mxGetN(prhs[12]);
m2 = mxGetM(prhs[12]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
VT=mxGetPr(prhs[13]);
m1 = mxGetN(prhs[13]);
m2 = mxGetM(prhs[13]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
tref=mxGetPr(prhs[14]);
m1 = mxGetN(prhs[14]);
m2 = mxGetM(prhs[14]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
Vth=mxGetPr(prhs[15]);
m1 = mxGetN(prhs[15]);
m2 = mxGetM(prhs[15]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
Vre=mxGetPr(prhs[16]);
m1 = mxGetN(prhs[16]);
m2 = mxGetM(prhs[16]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
Vlb=mxGetPr(prhs[17]);
m1 = mxGetN(prhs[17]);
m2 = mxGetM(prhs[17]);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");

// tausyn,V0,T,dt,maxns,Irecord 18-23
tausyn=mxGetPr(prhs[18]); // synaptic decay parameter (varies between synapses)
total_types_synapses = mxGetM(prhs[18]);
m2 = mxGetN(prhs[18]);
if(total_types_synapses==1 && m2!=1)
    total_types_synapses=m2;

s_en=mxGetScalar(prhs[29]); // synaptic decay parameter (varies between synapses)
// mexPrintf("Strength of ensheathment: %f\n", s_en);
// mexEvalString("drawnow");


v0 = mxGetPr(prhs[19]);       // initial conditions for voltage
N = mxGetM(prhs[19]);
m2 = mxGetN(prhs[19]);
if(N==1 && m2!=1)
    N=m2;
T = mxGetScalar(prhs[20]);   // total time of simulation
dt = mxGetScalar(prhs[21]);  // time step

maxns = ((int)mxGetScalar(prhs[22])); // max number of spikes before quitting

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
if(Nrecord==1 && m1!=1){
    Nrecord=m1;
    m1=1;
}
if(Nrecord==1 || m1!=1)
    mexErrMsgTxt("Irecord should be Nrecordx1 or 1xNrecord");

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

/******
 * Finished importing variables.
 *******/


N=Ne+Ni;

Nt=(int)round(T/dt);

/******
 * Now allocate new variables.
 * This is also mex specific.  Use malloc in C, etc.
 *****/

/* Allocate output vector */
plhs[0] = mxCreateDoubleMatrix(2, maxns, mxREAL);
s=mxGetPr(plhs[0]);

/* Allocate output vector */
plhs[1] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
etare=mxGetPr(plhs[1]);
/* Allocate output vector */
plhs[2] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
etari=mxGetPr(plhs[2]);


/* Allocate output vector */
plhs[3] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
vr=mxGetPr(plhs[3]);

/* Allocate membrane potential */
v=mxMalloc(N*sizeof(double));

// four types of inputs
eta_vec_next_e0=mxMalloc(N*sizeof(double));
eta_vec_next_e1=mxMalloc(N*sizeof(double));
eta_vec_next_i0=mxMalloc(N*sizeof(double));
eta_vec_next_i1=mxMalloc(N*sizeof(double));

eta_adjusted_e0=mxMalloc(N*sizeof(double));
eta_adjusted_e1=mxMalloc(N*sizeof(double));
eta_adjusted_i0=mxMalloc(N*sizeof(double));
eta_adjusted_i1=mxMalloc(N*sizeof(double));

refstate=mxMalloc(N*sizeof(int));

/*****
 * Finished allocating variables
 ****/

/* Inititalization*/
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
}

for(jj=0;jj<Nrecord;jj++){
    /* Find index into local variables */
    
    /* Find index into local variables */
	j=Irecord[jj];
	//j=(int)round(Irecord[jj]-1);
                
	if(j>=N || j<0)
		mexErrMsgTxt("Irecord contains out of bounds indices.");
  
	etare[jj+Nrecord*0]=eta_adjusted_e0[j]+eta_adjusted_e1[j];
	etari[jj+Nrecord*0]=eta_adjusted_i0[j]+eta_adjusted_i1[j];
	vr[jj+Nrecord*0]=v[j];
}


Ntref[0]=(int)round(tref[0]/dt);
Ntref[1]=(int)round(tref[1]/dt);

/* Number of x spikes in whole network at each time bin */
/*
nex=(int)round(((double)Ne)*rxe*dt);
nix=(int)round(((double)Ni)*rxi*dt);
if(nex<100 || nix<100)
    mexWarnMsgTxt("nex or nix is small.  Bad approx to Poisson for input spikes.");
*/

/* Initialize number of spikes */
ns=0;


flag=0;

for(i=0;i<Nt && ns<maxns;i++){
	
	
	// loop through each neuron
	for(j=0; j < N; j++){
	     				
		// Uses Euler's method to approximate the exponential decay     		 		
		eta_adjusted_e0[j] = eta_adjusted_e0[j]-eta_adjusted_e0[j]*(dt/tausyn[0]);
		eta_adjusted_e1[j] = eta_adjusted_e1[j]-eta_adjusted_e1[j]*(dt/tausyn[1]);
		eta_adjusted_i0[j] = eta_adjusted_i0[j]-eta_adjusted_i0[j]*(dt/tausyn[2]);
		eta_adjusted_i1[j] = eta_adjusted_i1[j]-eta_adjusted_i1[j]*(dt/tausyn[3]);

		// code for the excitatory cells
		if(j < Ne){
			
			// if the cell is not in the refractory state
			if(refstate[j]<=0){
				//v[j]=v[j]+(eta[j]-gl[0]*(v[j]-Vleak[0])
					//+gl[0]*DeltaT[0]*EXP((v[j]-VT[0])/DeltaT[0]))*dt/C[0];
				
				v[j]=v[j]+(eta_adjusted_e0[j]+eta_adjusted_e1[j]+eta_adjusted_i0[j]+eta_adjusted_i1[j]
					-gl[0]*(v[j]-Vleak[0])
					+gl[0]*DeltaT[0]*exp((v[j]-VT[0])/DeltaT[0]))*dt/C[0];
                
				// Add the correct noisy current for each subpopulation
				if(j < Ne1){
					v[j]=v[j]+Ix1e[i]*dt/C[0];
				}else{
					v[j]=v[j]+Ix2e[i]*dt/C[0];
				}
                
				v[j]=fmax(v[j],Vlb[1]);
				
			}else{
				// Hold the voltage at threshold for the refractory time
			    // and drop it to Vre at the last time point
				if(refstate[j]>1){
					v[j]=Vth[0]; /*-=(Vth[0]-Vre[0])/((double)Ntref[0]);*/
				}else{
					v[j]=Vre[0];
				}
							
				refstate[j]= refstate[j]-1;
			}
			
			
            // If a spike occurs, the neuron enters the refactory state for Ntref time steps
            if(v[j]>=Vth[0] && refstate[j]<=0 && ns<maxns){
                refstate[j]=Ntref[0];
                v[j]=Vth[0];       // reset membrane potential
                s[0+2*ns] = i*dt;  // spike time
                s[1+2*ns] = j+1;   // neuron index  
                ns=ns+1;           // update total number of spikes
                
                // For each postsynaptic target   
                id_current=id_start_firing[j];
                while(firing_neuron[id_current] ==j){
					
					// Checks whether the synapse in ensheathed
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
		}else{ //code for the inhibitory cells
            // if the cell is not refractory
            if(refstate[j]<=0){
                //v[j]=v[j]+(eta[j]-gl[1]*(v[j]-Vleak[1])
				//	+gl[1]*DeltaT[1]*EXP((v[j]-VT[1])/DeltaT[1]))*dt/C[1];
				
                v[j]=v[j]+(eta_adjusted_e0[j]+eta_adjusted_e1[j]+eta_adjusted_i0[j]+eta_adjusted_i1[j]
					-gl[1]*(v[j]-Vleak[1])
					+gl[1]*DeltaT[1]*exp((v[j]-VT[1])/DeltaT[1]))*dt/C[1];
                
                if(j< (Ne+Ni1)){
					v[j]=v[j]+Ix1i[i]*dt/C[1];

                }else{
					v[j]=v[j]+Ix2i[i]*dt/C[1];
                	
                }
				v[j]=fmax(v[j],Vlb[1]);
			}else{
                if(refstate[j]>1){
                	v[j]=Vth[1];
                }else{
                	v[j]=Vre[1];
                }
                refstate[j]=refstate[j]-1;
			}
			
            // If a spike occurs 
            if(v[j]>=Vth[1] && refstate[j]<=0 && ns<maxns){
				
                refstate[j]=Ntref[1];
                v[j]=Vth[1];      // reset membrane potential 
                s[0+2*ns] = i*dt;
                s[1+2*ns] = j+1;
                
                ns=ns+1; // update total number of spikes 
                
                // For each postsynaptic target             
                id_current=id_start_firing[j];
                while(firing_neuron[id_current] ==j){
					
					// Checks whether the synapse in ensheathed
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
	
    // Store recorded variables
    for(jj=0;jj<Nrecord;jj++){
        // Find index into local variables 
        j= Irecord[jj];
    
        if(j>=N || j<0){
        	mexErrMsgTxt("Bad index in Irecord.");
        }
            
   
        etare[jj+Nrecord*i]=eta_adjusted_e0[j]+eta_adjusted_e1[j];
		etari[jj+Nrecord*i]=eta_adjusted_i0[j]+eta_adjusted_i1[j];
        vr[jj+Nrecord*i]=v[j];    
    	
    }
	
    /* Propagate spikes */
	for(j=0;j<N;j++){
		
		eta_adjusted_e0[j] = eta_adjusted_e0[j]+eta_vec_next_e0[j]/tausyn[0];
		eta_adjusted_e1[j] = eta_adjusted_e1[j]+eta_vec_next_e1[j]/tausyn[1];
		eta_adjusted_i0[j] = eta_adjusted_i0[j]+eta_vec_next_i0[j]/tausyn[2];
		eta_adjusted_i1[j] = eta_adjusted_i1[j]+eta_vec_next_i1[j]/tausyn[3];
				
		eta_vec_next_e0[j] = 0;
		eta_vec_next_e1[j] = 0;
		eta_vec_next_i0[j] = 0;
		eta_vec_next_i1[j] = 0;
	}
		
	//prints out the progress at 5% intervals
	// frac_complete = i/((double) Nt);
	// if(i%11000==0 || i==0){
	// 	mexPrintf("Current timestep: %d; Percent done: %f \n", i,frac_complete);
	// 	mexEvalString("drawnow");
	// }
			       
}

/* Issue a warning if max number of spikes reached */
if(ns>=maxns)
   mexWarnMsgTxt("Maximum number of spikes reached, simulation terminated.");

/* Free allocated memory */
mxFree(v);
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




