#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <algorithm>


void sumo(double *im, double *imout, double *lam1, double* lam2, double *A, double *B,double sigmamin, double sigmamax, double d,double lambdamin,double ksfactor, int M, int N)
{
double l1,l2,xr,yr,gksum,sigmadim1,sigmadim2,gk;
int ksmax,sumoi,sumoj;

for (int i=0;i<M;i++)
    for (int j=0;j<N;j++)
    {
       l1 = lam1[i+M*j];
       l2 = lam2[i+M*j];
       sigmadim1 = std::min(sigmamax,(sigmamax-sigmamin)*exp(-d*(l1-lambdamin))+sigmamin);
       sigmadim2 = std::min(sigmamax,(sigmamax-sigmamin)*exp(-d*(l2-lambdamin))+sigmamin);
        
       // Creating kernels
       ksmax=(int)ceil(ksfactor*std::max(sigmadim1,sigmadim2));
       gksum = 0;
       for (int xm=-ksmax;xm<=ksmax;xm++)
            for (int ym=-ksmax;ym<=ksmax;ym++)
            {
                xr = A[i+M*j]*xm + B[i+M*j]*ym;
                yr = A[i+M*j]*ym - B[i+M*j]*xm;        
                gk=exp(-0.5*((xr/sigmadim1)*(xr/sigmadim1)+(yr/sigmadim2)*(yr/sigmadim2)));
                //gk=exp(-0.5*((xr/sigmadim1)*(xr/sigmadim1)+(yr/sigmadim2)*(yr/sigmadim2)))/2/M_PI/sigmadim1/sigmadim2;
                
                sumoi = i+ym;
                sumoj = j+xm;
                //gksum +=  gk;
                    
                if (sumoi>0 && sumoj>0 && sumoi<M && sumoj<N)
                {
                    gksum +=  gk;
                    imout[i+M*j]=imout[i+M*j] + gk*im[sumoi+M*sumoj];
                }    
            }
       imout[i+M*j]=imout[i+M*j]/gksum;
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *imout_m;
  const mwSize *dims;
  double *im, *imout, *A, *B, *lam1, *lam2;
  int  N, M;
  double sigmamin,sigmamax,d,lambdamin,ksfactor;
     
  if (nrhs<10)
      mexErrMsgIdAndTxt("structure_smooth_sum:fewinput","Input should be an image and nine other things.");
  //associate inputs/outputs
 
   im = mxGetPr(prhs[0]);
  lam1 = mxGetPr(prhs[1]);
  lam2 = mxGetPr(prhs[2]);
  A = mxGetPr(prhs[3]);
  B = mxGetPr(prhs[4]);
  sigmamin = mxGetScalar(prhs[5]);
  sigmamax = mxGetScalar(prhs[6]);
  d = mxGetScalar(prhs[7]);
  lambdamin = mxGetScalar(prhs[8]);
  ksfactor = mxGetScalar(prhs[9]);
 
  //figure out dimensions
  dims = mxGetDimensions(prhs[0]);
  M = (int)dims[0];
  N = (int)dims[1]; 

  imout_m = plhs[0] = mxCreateDoubleMatrix(M,N,mxREAL);
  imout = mxGetPr(imout_m);

 
  
  
  sumo(im,imout,lam1,lam2,A,B,sigmamin,sigmamax,d,lambdamin,ksfactor,M,N);
  
}
