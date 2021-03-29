#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mex.h"



void infoDistance(int M, int L, uint32_T mtx[], double CovIJ[], uint32_T sub[], uint32_T NN[]) {
#pragma omp parallel for
  
  int PDSize = 21;
  int N = NN[0];
  double zerocase = 0.;
  
  // Pre-calculate all frequencies
  int Ms[L][PDSize]; // Marginal probabilities

  int treepop[N];
  double Hpos[L];
  for (int i = 0; i < L; i++) {
    Hpos[i] = 0.;

    for (int n = 0; n < PDSize; n++) {
      Ms[i][n] = 0;
    }
  }
  for (int m = 0; m < N; m++) {
    treepop[m] = 0; //pre-allocate
  }
   
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      CovIJ[i+L*j] = (double) zerocase;
    }
  }

  // Counting
  for (int i = 0; i < M; i++) {
    treepop[sub[i]] += 1; // count
  }
  
  for (int i = 0; i < L; i++) {
    for (int k = 0; k < M; k++) {
      Ms[i][mtx[k+M*i]] += 1;
    }
  }

  //  Calculate marginal entropies

  for (int i = 0; i < L; i++) {
    for (int n = 0; n < PDSize; n++) {
      if (Ms[i][n] > 0) {
        double curp = (double) Ms[i][n]/M;
        Hpos[i] += -curp*log(curp)/log(2); // Calculate Hi
      }
    }
  }
  
  // Now calculate covariance using information theoretic methods
  for (int i = 0; i < L; ++i) {
    double Hi = Hpos[i];
    if (Hi > 0.) {

      int MPis[PDSize][N]; // Marginal probabilities for all subtrees
      for (int n = 0; n < PDSize; n++) {
        for (int k = 0; k < N; k++) {
          MPis[n][k] = 0;
        }
      }
      for (int k = 0; k < M; k++) {
        MPis[mtx[k+M*i]][sub[k]] += 1; // count
      }
    
      
      for (int j = 0; j <= i; ++j) {
        double Hj = Hpos[j];
        if (Hj > 0.) {
          int MPjs[PDSize][N]; // Marginal probabilities for all subtrees
          for (int n = 0; n < PDSize; n++) {
            for (int k = 0; k < N; k++) {
              MPjs[n][k] = 0;
            }
          }
          for (int k = 0; k < M; k++) {
            MPjs[mtx[k+M*j]][sub[k]] += 1; // count
          }
          
          
          double JH = 0., JHnull = 0.;
          
          double MPs[PDSize][PDSize]; // Marginal probabilities
          for (int n = 0; n < PDSize; n++) {
            for (int m = 0; m < PDSize; m++) {
              MPs[n][m] = 0.;
            }
          }
        
          double accessed;
          for (int s = 0; s < M; s++) {
            MPs[mtx[s+M*i]][mtx[s+M*j]] += 1.;
          }
        
          //Sum up the information where the JPD is nonzero
        
          for (int n = 0; n < PDSize; n++) {
            if (Ms[i][n] > 0) {
              for (int m = 0; m < PDSize; m++) {
                if (Ms[j][m] > 0) {
                  if (MPs[n][m] > 0.) {
                    double JProb = MPs[n][m]/M;
                    JH += -JProb*log(JProb)/log(2);
                  }
                
                  double JProbnull = 0.;
                
                  for (int k = 0; k < N; k++) {
                    JProbnull += (double) MPis[n][k]*MPjs[m][k]/treepop[k];
                  }
                
                  if (JProbnull > 0.) {
                    JHnull += -JProbnull*log(JProbnull/M)/log(2)/M;
                  }
                }        
              }
            }
          }
            
          if (JH == 0.)
            CovIJ[j+L*i] = CovIJ[i+L*j] = (double) zerocase;
          else if (JHnull == 0.)
            CovIJ[j+L*i] = CovIJ[i+L*j] = (double) zerocase;
          else
            CovIJ[j+L*i] = CovIJ[i+L*j] = (double) (2*JHnull - Hi - Hj)/JHnull-(2*JH - Hi - Hj)/JH;
            
        }
      }
    }
  }
  

}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int nrows;
  int ncols;
  uint32_T *inMatrix;               /* MxN input matrix */
  uint32_T *subtrees;              /* categorical variable */
  uint32_T *numtrees; /* number of unique trees in categorical variable */
  double *outMatrix;              /* output matrix */


  /* check for proper number of arguments */
  if(nrhs!=3) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Three inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    
    /* make sure the first input argument is type int */
    if( !mxIsUint32(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notInt","Input matrix must be type uint32.");
    }
    /* create a pointer to the real data in the input matrix  */
    inMatrix = (uint32_T *)mxGetData(prhs[0]);

    /* make sure the second input argument is type int */
    if( !mxIsUint32(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notInt","Categorical vector must be type uint32.");
    }
    /* create a pointer to the real data in the categorical matrix */
    subtrees = (uint32_T *)mxGetData(prhs[1]);

    /* make sure the third input argument is type int */
    if( !mxIsUint32(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notInt","Categorical vector must be type uint32.");
    }
    /* create a pointer to the number of data */
    numtrees = (uint32_T *) mxGetData(prhs[2]);
    

    /* get dimensions of the input matrix */
    nrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,(int)ncols*ncols,mxREAL);

    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    infoDistance(nrows,ncols,inMatrix,outMatrix,subtrees,numtrees);
}


