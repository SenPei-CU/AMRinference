  //Integration master equations of AMRO transmission.
/*********************************************************************
 * masterequation_p.c
 * Keep in mind:
 * <> Use 0-based indexing as always in C or C++
 * <> Indexing is column-based as in Matlab (not row-based as in C)
 * <> Use linear indexing.  [x*dimy+y] instead of [x][y]
 * Adapted from the code by Shawn Lankton (http://www.shawnlankton.com/2008/03/getting-started-with-mex-a-short-tutorial/)
 ********************************************************************/
#include <matrix.h>
#include <mex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

void mexFunction(int nlmxhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//declare variables
    mxArray *mxx, *mxnl, *mxpart, *mxpatlist, *mxdegree, *mxpara, *mxalpha, *mxreadmittedpat, *mxnewx;

    const mwSize *dims;
    double *x, *nl, *part, *patlist, *degree, *para, *alpha, *readmittedpat, *newx;
    int dimx, dimy, Nmax, N, num_ens, num_para, num_readmittedpat;
    //Nmax: number of total patients, N: number of currently hospitalized patients

//associate inputs
    mxx = mxDuplicateArray(prhs[0]);
    mxnl = mxDuplicateArray(prhs[1]);
    mxpart = mxDuplicateArray(prhs[2]);
    mxpatlist = mxDuplicateArray(prhs[3]);
    mxdegree = mxDuplicateArray(prhs[4]);
    mxpara = mxDuplicateArray(prhs[5]);
    mxalpha = mxDuplicateArray(prhs[6]);
    mxreadmittedpat = mxDuplicateArray(prhs[7]);
    

//figure out dimensions
    dims = mxGetDimensions(prhs[0]);//x
    dimy = (int)dims[0]; dimx = (int)dims[1]; 
    Nmax = dimy;//number of total patients
    num_ens = dimx;
    
    dims = mxGetDimensions(prhs[3]);//patlist
    N = (int)dims[0];//number of currently hospitalized patients
    
    dims = mxGetDimensions(prhs[5]);//para
    num_para = (int)dims[0];//number of parameters
    
    dims = mxGetDimensions(prhs[7]);//readmittedpat
    num_readmittedpat = (int)dims[0];//number of readmittedpat

//associate outputs
    mxnewx = plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    

//associate pointers
    x = mxGetPr(mxx);
    nl = mxGetPr(mxnl);
    part = mxGetPr(mxpart);
    patlist = mxGetPr(mxpatlist);
    degree = mxGetPr(mxdegree);
    para = mxGetPr(mxpara);
    alpha = mxGetPr(mxalpha);
    readmittedpat = mxGetPr(mxreadmittedpat);
    newx = mxGetPr(mxnewx);

//do something
    int i, j, node, nei, t, k;
    double v, gamma, beta, alphai, temp;
    
    
    for (k = 0; k < num_ens; k++){
        
        //parameters
        beta = para[Nmax+k*num_para]; 
//         printf("%f\n",beta);

        //assign states for patients who first appear in hospital
        for (i = 0; i < N; i++) {
            node = patlist[i]-1;//index from 0
            gamma = para[node+k*num_para]; 
            if (x[node+k*Nmax] == -1) {//first appear
                x[node+k*Nmax] = gamma;//colonized
            }
        }

        //copy x to newx
        for (i = 0; i < Nmax; i++)
            newx[i+k*Nmax] = x[i+k*Nmax];
        

    //     printf("1\n");

        //update nodes in the network  
        for (i = 0; i < N; i++) {
            node = (int)(patlist[i]-1);
            alphai = alpha[node+k*Nmax];//decolonization rate
            newx[node+k*Nmax] = x[node+k*Nmax] - alphai*x[node+k*Nmax];//decolonization
            
            temp = 1;//probability of no transmission
            //neighbors
            for (j = part[node]-1; j < part[node+1]-1; j++){
                nei = (int)(nl[j]-1);
//                 newx[node+k*Nmax] += beta*(1-x[node+k*Nmax])*x[nei+k*Nmax];
                temp = temp*(1-beta*(1-x[node+k*Nmax])*x[nei+k*Nmax]);

            }
            newx[node+k*Nmax] = newx[node+k*Nmax] + 1 - temp;
            //check bound
            if (newx[node+k*Nmax]>1)
                newx[node+k*Nmax]=1;
            if (newx[node+k*Nmax]<0)
                newx[node+k*Nmax]=0;
        }

    //     printf("2\n");

      //update nodes ever appear that outside hospitals
        for (i = 0; i < num_readmittedpat; i++) {
            node = readmittedpat[i]-1;
            alphai = alpha[node+k*Nmax];
            newx[node+k*Nmax] = x[node+k*Nmax] - alphai*x[node+k*Nmax];
            //check bound
            if (newx[node+k*Nmax]>1)
                newx[node+k*Nmax]=1;
            if (newx[node+k*Nmax]<0)
                newx[node+k*Nmax]=0;
        }

    //     printf("3\n");
        //copy newx to x
        for (i = 0; i < Nmax; i++)
            x[i+k*Nmax] = newx[i+k*Nmax];

    }

    return;
}
