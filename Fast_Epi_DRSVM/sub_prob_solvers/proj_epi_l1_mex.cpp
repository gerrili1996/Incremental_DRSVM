
#include <iostream>
#include <math.h>
#include <matrix.h>
#include <algorithm>
#include <stdlib.h>
#include <mex.h>
 
#define v_in  prhs[0] 
#define s_in prhs[1] 
#define x_out plhs[0] 
#define t_out plhs[1] 
 

void proj_epi_l1(double* x, double* t, double* v, double s, int n)
{
    double* y = (double*)malloc(n * sizeof(double)); // abs(v)
    double* y_upper = (double*)malloc(n * sizeof(double));// y(y>lam)
    double* y_lower = (double*)malloc(n * sizeof(double)); // y(y<lam)
    double a, b, lam,g,sum_y;
    a = -s;
    b = 1;
    sum_y = 0;
    int nsize = n;
    for (int i = 0; i < n; i++)
    {
        y[i] = fabs(v[i]);
        sum_y += y[i];
    }
    if (sum_y <= s){
        memcpy(x, v, n*sizeof(double));
        *t = s; 
        return; 
    }
    while (nsize != 0)
    {
        sum_y = 0; 
        int nk =0;  // length of y_upper
        int nc = 0; // length of y_lower
        lam = y[0]; 
        for (int i = 0; i < nsize; i++) 
        {
            if (y[i] > lam)
            {
                y_upper[nk] = y[i];
                sum_y += y[i];
                nk++; 
            }
        }
        g = a + sum_y - lam * (b +nk);
        if (g < 0)
        {
            a = a + sum_y +lam;
            b = b + nk + 1; 
           for (int i = 0; i < nsize; i++) 
           {
               if (y[i]<lam)
               { 
               y_lower[nc] = y[i];
                nc++;
               }
            }
            memcpy(y, y_lower, n * sizeof(double));
            nsize = nc;
        }
        else if (g > 0)
        {
            memcpy(y, y_upper, n * sizeof(double));
            nsize = nk;
        }
        else break;      
    }

    lam = a / b;
    *t = lam + s;
    for (int i = 0; i < n; i++)
    {
        x[i] = fmax(0, v[i] - lam) - fmax(0, -v[i] - lam);
    }
    free(y_upper);
    free(y_lower);
    free(y);
}

 
/*** interface ****/ 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double p; 
    /* Check Input */
    if(nrhs < 1 || nrhs > 2)
            mexErrMsgTxt("Must have either 1 or 2 input arguments.");
    if(nlhs > 2)
            mexErrMsgTxt("Too many output arguments.");
    if(mxIsComplex(v_in)  ||  mxGetNumberOfDimensions(v_in) != 2 || mxIsSparse(v_in)  ||  !mxIsDouble(v_in))
            mexErrMsgTxt("Sorry! v must be a real vector.");
    
  if(mxIsComplex(s_in)||!mxIsDouble(s_in) || mxGetNumberOfElements(s_in) != 1)
                mexErrMsgTxt("s must be a double scalar.");
 
    /*** Read the input ***/
    double s; 
    double *v; 
    int n; 
    s = mxGetScalar(s_in); 
    v = mxGetPr(v_in); // get pointer to v's data  
    n  = mxGetM(v_in);
    /*** Create the output  ***/
    double *x, *t; 
    x_out = mxCreateDoubleMatrix(n,1,mxREAL);
    x = mxGetPr(x_out);  
    t_out = mxCreateDoubleMatrix(1,1,mxREAL);
    t = mxGetPr(t_out); 
    proj_epi_l1(x,t,v,s,n); 
}
