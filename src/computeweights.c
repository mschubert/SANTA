#include <R.h>
#include <Rdefines.h>

SEXP computeweights(SEXP b, SEXP weight, SEXP Rnvertices, SEXP Rdiam)
{
    int i;
    
    int nvertices = asInteger(Rnvertices);
    int diam = asInteger(Rdiam);
    
    // create pointers to long vectors
    int *pb = INTEGER(PROTECT(b));
    double *pweight = REAL(weight);
    
    // create results vector
    SEXP result;
    PROTECT(result = allocVector(REALSXP, nvertices * diam));
    double *presult = REAL(result);
        
    for (i = 0; i < nvertices * diam; i++) {
        presult[i] = 0;
    }
    
    for (i = 0; i < nvertices * nvertices; i++) {
        presult[(pb[i] - 1) * nvertices + (i % nvertices)] += pweight[i / nvertices];
    }
    
    UNPROTECT(1);
    return(result);
}
