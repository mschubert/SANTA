#include <R.h>
#include <Rdefines.h>

SEXP computeweights(SEXP Rb, SEXP Rweight, SEXP Rnvertices, SEXP Rdiam)
{
    int i, nvertices, diam;
    
    int *b;
    b = INTEGER(Rb);
    
    double *weight;
    weight = REAL(Rweight);
    
    Rnvertices = coerceVector(Rnvertices, INTSXP);
    nvertices = INTEGER(Rnvertices)[0];
    
    Rdiam = coerceVector(Rdiam, INTSXP);
    diam = INTEGER(Rdiam)[0];
    
    int tot = nvertices * nvertices;
    SEXP result;
    PROTECT(result = allocVector(REALSXP, nvertices * diam));
    
    for (i = 0; i < nvertices * diam; i++) {
        REAL(result)[i] = 0;
    }
    
    for (i = 0; i < tot; i++) {
        REAL(result)[(b[i] - 1) * nvertices + (i % nvertices)] += weight[i / nvertices];
    }
    
    UNPROTECT(1);
    return result;
}
