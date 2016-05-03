#include <R.h>
#include <Rdefines.h>

SEXP computenodeAUK(SEXP Bv, SEXP vw, SEXP nverticesR, SEXP maxBR)
{
    // Bv: a vector version of B (by columns)
    // vw: the vertex weights 
    // nverticesR: the number of vertices in g
    // maxBR: the maximum value of B
    // could I speed up code, by deleting variables as I go along?
    
    // setup 
    SEXP nodeS, nodeAUK, avw;
    int nvertices = asInteger(PROTECT(nverticesR));
    int maxB = asInteger(PROTECT(maxBR));
    int i;
    double normc, mvw; // normc: the normalization character. mvw: mean vertex weight
    int *pBv = INTEGER(PROTECT(Bv)); // create pointer to the B vector
    double *pvw = REAL(PROTECT(vw)); // create pointer to the vertex weights 
    double *pnodeS, *pnodeAUK, *pavw; // pavw: pointer to the adjusted vertex weights 
    
    // compute the mean vertex weight
    mvw = 0;
    for (i = 0; i < nvertices; i++) mvw += pvw[i];
    mvw = mvw / nvertices;
    
    // adjust the vertex weights, by subtracting the mean vertex weight from each
    PROTECT(avw = allocVector(REALSXP, nvertices));
    memset(REAL(avw), 0, nvertices * sizeof(normc));
    pavw = REAL(avw); 
    for (i = 0; i < nvertices; i++) pavw[i] = pvw[i] - mvw;

    // create nodeS and set all values to 0
    PROTECT(nodeS = allocVector(REALSXP, nvertices * maxB));
    memset(REAL(nodeS), 0, nvertices * maxB * sizeof(normc));
    pnodeS = REAL(nodeS); 
        
    // cycle through B, adding as required
    for (i = 0; i < nvertices * nvertices; i++) pnodeS[(pBv[i] - 1) * nvertices + (i % nvertices)] += pavw[i / nvertices];
    
    // compute cumulative weights (nodeS)
    for (i = nvertices; i < nvertices * maxB; i++) pnodeS[i] += pnodeS[i - nvertices];
    
    // create nodeAUK
    PROTECT(nodeAUK = allocVector(REALSXP, nvertices));
    memset(REAL(nodeAUK), 0, nvertices * sizeof(normc));
    pnodeAUK = REAL(nodeAUK); 
    
    // sum columns then divide by the number of columns (nodeAUK)
    for (i = 0; i < nvertices * maxB; i++) pnodeAUK[i % nvertices] += pnodeS[i];
    for (i = 0; i < nvertices; i++) pnodeAUK[i] = pnodeAUK[i] / nvertices; 
    
    UNPROTECT(7);
    return(nodeAUK);
}
