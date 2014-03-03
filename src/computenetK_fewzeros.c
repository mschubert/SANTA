#include <R.h>
#include <Rdefines.h>

SEXP computenetK_fewzeros(SEXP B, SEXP vw, SEXP n_verticesR, SEXP maxBR)
{
    // B: a vector version of B (by columns)
    // vw: the vertex weights 
    // n_verticesR: the number of vertices in g
    // maxBR: the maximum value of B
    
    // setup 
    SEXP netK;
    int n_vertices = asInteger(n_verticesR);
    int maxB = asInteger(maxBR);
    int *p_B = INTEGER(PROTECT(B)); // create pointer to the B vector
    double *p_vw = REAL(PROTECT(vw)); // create pointer to the vertex weights 
    int i;
    double mean_vw, normalization_c; // mean_vw: mean vertex weight. normalization_c: the normalization character. 
    double *p_nodeS, *p_netK, *p_adjusted_vw; // p_adjusted_vw: pointer to the adjusted vertex weights 
    
    // compute the mean vertex weight 
    mean_vw = 0;
    for (i = 0; i < n_vertices; i++) mean_vw += p_vw[i];
    mean_vw = mean_vw / n_vertices;
    
    // adjust the vertex weights, by subtracting the mean vertex weight from each
    p_adjusted_vw = malloc(sizeof(double) * n_vertices);
    for (i = 0; i < n_vertices; i++) p_adjusted_vw[i] = p_vw[i] - mean_vw;

    // create nodeS and cycle through B, adding adjusted weights to nodeS
    p_nodeS = malloc(sizeof(double) * n_vertices * maxB);
    memset(p_nodeS, 0, sizeof(double) * n_vertices * maxB);
    for (i = 0; i < n_vertices * n_vertices; i++) p_nodeS[(p_B[i] - 1) * n_vertices + (i % n_vertices)] += p_adjusted_vw[i / n_vertices];
    
    // compute cumulative weights
    for (i = n_vertices; i < n_vertices * maxB; i++) p_nodeS[i] += p_nodeS[i - n_vertices];
    
    // compute netK, by multiplying rows by weights and summing columns
    PROTECT(netK = allocVector(REALSXP, maxB));
    memset(REAL(netK), 0, maxB * sizeof(double));
    p_netK = REAL(netK); 
    for (i = 0; i < n_vertices * maxB; i++) p_netK[i / n_vertices] += p_nodeS[i] * p_vw[i % n_vertices];
    
    // normalize netK 
    normalization_c = 2 / (mean_vw * mean_vw * n_vertices * n_vertices);
    for (i = 0; i < maxB; i++) p_netK[i] = p_netK[i] * normalization_c;
    
    free(p_adjusted_vw);
    free(p_nodeS);
    UNPROTECT(3);
    return(netK);
}
