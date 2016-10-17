#include <R.h>
#include <Rdefines.h>

SEXP computenetK_manyzeros(SEXP B, SEXP vw, SEXP n_verticesR, SEXP maxBR)
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
   
    int i, j, n_nonzero_vertices, B_value;
    int *p_which_vw_nonzero; // p_which_vw_nonzero: pointer to the vector of non-zero vertex weights
    double mean_vw, normalization_c; // mean_vw: mean vertex weight. normalization_c: the normalization character. 
    double *p_nodeS, *p_netK, *p_adjusted_vw; // p_adjusted_vw: pointer to the adjusted vertex weights 
    
   
    // compute the number of non zero vertex weights
    n_nonzero_vertices = 0;
    for (i = 0; i < n_vertices; i++) if (p_vw[i] != 0) n_nonzero_vertices++;

    // idenify those nodes with non zero weights
    p_which_vw_nonzero = R_alloc(n_nonzero_vertices, sizeof(int));
    j = 0;
    for (i = 0; i < n_vertices; i++) {
        if (p_vw[i] != 0) {
            p_which_vw_nonzero[j] = i;
            j++; 
        }
    }
    
    // compute the mean vertex weight 
    mean_vw = 0;
    for (i = 0; i < n_vertices; i++) mean_vw += p_vw[i];
    mean_vw = mean_vw / n_vertices;
    
    // adjust the vertex weights, by subtracting the mean vertex weight from each
    p_adjusted_vw = malloc(sizeof(double) * n_vertices);
    memset(p_adjusted_vw, 0.0, sizeof(double) * n_vertices);
    for (i = 0; i < n_vertices; i++) p_adjusted_vw[i] = p_vw[i] - mean_vw;

    // create nodeS and cycle through B, adding adjusted weights to nodeS
    p_nodeS = malloc(sizeof(double) * n_nonzero_vertices * maxB);
    memset(p_nodeS, 0, sizeof(double) * n_nonzero_vertices * maxB);
    for (i = 0; i < n_nonzero_vertices * n_vertices; i++) {  
        B_value = p_B[p_which_vw_nonzero[i % n_nonzero_vertices] + (i / n_nonzero_vertices) * n_vertices] - 1;
        p_nodeS[B_value * n_nonzero_vertices + (i % n_nonzero_vertices)] += p_adjusted_vw[i / n_nonzero_vertices];
    }
           
    // compute cumulative weights
    for (i = n_nonzero_vertices; i < n_nonzero_vertices * maxB; i++) p_nodeS[i] += p_nodeS[i - n_nonzero_vertices];
    
    // compute netK, by multiplying rows by weights and summing columns
    PROTECT(netK = allocVector(REALSXP, maxB));
    memset(REAL(netK), 0, maxB * sizeof(double));
    p_netK = REAL(netK);
    for (i = 0; i < n_nonzero_vertices * maxB; i++) p_netK[i / n_nonzero_vertices] += p_nodeS[i] * p_vw[p_which_vw_nonzero[i % n_nonzero_vertices]];
        
    // normalize netK 
    normalization_c = 2 / (mean_vw * mean_vw * n_vertices * n_vertices);
    for (i = 0; i < maxB; i++) p_netK[i] = p_netK[i] * normalization_c;
    
    free(p_adjusted_vw);
    free(p_nodeS);
    UNPROTECT(3);
    return(netK);
}
