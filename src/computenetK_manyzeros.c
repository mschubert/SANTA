#include <R.h>
#include <Rdefines.h>

SEXP computenetK_manyzeros(SEXP B, SEXP seed_vw, SEXP query_vw, SEXP n_verticesR, SEXP maxBR)
{
    // B: a vector version of B (by columns)
    // seed_vw: the seed vertex weights (the same as query_vw if Knet)
    // query_vw: the query vertex weights (the same as seed_vw if Knet)
    // n_verticesR: the number of vertices in g
    // maxBR: the maximum value of B
    
    // setup 
    SEXP netK;
    int n_vertices = asInteger(n_verticesR);
    int maxB = asInteger(maxBR);
    int *p_B = INTEGER(PROTECT(B)); // create pointer to the B vector
    double *p_seed_vw = REAL(PROTECT(seed_vw)); // create pointer to the seed vertex weights 
    double *p_query_vw = REAL(PROTECT(query_vw)); // create pointer to the seed vertex weights 
   
    int i, j, n_nonzero_seed, B_value;
    int *p_which_seed_nonzero; // p_which_seed_nonzero: pointer to the vector of non-zero vertex weights
    double mean_seed_vw, mean_query_vw, normalization_c; // mean_seed_vw: mean seed vertex weight. mean_query_vw: mean query vertex weight. normalization_c: the normalization character. 
    double *p_nodeS, *p_netK, *p_adjusted_query_vw; // p_adjusted_query_vw: pointer to the adjusted query vertex weights 
    
    // compute the number of non zero (seed) vertex weights
    n_nonzero_seed = 0;
    for (i = 0; i < n_vertices; i++) if (p_seed_vw[i] != 0) n_nonzero_seed++;

    // idenify those nodes with non zero seed weights
    p_which_seed_nonzero = R_alloc(n_nonzero_seed, sizeof(int));
    j = 0;
    for (i = 0; i < n_vertices; i++) {
        if (p_seed_vw[i] != 0) {
            p_which_seed_nonzero[j] = i;
            j++; 
        }
    }
    
    // compute the mean (seed) vertex weight
    mean_seed_vw = 0;
    for (i = 0; i < n_vertices; i++) mean_seed_vw += p_seed_vw[i];
    mean_seed_vw = mean_seed_vw / n_vertices;
    
    // compute the mean (query) vertex weight 
    mean_query_vw = 0;
    for (i = 0; i < n_vertices; i++) mean_query_vw += p_query_vw[i];
    mean_query_vw = mean_query_vw / n_vertices;
    
    // adjust the (query) vertex weights, by subtracting the mean vertex weight from each
    p_adjusted_query_vw = malloc(sizeof(double) * n_vertices);
    memset(p_adjusted_query_vw, 0.0, sizeof(double) * n_vertices);
    for (i = 0; i < n_vertices; i++) p_adjusted_query_vw[i] = p_query_vw[i] - mean_query_vw;

    // create nodeS and cycle through B, adding adjusted weights to nodeS
    p_nodeS = malloc(sizeof(double) * n_nonzero_seed * maxB);
    memset(p_nodeS, 0, sizeof(double) * n_nonzero_seed * maxB);
    for (i = 0; i < n_nonzero_seed * n_vertices; i++) {  
        B_value = p_B[p_which_seed_nonzero[i % n_nonzero_seed] + (i / n_nonzero_seed) * n_vertices] - 1;
        p_nodeS[B_value * n_nonzero_seed + (i % n_nonzero_seed)] += p_adjusted_query_vw[i / n_nonzero_seed];
    }
           
    // compute cumulative weights
    for (i = n_nonzero_seed; i < n_nonzero_seed * maxB; i++) p_nodeS[i] += p_nodeS[i - n_nonzero_seed];
    
    // compute netK, by multiplying rows by weights and summing columns
    PROTECT(netK = allocVector(REALSXP, maxB));
    memset(REAL(netK), 0, maxB * sizeof(double));
    p_netK = REAL(netK);
    for (i = 0; i < n_nonzero_seed * maxB; i++) p_netK[i / n_nonzero_seed] += p_nodeS[i] * p_seed_vw[p_which_seed_nonzero[i % n_nonzero_seed]];
        
    // normalize netK 
    normalization_c = 2 / (mean_seed_vw * mean_query_vw * n_vertices * n_vertices);
    for (i = 0; i < maxB; i++) p_netK[i] = p_netK[i] * normalization_c;
    
    free(p_adjusted_query_vw);
    free(p_nodeS);
    UNPROTECT(4);
    return(netK);
}
