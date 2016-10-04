//
//  main.cpp
//  test
//
//  Created by FANGYUAN GUO on 2/21/16.
//  Copyright © 2016 FANGYUAN GUO. All rights reserved.
//
#include <iostream>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <cmath>

double max(double a, double b){
    return (a > b ) ? a : b;
}

double ** do_euler_explicit(int N, int M, float r, double sigma, double K, double S_max, double T)
{
    
    
    double ** result_table;
    
    double v_j,v_j_1,v_j_f,v_begin,v_last;
    double x = log(S_max);
    
    result_table = new double *[M];
    for(int i =0; i<M; i++)
    {
        result_table[i] = new double[N];
        for(int j = 0; j<N; j++)
            result_table[i][j] = 0;
    }
    
    // Build Example Matrix
    gsl_matrix * t_m = gsl_matrix_alloc(N,N);
    gsl_matrix_set_zero(t_m);
    
    double delta_t = T/M;
    double h = 2*log(S_max)/N;
    double beta = r - 0.5 * sigma * sigma;
    
    v_j = 1 - sigma * sigma * delta_t / h / h - r * delta_t;
    v_j_1 = ( 0.5 * sigma * sigma / h / h + 0.5 * beta / h ) * delta_t;
    v_j_f = ( 0.5 * sigma * sigma / h / h - 0.5 * beta / h )* delta_t;
    
    
    v_begin = v_j + v_j_f;
    v_last = v_j + v_j_1;
    
    
    for(int i = 0; i < N; i++){
        if (i == 0) {
            gsl_matrix_set(t_m, i, i, v_begin);
        }
        if (i == N-1){
            gsl_matrix_set(t_m, i, i, v_last);
        }
        else if (i>0 && i<N-1)
        gsl_matrix_set(t_m, i, i, v_j);
    }
    
    
    for(int i = 0; i < N-1; i++)
        gsl_matrix_set(t_m, i+1, i, v_j_f);
    
    for(int i = 0; i < N-1; i++)
        gsl_matrix_set(t_m, i, i+1, v_j_1);
    
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
            std::cout<<gsl_matrix_get(t_m, i, j)<<" ";
        std::cout<<std::endl;
    }
    // Solution Vector
    gsl_vector * s_v_old = gsl_vector_alloc(N);
    gsl_vector * s_v_new = gsl_vector_alloc(N);
    
    gsl_vector_set_zero(s_v_old);
    gsl_vector_set_zero(s_v_new);
    
    //example initial value
    //gsl_vector_set(s_v_old, N/2, 2.5);
//    gsl_vector_set(s_v_old, 0, 0);
//    //gsl_vector_set(s_v_old, N/2 - 1, 1.5);
//    gsl_vector_set(s_v_old, N, 1.5);
//    gsl_vector_set(s_v_old, N/2 + 1, 1.5);
    
    ////////////////////added ///////////////
    for (int ii = 0 ; ii < N; ii++) {
        gsl_vector_set(s_v_old, ii, max ( exp(-x + (ii+1) * h) - K, 0));
    }
    //////////////////////////////////////////
    
    std::cout<<"New\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_new, i)<<" ";
    std::cout<<std::endl;
    
    std::cout<<"Old\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_old, i)<<" ";
    std::cout<<std::endl;
    
    std::cout<<"Starting\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_old, i)<<" ";
    std::cout<<std::endl;
    
    for(int t_j = 0; t_j < M; t_j++)
    {
        gsl_blas_dgemv(CblasNoTrans, 1.0, t_m, s_v_old, 0.0, s_v_new);
        
        
        for (int i = 0; i < N; i++) {
            if ( gsl_vector_get(s_v_new, i) < 0 )
                gsl_vector_set(s_v_new , i , 0);
        }
        gsl_vector_set(s_v_new , N-1 , S_max-K*exp(-r*delta_t*(t_j+1)));
        gsl_vector_set(s_v_new , 0 ,0);
        
        for(int i = 0; i < N; i++)
            std::cout<<gsl_vector_get(s_v_new, i)<<" ";
        std::cout<<std::endl;
        
        memcpy(result_table[t_j], s_v_new->data, s_v_new->size*sizeof(double));
        
        gsl_vector_swap(s_v_new, s_v_old);
    }
    gsl_matrix_free(t_m);
    gsl_vector_free(s_v_old);
    gsl_vector_free(s_v_new);
    
    
    return result_table;
};

double ** do_euler_implicitet(int N, int M, float r, double sigma, double K, double S_max, double T)
{

    
    double ** result_table;
    
    double v_j,v_j_1,v_j_f,v_begin,v_last;
    double x = log(S_max);
    
    result_table = new double *[M];
    for(int i =0; i<M; i++)
    {
        result_table[i] = new double[N];
        for(int j = 0; j<N; j++)
            result_table[i][j] = 0;
    }
    

    gsl_matrix * t_m = gsl_matrix_alloc(N,N);
    gsl_matrix_set_zero(t_m);
    
    
    double delta_t = T/M;
    double h = 2*log(S_max)/N;
    double beta = r - 0.5 * sigma * sigma;
    
    v_j = 1 + sigma * sigma * delta_t / h / h + r * delta_t;
    v_j_1 = ( - 0.5 * sigma * sigma / h / h - 0.5 * beta / h ) * delta_t;
    v_j_f = ( - 0.5 * sigma * sigma / h / h + 0.5 * beta / h ) * delta_t;
    
    v_begin = v_j + v_j_f;
    v_last = v_j + v_j_1;
    
    for(int i = 0; i < N; i++){
        if (i == 0) {
            gsl_matrix_set(t_m, i, i, v_begin);
        }
        if (i == N-1){
            gsl_matrix_set(t_m, i, i, v_last);
        }
        else if (i>0 && i<N-1)
            gsl_matrix_set(t_m, i, i, v_j);
    }
//    for(int i = 0; i < N; i++)
//        gsl_matrix_set(t_m, i, i, v_j);
    
    for(int i = 0; i < N-1; i++)
        gsl_matrix_set(t_m, i+1, i, v_j_f);
    
    for(int i = 0; i < N-1; i++)
        gsl_matrix_set(t_m, i, i+1, v_j_1);
    
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
            std::cout<<gsl_matrix_get(t_m, i, j)<<" ";
        std::cout<<std::endl;
    }
    // Solution Vector
    gsl_vector * s_v_old = gsl_vector_alloc(N);
    gsl_vector * s_v_new = gsl_vector_alloc(N);
    
    gsl_vector_set_zero(s_v_old);
    gsl_vector_set_zero(s_v_new);
    
//    gsl_vector_set(s_v_old, N/2, 2.5);
//    gsl_vector_set(s_v_old, N/2 - 1, 1.5);
//    gsl_vector_set(s_v_old, N/2 + 1, 1.5);
    
    ////////////////////added ///////////////
    for (int ii = 0 ; ii < N; ii++) {
        gsl_vector_set(s_v_old, ii, max ( exp(-x + (ii+1) * h) - K, 0));
    }
    //////////////////////////////////////////
    
    std::cout<<"New\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_new, i)<<" ";
    std::cout<<std::endl;
    
    std::cout<<"Old\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_old, i)<<" ";
    std::cout<<std::endl;
    
    // printing out the initial condition
    std::cout<<"Starting\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_old, i)<<" ";
    std::cout<<std::endl;
    
    int s;
    gsl_matrix * inverse_t_m = gsl_matrix_alloc(N,N);
    gsl_matrix_set_zero(inverse_t_m);
    gsl_permutation *p = gsl_permutation_alloc(N);
    gsl_linalg_LU_decomp(t_m, p, &s);
    gsl_linalg_LU_invert (t_m, p, inverse_t_m);
    /////////////////////////////////////////////////
    
    for(int t_j = 0; t_j < M; t_j++)
    {

        gsl_blas_dgemv(CblasNoTrans, 1.0, inverse_t_m, s_v_old, 0.0, s_v_new);
        
        
        for (int i = 0; i < N; i++) {
            if ( gsl_vector_get(s_v_new, i) < 0 )
                gsl_vector_set(s_v_new , i , 0);
        }
        gsl_vector_set(s_v_new , N-1 , S_max-K*exp(-r*delta_t*(t_j+1)));
        gsl_vector_set(s_v_new , 0 ,0);
        
        for(int i = 0; i < N; i++)
            std::cout<<gsl_vector_get(s_v_new, i)<<" ";
        std::cout<<std::endl;
        
        memcpy(result_table[t_j], s_v_new->data, s_v_new->size*sizeof(double));
        
        gsl_vector_swap(s_v_new, s_v_old);
    }
    gsl_matrix_free(t_m);
    gsl_vector_free(s_v_old);
    gsl_vector_free(s_v_new);
    

    return result_table;
};

double ** do_crank_nicolson(int N, int M, float r, double sigma, double K, double S_max, double T)
{
    double v_j_1,v_j_1_1,v_j_f_1,v_j_2,v_j_1_2,v_j_f_2,v_begin_1,v_begin_2,v_last_1,v_last_2;
    double x = log(S_max);
    ////////////////////////////////////////////
    double ** result_table;
    result_table = new double *[M];
    for(int i =0; i<M; i++)
    {
        result_table[i] = new double[N];
        for(int j = 0; j<N; j++)
            result_table[i][j] = 0;
    }
 
    gsl_matrix * t_m_1 = gsl_matrix_alloc(N,N);
    gsl_matrix * t_m_2 = gsl_matrix_alloc(N,N);
    
    gsl_matrix * t_m = gsl_matrix_alloc(N,N);


    gsl_matrix_set_zero(t_m_1);
    gsl_matrix_set_zero(t_m_2);
    gsl_matrix_set_zero(t_m);


    
    double delta_t = T/M;
    double h = 2*log(S_max)/N;
    double beta = r - 0.5 * sigma * sigma;
    
    v_j_1 = 1 / delta_t + 0.5 * sigma * sigma / h / h + 0.5 * r;
    v_j_1_1 = - 0.25 * sigma * sigma / h / h - 0.25 * beta / h;
    v_j_f_1 = - 0.25 * sigma * sigma / h / h + 0.25 * beta / h;
    
    v_j_2 = 1 / delta_t - 0.5 * sigma * sigma / h / h - 0.5 * r;
    v_j_1_2 = 0.25 * sigma * sigma / h / h + 0.25 * beta / h;
    v_j_f_2 = 0.25 * sigma * sigma / h / h - 0.25 * beta / h;

    v_begin_1 = v_j_1 + v_j_f_1;
    v_last_1 = v_j_1 + v_j_1_1;
    
    v_begin_2 = v_j_2 + v_j_f_2;
    v_last_2 = v_j_2 + v_j_1_2;
    
    
    for(int i = 0; i < N; i++){
        if (i == 0) {
            gsl_matrix_set(t_m_1, i, i, v_begin_1);
        }
        if (i == N-1){
            gsl_matrix_set(t_m_1, i, i, v_last_1);
        }
        else if (i>0 && i<N-1)
            gsl_matrix_set(t_m_1, i, i, v_j_1);
    }
//    for(int i = 0; i < N; i++)
//        gsl_matrix_set(t_m_1, i, i, v_j_1);
    
    for(int i = 0; i < N-1; i++)
        gsl_matrix_set(t_m_1, i+1, i, v_j_f_1);
    
    for(int i = 0; i < N-1; i++)
        gsl_matrix_set(t_m_1, i, i+1, v_j_1_1);
    
//    for(int i = 0; i < N; i++)
//        gsl_matrix_set(t_m_2, i, i, v_j_2);

    for(int i = 0; i < N; i++){
        if (i == 0) {
            gsl_matrix_set(t_m_2, i, i, v_begin_2);
        }
        if (i == N-1){
            gsl_matrix_set(t_m_2, i, i, v_last_2);
        }
        else if (i>0 && i<N-1)
            gsl_matrix_set(t_m_2, i, i, v_j_2);
    }
    
    for(int i = 0; i < N-1; i++)
        gsl_matrix_set(t_m_2, i+1, i, v_j_f_2);
    
    for(int i = 0; i < N-1; i++)
        gsl_matrix_set(t_m_2, i, i+1, v_j_1_2);
    
//    for(int i = 0; i < N; i++)
//    {
//        for(int j = 0; j < N; j++)
//            std::cout<<gsl_matrix_get(t_m, i, j)<<" ";
//        std::cout<<std::endl;
//    }

    
    // Solution Vector
    gsl_vector * s_v_old = gsl_vector_alloc(N);
    gsl_vector * s_v_new = gsl_vector_alloc(N);
    
    gsl_vector_set_zero(s_v_old);
    gsl_vector_set_zero(s_v_new);
    
    
    ////////////////////added ///////////////
    for (int ii = 0 ; ii < N; ii++) {
        gsl_vector_set(s_v_old, ii, max ( exp(-x + ii * h) - K, 0));
    }
    //////////////////////////////////////////
    
    std::cout<<"New\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_new, i)<<" ";
    std::cout<<std::endl;
    
    std::cout<<"Old\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_old, i)<<" ";
    std::cout<<std::endl;
    
    // printing out the initial condition
    std::cout<<"Starting\n";
    for(int i = 0; i < N; i++)
        std::cout<<gsl_vector_get(s_v_old, i)<<" ";
    std::cout<<std::endl;
    
    int s;
    gsl_matrix * inverse_t_m_1 = gsl_matrix_alloc(N,N);
    gsl_matrix_set_zero(inverse_t_m_1);
    gsl_permutation *p = gsl_permutation_alloc(N);
    gsl_linalg_LU_decomp(t_m_1, p, &s);
    gsl_linalg_LU_invert (t_m_1, p, inverse_t_m_1);
    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, inverse_t_m_1, t_m_2, 0.0, t_m);

    
    for(int t_j = 0; t_j < M; t_j++)
    {
        // This is what actually does the matrix vector multiplication.

        gsl_blas_dgemv(CblasNoTrans, 1.0, t_m, s_v_old, 0.0, s_v_new);
        
        
        for (int i = 0; i < N; i++) {
            if ( gsl_vector_get(s_v_new, i) < 0 )
                gsl_vector_set(s_v_new , i , 0);
        }
        /////////////////////////////////////////
        /////////////////////////////////////////
        gsl_vector_set(s_v_new , N-1 , S_max-K*exp(-r*delta_t*(t_j+1)));
        gsl_vector_set(s_v_new , 0 ,0);
        /////////////// important ///////////////////////
        
        for(int i = 0; i < N; i++)
            std::cout<<gsl_vector_get(s_v_new, i)<<" ";
        std::cout<<std::endl;
        
        memcpy(result_table[t_j], s_v_new->data, s_v_new->size*sizeof(double));
        
        gsl_vector_swap(s_v_new, s_v_old);
    }
    gsl_matrix_free(t_m);
    gsl_vector_free(s_v_old);
    gsl_vector_free(s_v_new);
    
    return result_table;
};

