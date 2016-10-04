//
//  main.cpp
//  test
//
//  Created by FANGYUAN GUO on 4/7/16.
//  Copyright © 2016 FANGYUAN GUO. All rights reserved.
//
#include <iostream>
#include <fstream>

#include "do_nm.h"

int main(void)
{
    std::ofstream output_file;
    output_file.open("out.csv");
    int N = 30;
    int M = 600;
    double r = 0.02;
    double sigma = .2;
    double K = 20;
    double S_max = 100;
    double T = 20;
    //double ** do_euler_explicit(int N, int M, float r, double sigma, double K, double S_max, double T)
    double ** result_table = NULL;
    result_table = do_euler_explicit(N, M, r, sigma, K, S_max, T);
    
    
    std::cout<<"The Result!\n";
    for(int i =0; i < M; i++)
    {
        for(int j = 0; j<N-1; j++)
        {
            std::cout<<result_table[i][j]<<" ";
            output_file<<result_table[i][j]<<',';
        }
        std::cout<<result_table[i][N-1]<<std::endl;
        output_file<<result_table[i][N-1]<<"\n";
    }
    
    std::cout<<"    ........................."<<std::endl;
    result_table = do_euler_implicitet(N, M, r, sigma, K, S_max, T);

    std::cout<<"The Result!\n";
    for(int i =0; i < M; i++)
    {
        for(int j = 0; j<N-1; j++)
        {
            std::cout<<result_table[i][j]<<" ";
            output_file<<result_table[i][j]<<',';
        }
        std::cout<<result_table[i][N-1]<<std::endl;
        output_file<<result_table[i][N-1]<<"\n";
    }

    result_table = do_crank_nicolson(N, M, r, sigma, K, S_max, T);
    
    std::cout<<"The Result!\n";
    for(int i =0; i < M; i++)
    {
        for(int j = 0; j<N-1; j++)
        {
            std::cout<<result_table[i][j]<<" ";
            output_file<<result_table[i][j]<<',';
        }
        std::cout<<result_table[i][N-1]<<std::endl;
        output_file<<result_table[i][N-1]<<"\n";
    }

    
}