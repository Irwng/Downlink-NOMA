/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/11/27
Description: include the functions of ML & MPA algorithms in recevier
***********************************************************/
#include "header.h"

/*************************
 *choose the ML algorithm
 *************************/

void ML(SymAfterFCMatrix symAfterFC[Nj][U], 
        SymAfterMPAMatrix* symAfterMPA, 
        SymAfterMapMatrix* masterConstell, 
        CSIMatrix* h, 
        MapMatrix& map){

    int j = 0;    
    double NormalBase = 0.0;                     /* the base to normalize */   
	double d[Mpoint] = {0.0};                    /* Euclidean dist of each star point*/   
    double f[Mpoint] = {0.0};                    /* prob of each star point */
    double Q_post[J][Mod] = {0.0};               /* post prob of each bit */
    SymAfterFCMatrix MasterConstellFixed[Mpoint];

    for(int u = 0; u < U; u++){
        /* adjust the master-constellation of each user */
        for(int i = 0; i < Mpoint; i++){                                       /*ideal master constellation coding*/\
            MasterConstellFixed[i] = h[u] * masterConstell[i];
        }
        #ifdef DebugMode
            if(u == 0){        
                cout<<"...................MasterConstellFixed..................."<<endl;
                for(int i = 0; i < Mpoint; ++i){
                    cout<<"...................Mpoint:"<< i<<" ..................."<<endl;
                    cout<<MasterConstellFixed[i]<<endl;
                }
            }        
        #endif 

        for(int nj = 0; nj < Nj; nj++){

            /* initialize */    
            for(j = 0; j < J; j++)
                for(int m = 0; m < Mod; ++m)
                    Q_post[j][m] = 0;
            
            for(int i = 0; i < Mpoint; i++){
                // bitset<J> bit(i);
                d[i] = 0;
                for(int nr = 0; nr < Nr; ++nr){
                    for(int m = 0; m < M; ++m){
                        d[i] += pow(abs(symAfterFC[nj][u](nr, m) - MasterConstellFixed[i](nr, m)), 2);
                    }
                }
            }
            #ifdef DebugMode
                if(nj == 0){
                    cout<<"...................the minDist:"<<MinDistPoint<<" of user: "<<u<<"..................."<<endl;
                    cout<<MasterConstellFixed[MinDistPoint]<<endl;
                }
            #endif 
            
            /* calculate the prob of each points and Marginal prob of each user */
            for(int i = 0; i < Mpoint; i++){
                bitset<J> bit(i);
                // f[i] = exp(- d[i]); 
                f[i] = exp(- d[i] / N_Var);  /* the format according to theory */ 
                for(j = 0; j < J ; ++j){
                    Q_post[j][ bit[j] ] += f[i]; 
                }
            }

            /* normalization */
            for(j = 0; j < J; j++){
                NormalBase = Q_post[j][0] + Q_post[j][1];
                Q_post[j][0] = Q_post[j][0]/NormalBase;
                Q_post[j][1] = Q_post[j][1]/NormalBase;
            }

            /* pass data */
            for(j = 0; j < J; j++){
                symAfterMPA[0](u, nj * J + j) = Q_post[j][0]; 
                symAfterMPA[1](u, nj * J + j) = Q_post[j][1]; 
            }
        }
    }

    #ifdef DebugMode
        cout<<endl<<".................ML:SymAfterMPA................."<<endl;
        cout<<symAfterMPA<<endl;
    #endif
}
