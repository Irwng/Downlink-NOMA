/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/11/27
Description: include the functions of ML & MPA algorithms in recevier
***********************************************************/
#include "header.h"

/*************************
 *choose the ML algorithm
 *************************/

void ML(SymAfterPPMatrix* symAfterPP, SymAfterMPAMatrix* symAfterMPA, 
        SymAfterMapMatrix* masterConstell){

    int j = 0;    
    double MinDist = 0.0;                        /* the minimum Euclidean dist */
    int MinDistPoint = 0.0;                      /* the minimum Euclidean dist point */
    double d[Mpoint] = {0.0};                    /* Euclidean dist of each star point*/   

    for(int u = 0; u < U; u++){

        MinDist = INT16_MAX;
        MinDistPoint = 0.0;                      /* the minimum Euclidean dist point */

        /* calculate the Euclidean distances and fix them */
        for(int i = 0; i < Mpoint; i++){
            d[i] = 0;
            for(int nr = 0; nr < Nr; ++nr){
                for(int m = 0; m < M; ++m){
                    d[i] += pow(abs(symAfterPP[u](nr, m) - masterConstell[i](nr, m)), 2);
                }
            }
            if(d[i] < MinDist){
                MinDist = d[i];
                MinDistPoint = i;
            }
        }
    #ifdef DebugMode
        if(u==1){
            cout<<"MinDistpoint: "<<MinDistPoint<<endl;
            bitset<J> tmpbit(MinDistPoint);
            for(j = 0; j < J; j++){
                cout<<tmpbit[j]<<" ";
            }
            cout<<endl;
            cout<<"MinDistpoint: "<<masterConstell[MinDistPoint]<<endl;
        } 
    #endif
        bitset<J> bits(MinDistPoint);

        for(j = 0; j < J; j++){
            symAfterMPA[1](u, j) = bits[j];
            symAfterMPA[0](u, j) = 1 - static_cast<double>(bits[j]); 
        }
    }
}
