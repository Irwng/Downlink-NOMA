/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/11/27
Description: include the functions of ML & MPA algorithms in recevier
***********************************************************/
#include "header.h"

/*************************
 *choose the ML algorithm
 *************************/

void ML(SymAfterPPMatrix symAfterPP[Nj][U], SymAfterMPAMatrix* symAfterMPA, 
        SymAfterMapMatrix* masterConstell){

    int j = 0;    
    double NormalBase[U][NJ] = {0.0};            /* the base to normalize */   
	double d[Mpoint] = {0.0};                    /* Euclidean dist of each star point*/   
    double f[Mpoint] = {0.0};                    /* prob of each star point */
    double Q_post[J][Mod] = {0.0};               /* post prob of each bit */
    double MinDist = 0.0;                        /* the minimum Euclidean dist */
    int MinDistPoint = 0;                        /* the minimum Euclidean dist point */

    /* do ML detection per Space-Time block */
    for(int u = 0; u < U; u++){
        for(int nj = 0; nj < Nj; nj++){

            /* initialize */    
            for(j = 0; j < J; j++)
                for(int mod = 0; mod < Mod; ++mod)
                    Q_post[j][mod] = 0;

            MinDist = INT16_MAX;
            MinDistPoint = 0;

            /* calculate the Euclidean distances */
            for(int i = 0; i < Mpoint; i++){
                d[i] = 0;
                for(int nt = 0; nt < Nt; ++nt){
                    for(int m = 0; m < M; ++m){
                        d[i] += pow(abs(symAfterPP[nj][u](nt, m) - masterConstell[i](nt, m)), 2);
                    }
                }                
                if(d[i] < MinDist){
                    MinDist = d[i];
                    MinDistPoint = i;
                }
            }

        #ifdef NoCode 
            /* completely correct */
            bitset<J> bits(MinDistPoint);

            for(j = 0; j < J; j++){
                symAfterMPA[0](u, nj * J + j) = 1 - bits[j]; 
                symAfterMPA[1](u, nj * J + j) = bits[j]; 
            }
            
        #else
            /* calculate the prob of each star point and marginal prob of each bit */
            for(int i = 0; i < Mpoint; i++){
                bitset<6> bit(i);
                f[i] = exp(-0.5*d[i]);//   / N_Var
                for(j = 0; j < J ; ++j){
                    Q_post[j][bit[j]] += f[i]; 
                }
            }
            
            /* normalization & pass the normalized prob of each bit to Viterbi Decoder */
            for(j = 0; j < J; j++){
                NormalBase[u][nj * J + j] = Q_post[j][0] + Q_post[j][1];
                symAfterMPA[0](u, nj * J + j) = Q_post[j][0]/NormalBase[u][nj * J + j];
                symAfterMPA[1](u, nj * J + j) = Q_post[j][1]/NormalBase[u][nj * J + j];
            }
        #endif

        #ifdef DebugMode
            if(nj == 0 && u ==0){
                cout<<"the minDistPoint of user0 is:"<<MinDistPoint<<endl;
                cout<<masterConstell[MinDistPoint]<<endl;
                bitset<J> bit(MinDistPoint);
                for(j = 0; j < J; ++j) cout<<bit[j]<<" ";
                cout<<endl;
                cout<<"distances"<<endl;
                for(int i = 0; i < Mpoint; i++){
                    bitset<J> bit(i);
                    for(j = 0; j < J; ++j) cout<<bit[j];
                    cout<<": "<<d[i]<<"; ";
                    if((i+1)%4==0)cout<<endl;
                }
                
                cout<<"symAfterMPA"<<endl;
                for(int mod = 0; mod < Mod ; ++mod){
                    cout<<"post of "<<mod<<": ";
                    for(j = 0; j < J ; ++j){
                        cout<<symAfterMPA[mod](0, j)<<" "; 
                    }
                    cout<<endl;
                }
            }
        #endif 
        }
    }
}

/*********************************************
 *choose the ML algorithm with Log-likelihood
 *********************************************/

void MLLL(SymAfterPPMatrix symAfterPP[Nj][U], SymAfterMPAMatrix* symAfterMPA, 
          SymAfterMapMatrix* masterConstell){

    int j = 0;    
	double d[Mpoint] = {0.0};                    /* Euclidean dist of each star point*/   
    double MinDist = 0.0;                        /* the minimum Euclidean dist */
    int MinDistPoint = 0;                        /* the minimum Euclidean dist point */

    /* do ML detection per Space-Time block */
    for(int u = 0; u < U; u++){
        for(int nj = 0; nj < Nj; nj++){

            vector<double> tmpd[J][Mod];
            MinDist = INT16_MAX;
            MinDistPoint = 0;

            /* calculate the Euclidean distances */
            for(int i = 0; i < Mpoint; i++){
                d[i] = 0;
                for(int nt = 0; nt < Nt; ++nt){
                    for(int m = 0; m < M; ++m){
                        d[i] += pow(abs(symAfterPP[nj][u](nt, m) - masterConstell[i](nt, m)), 2);
                    }
                }                
                if(d[i] < MinDist){
                    MinDist = d[i];
                    MinDistPoint = i;
                }
            }
        #ifdef NoCode 
            /* completely correct */
            bitset<J> bits(MinDistPoint);

            for(j = 0; j < J; j++){
                symAfterMPA[0](u, nj * J + j) = 1 - bits[j]; 
                symAfterMPA[1](u, nj * J + j) = bits[j]; 
            }
            
        #else
            /* calculate the prob of each star point and marginal prob of each bit */
            for(int i = 0; i < Mpoint; i++){
                bitset<6> bit(i);
                for(j = 0; j < J ; ++j){
                    tmpd[j][bit[j]].push_back(-d[i]/N_Var);// 
                }
            }

            for(j = 0; j < J ; ++j){
                for(int mod = 0; mod < Mod ; ++mod){
                    symAfterMPA[mod](u, nj * J + j) = fmax(tmpd[j][mod]); 
                }
            }
        #endif

        #ifdef DebugMode
            if(nj == 0 && u ==0){
                cout<<"the minDistPoint of user0 is:"<<MinDistPoint<<endl;
                cout<<masterConstell[MinDistPoint]<<endl;
                bitset<J> bit(MinDistPoint);
                for(j = 0; j < J; ++j) cout<<bit[j]<<" ";
                cout<<endl;
                cout<<"distances"<<endl;
                for(int i = 0; i < Mpoint; i++){
                    bitset<J> bit(i);
                    for(j = 0; j < J; ++j) cout<<bit[j];
                    cout<<": "<<d[i]<<"; ";
                    if((i+1)%4==0)cout<<endl;
                }
                cout<<"distances associated with the c[5]= 1"<<endl;
                cout<<"the size of the vector:"<<tmpd[5][1].size()<<endl;
                for(auto tmpdist: tmpd[5][1]){
                    cout<<tmpdist<<" ";
                }
                cout<<endl;
                
                cout<<"symAfterMPA"<<endl;
                for(int mod = 0; mod < Mod ; ++mod){
                    cout<<"post of "<<mod<<": ";
                    for(j = 0; j < J ; ++j){
                        cout<<symAfterMPA[mod](0, j)<<" "; 
                    }
                    cout<<endl;
                }
            }
        #endif 
        }
    }
}
