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
    double MinDist = 0.0;                        /* the minimum Euclidean dist */
    int MinDistPoint = 0.0;                      /* the minimum Euclidean dist point */
    double d[Mpoint] = {0.0};                    /* Euclidean dist of each star point*/   

    for(int nj = 0; nj < Nj; nj++){
        for(int u = 0; u < U; u++){

            MinDist = INT16_MAX;
            MinDistPoint = 0.0;                      /* the minimum Euclidean dist point */

            /* calculate the Euclidean distances and fix them */
            for(int i = 0; i < Mpoint; i++){
                d[i] = 0;
                for(int nr = 0; nr < Nr; ++nr){
                    for(int m = 0; m < M; ++m){
                        d[i] += pow(abs(symAfterPP[nj][u](nr, m) - masterConstell[i](nr, m)), 2);
                    }
                }
                if(d[i] < MinDist){
                    MinDist = d[i];
                    MinDistPoint = i;
                }
            }
            
            bitset<J> bits(MinDistPoint);

            for(j = 0; j < J; j++){
                symAfterMPA[1](u, nj * J + j) = bits[j];
                symAfterMPA[0](u, nj * J + j) = 1 - static_cast<double>(bits[j]); 
            }
        }
    }
}

/**************************
 *choose the MPA algorithm
***************************/

void MPA(SymAfterPPMatrix symAfterPP[Nj][U], SymAfterMPAMatrix* symAfterMPA, SubConstellMatrix& subConstell, char *argv[]){
    
    int i;                                                  /* default indicator */
    int j;                                                  /* indicator for J */
    int k;                                                  /* indicator for resources */
    int m;                                                  /* indicator for modulation */
    int p,e1,e2,flag,found,counter;					        /* frequency indicator */
    int Imax = 10;									        /* max iteration */
	int t,g; 
     /* function nodes connected with variable */
    int Vc[J][Du] = {{0, 1},
                     {0, 2},
                     {0, 3},
                     {1, 2},
                     {1, 3},
                     {2, 3}};                           
    /* variable nodes connected with function */
    int Fc[Nt*M][Dr] = {{0, 1, 2},
                        {0, 3, 4},
                        {1, 3, 5},
                        {2, 4, 5}};  

    switch (*argv[2]){
        case '1':
        case '3':
            break;

        case '2':
        default:
            cout<<"irregular codebook structure!"<<endl;
            return;
    }

    int VNode1,VNode2,FNode1,FNode2,FNode,VNode;
    int fnode[Mod],vnode[Mod];
	double temp;
    double MinDist = 0.0;                                   /* the minimum Euclidean dist */
    double d[Spoint] = {0.0};                               /* Euclidean distance of each point */
	double TH = 0.0001;                                     /* threshold value */
	double sum;
	double unorm[Mod],unorm0[Spoint];                       /* u-normalization probability */
	double ap[J][Mod] = {0.0};                              /* prior probability */
	double Phi[Nt*M][Spoint] = {0.0};                       /* (int)pow(M,Dr)=8 */
	double I_g2v[Nt*M][Dr][Mod];                            /* message passing from FN to its neighboring VNs */
	double I_v2g[J][Du][Mod];					            /* message passing from VN to its neighboring FNs */
    double Q_prev[J][Mod] = {0};                            /* final probability before iterative */
    double Q_post[J][Mod] = {0};                            /* final probability after iterative */


    for(int nj = 0; nj < Nj; nj++){
        for(int u = 0; u < U; u++){
            /* Initialize */
            for(j = 0; j < J; j++){
                for(m = 0; m < Mod; m++){
                    Q_post[j][m] = 0;
                    Q_prev[j][m] = 0;
                }
            }
           
            /* prior probability */
            for(j = 0; j < J; j++){
                for(m = 0; m < Mod; m++){
                    ap[j][m] = (double)1/Mod;
                }
            }

            /* Phi_n */
            for(int nr = 0; nr < Nr; ++nr){
                for(int m = 0; m < M; ++m){
                    sum = 0;
                    MinDist = INT16_MAX;
                    for(i = 0; i < Spoint; i++){
                        d[i] = pow(abs(symAfterPP[nj][u](nr, m) - subConstell(nr*Nt + m, i)), 2);
                        if(d[i] < MinDist){
                            MinDist = d[i];
                        }
                    }

                    for(i = 0; i < Spoint; i++){
                        unorm0[i] = exp(-d[i]/N_Var);
                        sum += unorm0[i];
                    }

                    for(i = 0; i < Spoint; i++){
                        Phi[nr*Nt + m][i] = unorm0[i]/sum;
                    } 
                }
            }
           
            #ifdef DebugMode
                if(nj==0 && u==0){
                    cout<<".................Phi................."<<u<<endl;
                    for(auto &tmp: Phi){
                        for(auto &tmp1: tmp){
                            cout<<tmp1<<" ";
                        }
                        cout<<endl;
                    }
                    cout<<endl;
                }
            #endif

            /* I_v2g[J][Du][M] */
            for(j = 0; j < J; j++){
                for(i = 0; i < Du; i++){
                    for(m = 0; m < Mod; m++){
                        I_v2g[j][i][m] = (double)1/Mod;
                    }
                }
            }
           
            /* iterative */
            for(int I = 0; I < Imax; I++){
                counter = 0;
                /* FN update */
                for(k = 0; k < Nt*M; k++){
                    for(i = 0; i < Dr; i++){
                        
                        if(i == 0){
                            VNode1 = Fc[k][1];					/* first node to pass external message */
                            VNode2 = Fc[k][2];					/* second node to pass external message */
                            for(p = 0; p < Du; p++){
                                if(Vc[VNode1][p] == k){
                                    FNode1 = p;                 /* calculate the target function node which get the message*/
                                }
                                if(Vc[VNode2][p] == k){
                                    FNode2 = p;
                                }
                            }
                        }
                        
                        else if(i == 1){
                            VNode1 = Fc[k][0];
                            VNode2 = Fc[k][2];
                            for(p = 0; p < Du; p++){
                                if(Vc[VNode1][p] == k){
                                    FNode1 = p;
                                }
                                if(Vc[VNode2][p] == k){
                                    FNode2 = p;
                                }
                            }
                        }

                        else{
                            VNode1 = Fc[k][0];
                            VNode2 = Fc[k][1];
                            for(p = 0; p < Du; p++){
                                if(Vc[VNode1][p] == k){
                                    FNode1 = p;
                                }
                                if(Vc[VNode2][p] == k){
                                    FNode2 = p;
                                }
                            }
                        }

                        for(m = 0; m < Mod; m++){                         /* target user*/
                            temp = 0;
                            for(e1 = 0; e1 < Mod; e1++){                  /* accumulation of first user */
                                for(e2 = 0; e2 < Mod; e2++){              /* accumulation of second user */
                                    /* symbol in SC points decided by m,e1,e2*/
                                    if(i == 0){
                                        t =	m*(int)pow(Mod,Dr-1) + e1*(int)pow(Mod,Dr-2) + e2*(int)pow(Mod,Dr-3);
                                    }
                                    else if(i == 1){
                                        t =	e1*(int)pow(Mod,Dr-1) + m*(int)pow(Mod,Dr-2) + e2*(int)pow(Mod,Dr-3);
                                    }
                                    else{
                                        t =	e1*(int)pow(Mod,Dr-1) + e2*(int)pow(Mod,Dr-2) + m*(int)pow(Mod,Dr-3);
                                    }
                                    temp += Phi[k][t] * I_v2g[VNode1][FNode1][e1] * I_v2g[VNode2][FNode2][e2];
                                }//UE2
                            }//UE1
                            I_g2v[k][i][m] = temp;
                        }//Self's 3 symbols
                    }//Dr
                }//K

                /* VN update */
                for(j = 0; j < J; j++){
                    for(i = 0; i < Du; i++){
                        sum = 0;
                        g = (i == 0)?1:0;
                        FNode = Vc[j][g];
                        for(k = 0; k < Dr; k++){
                            if(Fc[FNode][k] == j){
                                VNode = k;
                            }
                        }

                        for(m = 0; m < Mod; m++){
                            unorm[m] = ap[j][m]*I_g2v[FNode][VNode][m];         /* calculate the normalization factor*/
                            sum += unorm[m];
                        }
                        for(m = 0; m < M; m++){
                            I_v2g[j][i][m] = unorm[m]/sum;
                        }
                    }//Du
                }//VN update J

                /* probability accounts */
                for(j = 0; j < J; j++){
                    found = 0;
                    flag = 0;
                    for(k = found; k < Nt*M; k++){
                        for(i = 0; i < Dr; i++){
                            if(Fc[k][i] == j){
                                fnode[flag] = k;
                                vnode[flag] = i;
                                found = found + 1;
                                flag ++;
                                break;
                            }
                        }
                        if(flag == Du)
                            break;
                    }

                    for(m = 0; m < Mod; m++){
                        Q_post[j][m] = ap[j][m] * I_g2v[ fnode[0] ][ vnode[0] ][m] 
                                                * I_g2v[ fnode[1] ][ vnode[1] ][m];
                        if((Q_post[j][m] - Q_prev[j][m]) > TH){
                            Q_prev[j][m] = Q_post[j][m];
                        }
                        else
                            counter++;
                    }
                }

                /*  threshold judgment */
                if(counter == Mod*J)
                    break;
            }//Iterative

            /* normalization */
            double temp = 0.0;
            for(j = 0; j < J; j++){
                temp = Q_post[j][0] + Q_post[j][1];
                Q_post[j][0] = Q_post[j][0]/temp;
                Q_post[j][1] = Q_post[j][1]/temp;
            }

            #ifdef DebugMode
                if(nj==0 && u==0){
                    cout<<".................MPA:Q_post................."<<endl;
                    for(j = 0; j < J; j++){
                        cout<<"("<<Q_post[j][0]<<","<<Q_post[j][1]<<")"<<" "; 
                    }
                    cout<<endl;
                }
            #endif

            for(j = 0; j < J; j++){
                symAfterMPA[0](u, nj * J + j) = Q_post[j][0];
                symAfterMPA[1](u, nj * J + j) = Q_post[j][1]; 
            }
        }
    }
    #ifdef DebugMode
        cout<<".................MPA:SymAfterMPA................."<<endl;
        cout<<symAfterMPA<<endl;
    #endif
}
