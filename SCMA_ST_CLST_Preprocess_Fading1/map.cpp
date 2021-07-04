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

    for(int u = 0; u < U; u++){
        for(int nj = 0; nj < Nj; nj++){

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
        #ifdef DebugMode
            if(u==1&&nj==0){
                cout<<"MinDistpoint: "<<MinDistPoint<<endl;
                cout<<"MinDistpoint: "<<masterConstell[MinDistPoint]<<endl;
            } 
        #endif
            bitset<J> bits(MinDistPoint);

            for(j = 0; j < J; j++){
                symAfterMPA[1](u, nj * J + j) = bits[j];
                symAfterMPA[0](u, nj * J + j) = 1 - static_cast<double>(bits[j]); 
            }
        }
    }
}


// void Receiver_OSIC(SymAfterFCMatrix symAfterFC[Nj][U], SymAfterMapMatrix* symAfterMap, 
//               CSIMatrix* h, SymAfterPPMatrix symAfterPP[Nj][U], 
//               PPMatrix* v, char* argv[]){

//     /* AWGN */
//     SymAfterFCMatrix tmp;
//     for(int nr = 0; nr < Nr; ++nr){
//         tmp(nr) = AWGN(N_Var);
//     }
//     symAfterFC = h * symAfterMap + tmp;//  
    
//     bool index_array[Nt]; // index of the decoded signal
//     for(auto& a: index_array) a = false;
//     for(int nt = 0; nt < Nt; ++nt){
//         /* update the channel gain matrix */
//         ComplexD* h_buff = new ComplexD[Nt*(Nt-nt)];
        
//         int column = 0;
//         for(int col = 0; col < Nt; ++col){
//             if(index_array[col] == true) continue;
//             for(int row = 0; row < Nt; ++row){
//                 h_buff[column*Nt + row] = h(row, col);
//             }        
//             column++;
//         }
        
//         Matrix<ComplexD, Dynamic, Dynamic> h_temp;
//         h_temp = Map<Matrix<ComplexD, Dynamic, Dynamic>>(h_buff, Nt, Nt-nt);                
        
//     #ifdef DebugMode
//         cout<<"nt: "<<nt<<endl;
//         cout<<"h_buff: "<<endl;
//         for(int k = 0; k < Nt*(Nt-nt); ++k){
//             cout<<h_buff[k]<<" ";
//         }
//         cout<<endl;
//         cout<<"h_temp:"<<endl<<h_temp<<endl;
//     #endif
//         delete[] h_buff;

//         /* calculate the h_temp's conjugate transpose matrix */
//         Matrix<ComplexD, Dynamic, Dynamic> ConjTrans_H;
//         ConjTrans_H = h_temp.conjugate().transpose();

//         /* calculate the base matrix */
//         Matrix<ComplexD, Dynamic, Dynamic> Denomiator_H;
//         Denomiator_H = ConjTrans_H * h_temp;
//         for(int i = 0; i < Nt-nt; ++i) Denomiator_H(i, i) += ComplexD(N_Var, 0);
        
//         /* calculate the preprocess matrix */
//         Matrix<ComplexD, Dynamic, Dynamic> V_temp;
//         V_temp = Denomiator_H.inverse() * ConjTrans_H;

//         Matrix<ComplexD, Dynamic, Dynamic> VzfH;
//         VzfH = V_temp * h_temp;

//         /* calculate the parameters need by SINR */
//         double powersum[Nt-nt];
//         double noisesum[Nt-nt];
//         for(int i = 0; i < Nt-nt; ++i){
//             powersum[i] = 0;
//             noisesum[i] = 0;
//             for(int j = 0; j < Nt-nt; ++j)
//                 powersum[i] += pow(abs(VzfH(i, j)), 2);

//             for(int j = 0; j < Nt; ++j)
//                 noisesum[i] += pow(abs(V_temp(i, j)), 2);
//         }
        
//         float SINR[Nt-nt] = {0};
//         float denominator = 0.0;
        
//         /* calculate the maximum SINR */        

//         for(int nr = 0; nr < Nr-nt; ++nr){
//             denominator = powersum[nr] - pow(abs(VzfH(nr,nr)), 2) + N_Var * noisesum[nr];
//             SINR[nr] = pow(abs(VzfH(nr,nr)), 2) / denominator;
//         }
//     #ifdef DebugMode
//         cout<<"SINR"<<endl;
//         for(auto b: SINR) cout<<b<<" ";
//         cout<<endl;
//     #endif

//         /* find the max SINR subchannel */
//         int index_max = 0;
//         for(int i = 1; i < Nt-nt; ++i){
//             if(SINR[i] > SINR[index_max]) index_max = i;         
//         }

//         /* find the order in ture line*/
//         /* initialize */
//         for(int i = 0; i < Nt; ++i){
//             if(index_array[i] == false){
//                 column = i;
//                 break;    
//             }
//         }

//         for(int i = 0; i < index_max; ++i){
//             column++;
//             while(column < Nt&& index_array[column] == true){
//                 column++;
//             }
//         }

//         if(index_array[column] == true) cout<<"error!";
//         index_array[column] = true;

//     #ifdef DebugMode
//         cout<<"index_max: "<<index_max<<endl;
//         cout<<"column: "<<column<<endl;
//         cout<<"index_array:"<<endl;
//         for(auto b: index_array) cout<<b<<" ";
//         cout<<endl;
//     #endif

//         /* decode the signal on the max SINR subchannel in this loop */
//         ComplexD x_temp(0,0);
//         for(int nr = 0; nr < Nt; ++nr) x_temp += V_temp(index_max, nr) * symAfterFC(nr);
//         /* check the performance of the post processing */
//     #ifdef BPSK
//         decode(column) = static_cast<int>(x_temp.real() > 0);
//         x_temp = ComplexD(static_cast<double>(2*decode(column) - 1),0);
//     #endif
        
//         /* cancell the interference */
//         for(int nr = 0; nr < Nt; ++nr) symAfterFC(nr) -= h_temp(nr, index_max) *sqrt(1.0/Nt) *x_temp;
//     }
 
// }