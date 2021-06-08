/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: include the basic functions of transmitter & recevier & channels
***********************************************************/
#include "header.h"

double N_Var;                                       /* variance of Noise */
double BER_TOTAL_Nocode = 0;                               /* total number of error symbols */
double BER_Nocode = 0;
double BER_TOTAL_Code = 0;                               /* total number of error symbols */
double BER_Code = 0;

SubConstellMatrix SubConstell;                      /* sub-constelltion point */
SymAfterMapMatrix MasterConstell[Mpoint];           /* Master-constelltion point */
SourceMatrix Source;                                /* source codewords */
CodeMatrix Code;                                    /* codewords after coding */
ModuMatrix Modu;                                    /* symbols after modulation */
MapMatrix F;                                        /* Mapping matrix */
SymAfterMapMatrix SymAfterMap[Nj];                  /* signals after mapping, Nj*Nt*M */
CSIMatrix H[U][Nj];                                     /* channel parameters , U*Nt*Nr */
CSIMatrix WeightedIdentityMatrix;                   /* MMSE assistance matrix */
PPMatrix V[U];                                      /* postprocessing matrix, Nj*Nr*M */
SymAfterBFMatrix SymAfterBF[Nj];                    /* receiving signals , Nj*Nr*M */
SymAfterFCMatrix SymAfterFC[Nj][U];
SymAfterPPMatrix SymAfterPP[Nj][U];                 /* receiving signals after postprocessing, Nj*Nr*M */
SymAfterMPAMatrix SymAfterMPA[Mod];
CodeMatrix CodeEsti[U];
SourceMatrix SourceEsti[U]; 
PPMatrix Q;

ComplexD p1;
ComplexD p2;
ComplexD p3;


void Initialize(char* argv[]){

    cout<<"SCMA_ST_Fadding_Detection_Latin"<<endl;
    outfile.open("SCMA_ST_Fadding_Detection_Latin.txt", ios::out|ios::app);
    outfile<<"....................................................."<<endl;
    outfile<<"SCMA_ST_Fadding_Detection_Latin"<<endl;   

    #ifdef OSIC
        switch (*argv[4]){
            case '1':
                cout<<"OSIC-SINR"<<endl;
                outfile<<"OSIC-SINR"<<endl;
                break;
            
            case '2':
                cout<<"OSIC-SNR"<<endl;
                outfile<<"OSIC-SNR"<<endl;
                break;

            default:
                break;
        }
    #endif


    switch (*argv[3]){
        case '1':
            cout<<"ZF"<<endl;
            outfile<<"ZF"<<endl;
            break;
        
        case '2':
            cout<<"MMSE"<<endl;
            outfile<<"MMSE"<<endl;
            break;

        default:
            break;
    }

    #ifdef LL
        cout<<"ML-LL"<<endl;
        outfile<<"ML-LL"<<endl;
    #else
        cout<<"ML"<<endl;
        outfile<<"ML"<<endl;
    #endif

#ifdef NoCode
    cout<<"NoCode, block length="<<NJ<<endl;
    outfile<<"NoCode, block length="<<NJ<<endl;
#else
    cout<<"Convolution Code & soft Viterbi decoder: rate=1/"<<ReciRate<<", block length="<<NJ<<endl;
    outfile<<"Convolution Code & soft Viterbi decoder: rate=1/"<<ReciRate<<", block length="<<NJ<<endl;
#endif
}


void InitMapMatrix(MapMatrix& map, SubConstellMatrix& subConstell, 
                   SymAfterMapMatrix* masterConstell, 
                   char* argv[]){
    
    double weight;
    double RF;

    /* initialize the mapping matrix's codewords */
    switch(*argv[1]){
        /* UDM */
        case '1':
            /* 202,2,4,9 */
            cout<<"UDM-202,2,4,9"<<endl;
            outfile<<"UDM-202,2,4,9"<<endl;

            weight = sqrt(202);
            p1 = ComplexD(2/weight, 9/weight);
            p2 = ComplexD(4/weight, 4/weight);
            p3 = ComplexD(9/weight, 2/weight);
            break;

        case '2':
            /* 42, 1, 2, 4 */
            cout<<"UDM-42,1,2,4"<<endl;
            outfile<<"UDM-42,1,2,4"<<endl;

            weight = sqrt(42);
            p1 = ComplexD(1/weight, 4/weight);
            p2 = ComplexD(2/weight, 2/weight);
            p3 = ComplexD(4/weight, 1/weight);
            break;

        case '3':
            /* 178, 3, 4, 8 */
            cout<<"UDM-178,3,4,8"<<endl;
            outfile<<"UDM-178,3,4,8"<<endl;

            weight = sqrt(178);
            p1 = ComplexD(3/weight, 8/weight);
            p2 = ComplexD(4/weight, 4/weight);
            p3 = ComplexD(8/weight, 3/weight);
            break;

        case '4':
            /* DT-1,1+1i,2i */
            cout<<"UDM-1,1+i,2i"<<endl;
            outfile<<"UDM-1,1+i,2i"<<endl;

            weight = sqrt(7);
            p1 = ComplexD(1/weight, 0/weight);
            p2 = ComplexD(1/weight, 1/weight);
            p3 = ComplexD(0/weight, 2/weight);
            break;


        /* PhaseRotation */
        case '6':
            /* RF=pi/4 */
            cout<<"PhaseRotation-pi/4"<<endl;
            outfile<<"PhaseRotation-pi/4"<<endl;

            weight = sqrt(3);
            RF = PI/4;
            p1 = ComplexD(cos(0*RF)/weight, sin(0*RF)/weight);
            p2 = ComplexD(cos(1*RF)/weight, sin(1*RF)/weight);
            p3 = ComplexD(cos(2*RF)/weight, sin(2*RF)/weight);    
            break;

        case '7':
            /* RF=pi/5 */
            cout<<"PhaseRotation-pi/5"<<endl;
            outfile<<"PhaseRotation-pi/5"<<endl;

            weight = sqrt(3);
            RF = PI/5;
            p1 = ComplexD(cos(0*RF)/weight, sin(0*RF)/weight);
            p2 = ComplexD(cos(1*RF)/weight, sin(1*RF)/weight);
            p3 = ComplexD(cos(2*RF)/weight, sin(2*RF)/weight);    
            break;

        case '8':
            /* RF=pi/6 */
            cout<<"PhaseRotation-pi/6"<<endl;
            outfile<<"PhaseRotation-pi/6"<<endl;

            weight = sqrt(3);
            RF = PI/6;
            p1 = ComplexD(cos(0*RF)/weight, sin(0*RF)/weight);
            p2 = ComplexD(cos(1*RF)/weight, sin(1*RF)/weight);
            p3 = ComplexD(cos(2*RF)/weight, sin(2*RF)/weight);    
            break;

        default:
            break;
    }

    /* initialize the mapping matrix's structure */    
    switch(*argv[2]){
        
        /* Latin */
        case '1':
            cout<<"Latin"<<endl;
            outfile<<"Latin"<<endl;

              /* Nt= 0 Nt= 1 */
            map<<p1,p2,p0,p0,
                 p2,p0,p3,p0,
                 p3,p0,p0,p1,
                 p0,p3,p1,p0,
                 p0,p1,p0,p2,
                 p0,p0,p2,p3;
            break;

        /* CLST */
        case '2':
            cout<<"CLST"<<endl;
            outfile<<"CLST"<<endl;

              /* Nt= 0 Nt= 1 */
            map<<p1,p0,p0,p0,
                 p2,p1,p0,p0,
                 p3,p2,p1,p0,
                 p0,p3,p2,p1,
                 p0,p0,p3,p2,
                 p0,p0,p0,p3;
            break;
        
        /* Orthogonal */
        case '3':
            cout<<"Orthogonal"<<endl;
            outfile<<"Orthogonal"<<endl;

            /* Nt= 0 Nt= 1 */
            map<<p1,p0,p0,p1,
                 p2,p0,p0,p2,
                 p3,p0,p0,p3,
                 p0,p1,p1,p0,
                 p0,p2,p2,p0,
                 p0,p3,p3,p0;
            break;
        
        /* half-Latin */
        case '4':
            cout<<"half-Latin"<<endl;
            outfile<<"half-Latin"<<endl;

              /* Nt= 0 Nt= 1 */
            map<<p1,p1,p0,p0,
                 p2,p0,p1,p0,
                 p3,p0,p0,p1,
                 p0,p2,p2,p0,
                 p0,p3,p0,p2,
                 p0,p0,p3,p3;
            break;
        default:
            break;
    }

    cout<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"M: "<<M<<setw(10)<<"J: "<<J<<setw(10)<<"NLoop: "<<NLoop<<endl;   
    outfile<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"M: "<<M<<setw(10)<<"J: "<<J<<setw(10)<<"NLoop: "<<NLoop<<endl;   
    
    cout<<"EbN0dB"<<setw(20)<<"BER_Nocode"<<setw(20)<<"BER_Code"<<endl;
    outfile<<"EbN0dB"<<setw(20)<<"BER_Nocode"<<setw(20)<<"BER_Code"<<endl;

    /* initialize the subconstellation */
    int a,b;
    double subMsg[Spoint][Dr] = {0.0};

    for(int i = 0; i < Spoint; i++){
        for(int j = Dr; j > 0; j--){
            a = (int)pow(2,-j+Dr+1);
            b = (int)pow(2,-j+Dr);
            if(i%a >= b)
                subMsg[i][j-1] = 1;                             /* 1 - -1,0 - 1 */
            else
                subMsg[i][j-1] = -1;
        }
    } 

    int flag = 0;
    for(int k = 0; k < Nt*M; k++){
		for(int i = 0; i < Spoint; i++){
			flag = 0;
			for(int index = 0; index < J; index++){
				if(map(index,k) != p0){
					subConstell(k,i) += map(index, k) * subMsg[i][flag];
					flag++;
				}
			}
		}
	}

    /* initialize the Masterconstellation */
    for(int i = 0; i < Mpoint; ++i){                                     
        bitset<J> bit(i);
        for(int nt = 0; nt < Nt; ++nt){
            for(int m = 0; m < M; ++m){
                masterConstell[i](nt, m) = p0;
                for(int j = 0; j < J; ++j){
                    masterConstell[i](nt, m) += map(j, nt * Nt + m) * static_cast<double>(2 * bit[j] - 1);
                }
            }
        }
    }

    for(int nt = 0; nt < Nt; ++nt){
        for(int nt = 0; nt < Nt; ++nt){
            Q(nt, nt) = p0;
        }
    }

    for(int mpoint = 0; mpoint < Mpoint; ++mpoint)
        Q += (masterConstell[mpoint].conjugate().transpose()) * masterConstell[mpoint];
    Q /= static_cast<double>(Mpoint);
}


/* key point */
void ChannelInitialize(int ebN0dB){

    double ebN0 = pow(10, (double)ebN0dB/10);

#ifdef NoCode
    double snr = (static_cast<double>(J)/(Nt*M)) * ebN0;
#else
    double snr = (static_cast<double>(J)/(Nt*M*ReciRate))*ebN0;//
#endif

    double N0 = power*Nt / snr;
    N_Var = N0/2;
    for(int nt = 0; nt < Nt; ++nt){
        WeightedIdentityMatrix(nt, nt) = ComplexD(N_Var, 0);
    }
    BER_TOTAL_Nocode = 0;
    BER_TOTAL_Code = 0;
    #ifdef DebugMode
        cout<<"N_Var: "<<N_Var<<endl;
    #endif
}


void BitSource(SourceMatrix& source){

    /* random number generator: 0 or 1 */
    for(int i = 0; i < LenBit; i++){
        source(i)= rand()%2;
    }
    #ifdef DebugMode
        cout<<".................source................."<<endl;
        cout<<source<<endl;
    #endif
}


void ConvEncoder(SourceMatrix& source, CodeMatrix& code){

    /* (5,7) in rate = 1/2 */
    int reg1=0, reg2=0;
    for(int i = 0; i < LenBit; i++){
        code(ReciRate*i + 0) = source(i) xor reg2;
        code(ReciRate*i + 1) = source(i) xor reg1 xor reg2;

        reg2 = reg1;
        reg1 = source(i);
    }

    #ifdef DebugMode
        cout<<".................code................."<<endl;
        for(int j = 0; j < J; j++)
            cout<<code(j)<<" ";
        cout<<endl;
    #endif
}


void Modulation(CodeMatrix& code, ModuMatrix& modu){

    /* BPSK: 0->-1, 1->1 */
    for(int nj = 0; nj < Nj; nj++){
        for(int j = 0; j < J; j++){
            modu(nj, j) = 2 * code(nj * J + j) - 1;
        }
    }
    #ifdef DebugMode
        cout<<".................bpsk signal................."<<endl;
        for(int j = 0; j < J; j++)
            cout<<modu(0, j)<<" ";
        cout<<endl;
    #endif
}


void Mapper(ModuMatrix& modu, SymAfterMapMatrix* symAfterMap, MapMatrix& map){

    for(int nj = 0; nj < Nj; nj++){
        for(int nt = 0; nt < Nt; nt++){
            for(int m = 0; m < M; m++){
                symAfterMap[nj](nt, m) = p0;
                for(int j = 0; j < J; j++){
                    symAfterMap[nj](nt, m) +=  map(j, m + nt * Nt) * modu(nj, j);
                }
            }
        }
#ifdef DebugMode
    if(nj == 0){
        cout<<"................symaftermap[0]................ "<<endl;
        cout<<symAfterMap[0]<<endl; 
    }
#endif
    }
}

void FadingChannel(CSIMatrix h[U][Nj]){

    /* channel parameters */
    for(int u = 0; u < U; ++u){
        for(int nj_L = 0; nj_L < Nj_L; ++nj_L){
            for(int nr = 0; nr < Nr; nr++){
                for(int nt = 0; nt < Nt;nt++){
                    double GG = sqrt(-2*log(randN() + 0.000000001));
                    double B = randN();
                    for(int l = 0; l < L; ++l){
                        h[u][nj_L*L + l](nr, nt) = ComplexD(sqrt(0.5) * GG * cos(2*PI*B),
                                                            sqrt(0.5) * GG * sin(2*PI*B));
                    }
                }
            }
        }
    }
    // #ifdef DebugMode
    //     cout<<"...........h of user0........... "<<endl;
    //     for(int nj = 0; nj < Nj; ++nj){
    //         cout<<h[0][nj]; 
    //     }
    //     cout<<endl;
    // #endif
}


ComplexD AWGN(double nvar){

    double GG = sqrt(-2*log(randN() + 0.000000001));
    double B = randN();
    ComplexD GuassN(sqrt(nvar)* GG * cos(2*PI*B), sqrt(nvar)* GG * sin(2*PI*B));
    return GuassN;
}


void Receiver(SymAfterFCMatrix symAfterFC[Nj][U], SymAfterBFMatrix* symAfterBF, 
              CSIMatrix h[U][Nj], SymAfterPPMatrix symAfterPP[Nj][U], 
              PPMatrix* v, char* argv[]){

    CSIMatrix ConjTrans_H;
    CSIMatrix Denomiator_H;

    for(int nj = 0; nj < Nj; nj++){
        for(int u = 0; u < U; u++){
            
            /* AWGN */
            SymAfterFCMatrix tmp;
            for(int nr = 0; nr < Nr; nr++){
                for(int m = 0; m < M; m++){
                    tmp(nr, m) = AWGN(N_Var);
                }
            }  

            symAfterFC[nj][u] = h[u][nj] * symAfterBF[nj] + tmp; //

            /* Prepeocess */
            ConjTrans_H = h[u][nj].conjugate().transpose();
            
            switch (*argv[3]){
                case '1':
                    Denomiator_H = ConjTrans_H * h[u][nj];
                    v[u] = Denomiator_H.inverse() * ConjTrans_H;
                    break;
                
                case '2':
                    Denomiator_H = ConjTrans_H * h[u][nj] + WeightedIdentityMatrix;
                    v[u] = Denomiator_H.inverse() * ConjTrans_H;
                    // Denomiator_H = ConjTrans_H * Q * h[u][nj] + WeightedIdentityMatrix;
                    // v[u] = Denomiator_H.inverse() * ConjTrans_H * Q;
                    break;

                default:
                    break;
            }
            symAfterPP[nj][u] = v[u] * symAfterFC[nj][u];
        }
    }
}


void Receiver_OSIC(SymAfterFCMatrix symAfterFC[Nj][U], SymAfterBFMatrix* symAfterBF,
                   CSIMatrix h[U][Nj], SymAfterPPMatrix symAfterPP[Nj][U],
                   PPMatrix* v, SubConstellMatrix& subConstell, char* argv[]){

    for(int nj = 0; nj < Nj; nj++){
        for(int u = 0; u < U; u++){
            
            /* AWGN */
            SymAfterFCMatrix tmp;
            for(int nr = 0; nr < Nr; ++nr){
                tmp(nr) = AWGN(N_Var);
            }
            for(int nr = 0; nr < Nr; ++nr)
                for(int m = 0; m < M; ++m)
                    symAfterPP[nj][u](nr, m) = p0;
            
            symAfterFC[nj][u] = h[u][nj] * symAfterBF[nj] + tmp; // 

            bool index_array[Nt]; // index of the decoded signal
            for(auto& a: index_array) a = false;
            
            for(int nt = 0; nt < Nt; ++nt){
                /* update the channel gain matrix */
                ComplexD* h_buff = new ComplexD[Nt*(Nt-nt)];
                
                int column = 0;
                for(int col = 0; col < Nt; ++col){
                    if(index_array[col] == true) continue;
                    for(int row = 0; row < Nt; ++row){
                        h_buff[column*Nt + row] = h[u][nj](row, col);
                    }        
                    column++;
                }
                
                Matrix<ComplexD, Dynamic, Dynamic> h_temp;
                h_temp = Map<Matrix<ComplexD, Dynamic, Dynamic>>(h_buff, Nt, Nt-nt);                
                
            #ifdef DebugMode
                cout<<"nt: "<<nt<<endl;
                cout<<"h_buff: "<<endl;
                for(int k = 0; k < Nt*(Nt-nt); ++k){
                    cout<<h_buff[k]<<" ";
                }
                cout<<endl;
                cout<<"h_temp:"<<endl<<h_temp<<endl;
            #endif
                
                delete[] h_buff;

                /* calculate the h_temp's conjugate transpose matrix */
                Matrix<ComplexD, Dynamic, Dynamic> ConjTrans_H;
                ConjTrans_H = h_temp.conjugate().transpose();

                /* calculate the base matrix */
                Matrix<ComplexD, Dynamic, Dynamic> Denomiator_H;
                Denomiator_H = ConjTrans_H * h_temp;
                
                if(*argv[3] == '2')
                    for(int i = 0; i < Nt-nt; ++i) 
                        Denomiator_H(i, i) += ComplexD(N_Var, 0);

                /* calculate the preprocess matrix */
                Matrix<ComplexD, Dynamic, Dynamic> V_temp;
                V_temp = Denomiator_H.inverse() * ConjTrans_H;

                Matrix<ComplexD, Dynamic, Dynamic> VzfH;
                VzfH = V_temp * h_temp;

                /* calculate the parameters need by SINR or SNR*/
                double powersum[Nt-nt];
                double noisesum[Nt-nt];
                for(int i = 0; i < Nt-nt; ++i){
                    powersum[i] = 0;
                    noisesum[i] = 0;
                    for(int j = 0; j < Nt-nt; ++j)
                        powersum[i] += pow(abs(VzfH(i, j)), 2);

                    for(int j = 0; j < Nt; ++j)
                        noisesum[i] += pow(abs(V_temp(i, j)), 2);
                }
                
                float SINR[Nt-nt] = {0};
                float denominator = 0.0;
                
                /* calculate the maximum SINR */        
                switch (*argv[4]){
                    case '1':
                        for(int nr = 0; nr < Nr-nt; ++nr){
                            denominator = powersum[nr] - pow(abs(VzfH(nr,nr)), 2) + N_Var * noisesum[nr];
                            SINR[nr] = pow(abs(VzfH(nr,nr)), 2) / denominator;
                        }
                        break;
                    case '2':
                        for(int nr = 0; nr < Nr-nt; ++nr){
                            denominator = N_Var * noisesum[nr];
                            SINR[nr] = 1 / denominator;
                        }
                        break;                    
                    default:
                        break;
                }
            #ifdef DebugMode
                cout<<"SINR"<<endl;
                for(auto b: SINR) cout<<b<<" ";
                cout<<endl;
            #endif

                /* find the max SINR subchannel */
                int index_max = 0;
                for(int i = 1; i < Nt-nt; ++i){
                    if(SINR[i] > SINR[index_max]) index_max = i;         
                }

                /* find the order in ture line*/
                /* initialize */
                column = 0;
                for(int i = 0; i < Nt; ++i){
                    if(index_array[i] == false){
                        column = i;
                        break;    
                    }
                }

                for(int i = 0; i < index_max; ++i){
                    column++;
                    while(column < Nt&& index_array[column] == true){
                        column++;
                    }
                }

                if(index_array[column] == true) cout<<"error!";
                index_array[column] = true;

            #ifdef DebugMode
                cout<<"index_max: "<<index_max<<endl;
                cout<<"column: "<<column<<endl;
                cout<<"index_array:"<<endl;
                for(auto b: index_array) cout<<b<<" ";
                cout<<endl;
            #endif
                
                /* check the performance of the post processing */
                for(int m = 0; m < M; ++m){
                    
                    int SubMinDistPoint = 0;
                    for(int nr = 0; nr < Nr; ++nr) 
                        symAfterPP[nj][u](column, m) += V_temp(index_max, nr) * symAfterFC[nj][u](nr, m);
                    
                    double SubMinDist = abs(symAfterPP[nj][u](column, m) - subConstell(0, 0));
                    double tmpDist;  
                    /* ML */
                    for(int spoint = 1; spoint < Spoint; ++spoint){
                        tmpDist = abs(symAfterPP[nj][u](column, m) - subConstell(0, spoint));
                        if(tmpDist < SubMinDist ){
                            SubMinDistPoint = spoint;
                            SubMinDist = tmpDist; 
                        }
                    }
                    /* calculate the Interference signals */
                    for(int nr = 0; nr < Nr; ++nr) 
                        symAfterFC[nj][u](nr, m) -= h_temp(nr, index_max) * subConstell(0, SubMinDistPoint);
                }
            }
        }
    }
}


void Compare(CodeMatrix* codeEsti, CodeMatrix& code, 
             SourceMatrix* sourceEsti, SourceMatrix& source){

    for(int u = 0; u < U; u++){
    #ifdef NoCode
        for(int nj = 0; nj < NJ; nj++){
            if(codeEsti[u](nj) != code(nj)) BER_TOTAL_Nocode++;
        }
    #else
        for(int nj = 0; nj < LenBit; nj++){
            if(sourceEsti[u](nj) != source(nj)) BER_TOTAL_Code++;
        }
        for(int nj = 0; nj < NJ; nj++){
            if(codeEsti[u](nj) != code(nj)) BER_TOTAL_Nocode++;
        }
    #endif            
    }
}


double fmax(vector<double>& seq){
    int n = seq.size();
    double tmp = seq[0];
    double result = 0.0;
    
    for(int i = 0; i < n-1; ++i){
        result = max(tmp, seq[i+1]) + log(1 + exp(-abs(tmp-seq[i+1])));
        tmp = result;
    }
    return result;
}
