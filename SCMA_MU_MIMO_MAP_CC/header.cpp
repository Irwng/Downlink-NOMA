/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: include the basic functions of transmitter & recevier & channels
***********************************************************/
#include "header.h"

double N_Var;                                       /* variance of Noise */
double BER_TOTAL = 0;                               /* total number of error symbols */
double BER = 0;
double MaxExponent = 0.0;                           /* Max exponent keep the exp(-alpha/N_var)!=0 */
double MaxExponentUpper = 0.0;

SubConstellMatrix SubConstell;                      /* sub-constelltion point */
SymAfterMapMatrix MasterConstell[Mpoint];            /* Master-constelltion point */
SourceMatrix Source;                                /* source codewords */
CodeMatrix Code;                                    /* codewords after coding */
ModuMatrix Modu;                                    /* symbols after modulation */
MapMatrix F;                                        /* Mapping matrix */
SymAfterMapMatrix SymAfterMap[Nj];                  /* signals after mapping, Nj*Nt*M */
CSIMatrix H[U];                                     /* channel parameters , U*Nt*Nr */
CSIMatrix WeightedIdentityMatrix;                   /* MMSE assistance matrix */
PPMatrix V[U];                                      /* postprocessing matrix, Nj*Nr*M */
SymAfterBFMatrix SymAfterBF[Nj];                    /* receiving signals , Nj*Nr*M */
SymAfterFCMatrix SymAfterFC[Nj][U];
SymAfterPPMatrix SymAfterPP[Nj][U];                 /* receiving signals after postprocessing, Nj*Nr*M */
SymAfterMPAMatrix SymAfterMPA[Mod];
CodeMatrix DecodeBuff[U];
DecodeMatrix Decode[U];  


void NormalIO(){

#ifdef NoCode
    cout<<"NoCode"<<endl;
    outfile<<"NoCode"<<endl;
#else
    cout<<"Convolution Code & soft Viterbi decoder: rate=1/"<<ReciRate<<", block length="<<NJ<<endl;
    outfile<<"Convolution Code & soft Viterbi decoder: rate=1/"<<ReciRate<<", block length="<<NJ<<endl;
#endif
    cout<<"ML"<<endl;
    outfile<<"ML"<<endl;

    #ifdef CLST
        cout<<"CLST"<<endl;
        outfile<<"CLST"<<endl;
    #else
        cout<<"Latin"<<endl;
        outfile<<"Latin"<<endl;
    #endif

    #ifdef UDM
        cout<<"UDM"<<endl;
        outfile<<"UDM"<<endl;
    #else
        cout<<"PhaseRotation"<<endl;
        outfile<<"PhaseRotation"<<endl;
    #endif

    #ifdef FaddingChannel
        cout<<"FaddingChannel"<<endl;
        outfile<<"FaddingChannel"<<endl;
    #else
        cout<<"AWGN"<<endl;
        outfile<<"AWGN"<<endl;
    #endif
    cout<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"M: "<<M<<setw(10)<<"J: "<<J<<endl;   
    outfile<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"M: "<<M<<setw(10)<<"J: "<<J<<endl;   
    
    cout<<"EbN0dB"<<setw(15)<<"BER"<<endl;
    outfile<<"EbN0dB"<<setw(15)<<"BER"<<endl;
}


void InitMapMatrix(MapMatrix& map, SubConstellMatrix& subConstell, SymAfterMapMatrix* masterConstell){

    //initialize the mapping matrix
    #ifdef CLST
      /* Nt= 0 Nt= 1 */
    map<<p1,p0,p0,p0,
         p2,p1,p0,p0,
         p3,p2,p1,p0,
         p0,p3,p2,p1,
         p0,p0,p3,p2,
         p0,p0,p0,p3;

    #else//Latin
      /* Nt= 0 Nt= 1 */
    map<<p1,p1,p0,p0,
         p2,p0,p1,p0,
         p3,p0,p0,p1,
         p0,p2,p2,p0,
         p0,p3,p0,p2,
         p0,p0,p3,p3;
    // map<<p1,p0,p0,p1,
    //      p2,p0,p0,p2,
    //      p3,p0,p0,p3,
    //      p0,p1,p1,p0,
    //      p0,p2,p2,p0,
    //      p0,p3,p3,p0;
    #endif

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
    /* maybe can simplified by bitset instead of this complex translation */
    for(int i = 0; i < Mpoint; ++i){                                       /*ideal master constellation coding*/\
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

    // #ifdef DebugMode
    //     cout<<"...................MasterConstell..................."<<endl;
    //     for(int i = 0; i < Mpoint; ++i){
    //         cout<<masterConstell[i]<<endl;
    //     }
    // #endif  
}


/* key point */
void ChannelInitialize(int ebN0dB){

    double ebN0 = pow(10, (double)ebN0dB/10);

#ifdef NoCode
    double snr = (static_cast<double>(J)/(Nt*M)) * ebN0;
#else
    double snr = (static_cast<double>(J)/(Nt*M*ReciRate))*ebN0;
#endif

    double N0 = 1 / snr;
    N_Var = N0/2;
    MaxExponent = Alpha * N_Var; // 
    MaxExponentUpper = Beta * N_Var; // 
    WeightedIdentityMatrix << 2*N_Var, 0, 0, 2*N_Var;
    BER_TOTAL = 0;
    #ifdef DebugMode
        cout<<"N_Var: "<<N_Var<<endl;
    #endif
}


void BitSource(SourceMatrix& source){

    //random number generator: 0 or 1
    for(int i = 0; i < LenBit; i++){
        source(i)= rand()%2;
    }
    #ifdef DebugMode
        cout<<".................source................."<<endl;
        cout<<source<<endl;
    #endif
}


void ConvEncoder(SourceMatrix& source, CodeMatrix& code){

    //(3,1,3)
#ifndef NoCode 

        /* 5 7 7 */
    int reg1=0,reg2=0;
    for(int i = 0; i < LenBit; i++){
        code(ReciRate*i + 0) = source(i) xor reg2; 
        code(ReciRate*i + 1) = source(i) xor reg1 xor reg2;
        code(ReciRate*i + 2) = source(i) xor reg1 xor reg2;

        reg2 = reg1;
        reg1 = source(i);
    }

    // int reg1=0, reg2=0;
    // for(int i = 0; i < LenBit; i++){
    //     code(ReciRate*i + 0) = source(i);
    //     code(ReciRate*i + 1) = source(i) xor reg1 xor reg2;
    //     code(ReciRate*i + 2) = source(i) xor reg2;

    //     reg2 = reg1;
    //     reg1 = source(i);
    // }
#else //Repetition
    for(int j = 0; j < ReciRate; j++){
        for(int i = 0; i < LenBit; i++){
            code(LenBit*j + i) = source(i);
        }
    }
#endif

#ifdef DebugMode
    cout<<".................code................."<<endl;
    cout<<code<<endl;
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
        cout<<modu<<endl;
    #endif
}


void Mapper(ModuMatrix& modu, SymAfterMapMatrix* symaftermap, MapMatrix& map){

    //Codebook:F
    for(int nj = 0; nj < Nj; nj++){
        for(int nt = 0; nt < Nt; nt++){
            for(int m = 0; m < M; m++){
                symaftermap[nj](nt, m) = p0;
                for(int j = 0; j < J; j++){
                    symaftermap[nj](nt, m) += map(j, m + nt * Nt) * modu(nj, j);
                }
            }
        }
        #ifdef DebugMode
            if(nj == 0){
                cout<<"................symaftermap[0]................ "<<endl;
                cout<<symaftermap[0]<<endl; 
            }
        #endif
    }
}


void FadingChannel(CSIMatrix* h){

    //channel parameters
    for(int u = 0; u<U; ++u){
        for(int nr = 0; nr<Nr; nr++){
            for(int nt = 0; nt<Nt;nt++){
                double GG = sqrt(-2*log(randN() + 0.000000001));
                double B = randN();
                h[u](nr, nt) = ComplexD(sqrt(0.5) * GG * cos(2*PI*B),
                                        sqrt(0.5) * GG * sin(2*PI*B));
            }
        }
    }
    #ifdef DebugMode
        for(int u = 0; u<U; ++u){
            cout<<"................h["<<u<<"].................. "<<endl;
            cout<<h[U]<<endl; 
        }
    #endif
}


ComplexD AWGN(double nvar){

    double GG = sqrt(-2*log(randN() + 0.000000001));
    double B = randN();
    ComplexD GuassN(sqrt(nvar)* GG * cos(2*PI*B), sqrt(nvar)* GG * sin(2*PI*B));
    return GuassN;
}


void Receiver(SymAfterFCMatrix symAfterFC[Nj][U], 
              SymAfterBFMatrix* symafterbf, 
              CSIMatrix* h){

    for(int nj = 0; nj < Nj; nj++){
        for(int u = 0; u < U; u++){
            /* AWGN */
            SymAfterFCMatrix tmp;
            for(int nr = 0; nr < Nr; nr++){
                for(int m = 0; m < M; m++){
                    tmp(nr, m) = AWGN(N_Var);
                }
            }

            symAfterFC[nj][u] = h[u] * symafterbf[nj] + tmp;//

            #ifdef DebugMode
                if(nj == 0){
                    // cout<<".................AWGN................."<<endl;
                    // cout<<tmp<<endl; 
                    
                    cout<<".................symAfterFC["<<nj<<"]["<<u<<"]................."<<endl;
                    cout<<symAfterFC[nj][u]<<endl; 
                }    
            #endif
        }
    }
}


void DirectDecoder(SymAfterMPAMatrix* symAfterMPA, CodeMatrix* decodebuff, 
                   CodeMatrix& code){
    
    for(int u = 0; u < U; u++){        
        for(int nj = 0; nj < NJ; nj++){
            // decodebuff[u](nj) = symAfterMPA[1](u,nj) > 0.5 ? 1:0;
            decodebuff[u](nj) = symAfterMPA[1](u,nj) > symAfterMPA[0](u,nj) ? 1:0;
        }
    }
}


void Compare(CodeMatrix* decodebuff, CodeMatrix& code, 
             SourceMatrix& source, SourceMatrix* decode){

    for(int u = 0; u < U; u++){
        #ifdef NoCode
            for(int nj = 0; nj < NJ; nj++){
                if(decodebuff[u](nj) != code(nj)) BER_TOTAL++;
            }
        #else
            for(int nj = 0; nj < LenBit; nj++){
                if(decode[u](nj) != source(nj)) BER_TOTAL++;
            }
        #endif            
    }
}


void IterationHelper(CodeMatrix* decodebuff, SourceMatrix* decode, 
                     SymAfterFCMatrix symAfterFC[Nj][U], CSIMatrix* h,
                     SymAfterPPMatrix symAfterPP[Nj][U]){

    SymAfterMapMatrix tmpSymAfterMap[Nj];
    
    for(int u = 0; u < U; u++){
        #ifdef NoCode
            Modulation(decodebuff[u], Modu);
        #else
            ConvEncoder(decode[u], Code);
            Modulation(Code, Modu);
        #endif
        Mapper(Modu, tmpSymAfterMap, F);

        for(int nj = 0; nj < Nj; ++nj){
            for(int u = 0; u < U; ++u){
                for(int m = 0; m < M; m++){
                    symAfterPP[nj][u](0, m) = ( symAfterFC[nj][u](m) - h[u](1) * tmpSymAfterMap[nj](1, m) ) / h[u](0);
                    symAfterPP[nj][u](1, m) = ( symAfterFC[nj][u](m) - h[u](0) * tmpSymAfterMap[nj](0, m) ) / h[u](1);
                }
            }
        }
    }
}


void ML_Iter(SymAfterPPMatrix symAfterPP[Nj][U], SymAfterMPAMatrix* symAfterMPA, 
             SymAfterMapMatrix* masterConstell){

    int j = 0;    
    double NormalBase = 0.0;                     /* the base to normalize */   
    // double MaxDist = 0.0;                        /* the maximum Euclidean dist */
    // double MinDist = 0.0;                        /* the minimum Euclidean dist */
    // int MinDistPoint = 0;
    // double MaxDistSubMaxExpo = 0.0;              /* dmax - alpha x Alpha = dmax - MaxExponent*/
	double d[Mpoint] = {0.0};                    /* Euclidean dist of each star point*/   
    double f[Mpoint] = {0.0};                    /* prob of each star point */
    double Q_post[J][Mod] = {0.0};               /* post prob of each bit */

    for(int nj = 0; nj < Nj; ++nj){
        for(int u = 0; u < U; ++u){

            /* reshape the recevier signals */
            ComplexD rx[4] = {symAfterPP[nj][u](0, 0), symAfterPP[nj][u](0, 1), 
                              symAfterPP[nj][u](1, 0), symAfterPP[nj][u](1, 1)};
            
            /* initialize */    
            for(j = 0; j < J; j++)
                for(int m = 0; m < Mod; m++)
                    Q_post[j][m] = 0;
            
            // MaxDist = 0.0;
            // MinDist = INT16_MAX;
            // MinDistPoint = 0;

            /* calculate the Euclidean distances and fix them */
            for(int i = 0; i < Mpoint; i++){
                d[i] = 0;
                for(int nt = 0; nt < Nt; ++nt){
                    for(int m = 0; m < M; ++m){
                        d[i] += pow(abs(rx[nt * Nt + m] - masterConstell[i](nt, m)), 2);
                    }
                }
                // if(d[i] > MaxDist){
                //     MaxDist = d[i];
                // }
                // if(d[i] < MinDist){
                //     MinDist = d[i];
                    // MinDistPoint = i;
                // }
            }
            
            /* fix the overload data */
            // if(MaxDist >= MaxExponent){
            //     MaxDistSubMaxExpo = MaxDist - MaxExponent;
            //     if(MaxDistSubMaxExpo > (MaxExponentUpper + MinDist)){
            //         MaxDistSubMaxExpo = MaxExponentUpper + MinDist;
            //     }
            //     for(int i = 0; i < Mpoint; i++){
            //         d[i] -= MaxDistSubMaxExpo;
            //     }
            // }
            
            /* calculate the prob of each points and Marginal prob of each user */
            for(int i = 0; i < Mpoint; i++){
                bitset<6> bit(i);
                // f[i] = exp(- d[i]); 
                f[i] = exp(- d[i] / N_Var);  /* the format according to theory */ 
                for(j = 0; j < J; j++){
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

            // bitset<J> bits(MinDistPoint);

            // for(j = 0; j < J; j++){
            //     symAfterMPA(u, nj * J + j) = bits[j]; 
            // }      
        }
    }

    #ifdef DebugMode
        cout<<endl<<".................ML:SymAfterMPA................."<<endl;
        cout<<symAfterMPA<<endl;
    #endif       
}


