/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: include the basic functions of transmitter & recevier & channels
***********************************************************/
#include "header.h"

double N_Var;                                       /* variance of Noise */
double BER_TOTAL = 0;                               /* total number of error symbols */
double BER = 0;

SymAfterMapMatrix MasterConstell[Mpoint];           /* Master-constelltion point */
SourceMatrix Source;                                /* source codewords */
ModuMatrix Modu;                                    /* symbols after modulation */
MapMatrix F;                                        /* Mapping matrix */
SymAfterMapMatrix SymAfterMap[Nj];                  /* signals after mapping, Nj*Nt*M */
CSIMatrix H[U];                                     /* channel parameters , U*Nt*Nr */
CSIMatrix WeightedIdentityMatrix;                   /* MMSE assistance matrix */
PPMatrix V[U];                                      /* postprocessing matrix, Nj*Nr*M */
SymAfterFCMatrix SymAfterFC[Nj][U];
SymAfterPPMatrix SymAfterPP[Nj][U];                 /* receiving signals after postprocessing, Nj*Nr*M */
SymAfterMPAMatrix SymAfterMPA[Mod];
SourceMatrix Decode[U];  

ComplexD p1;
ComplexD p2;
ComplexD p3;
ComplexD p4;
ComplexD p5;

void Initialize(char* argv[]){

    cout<<"SCMA_ST_CLST_Fadding"<<endl;
    outfile.open("SCMA_ST_CLST_Fadding.txt", ios::out|ios::app);
    outfile<<"....................................................."<<endl;
    outfile<<"SCMA_ST_CLST_Fadding"<<endl;

#ifdef CompareOfStructure
    cout<<"CompareOfStructure"<<endl;
    outfile<<"CompareOfStructure"<<endl;
#endif

#ifdef CompareOfBlock
    cout<<"CompareOfBlock"<<endl;
    cout<<"Don't use NOS!"<<endl;
    outfile<<"CompareOfBlock"<<endl;
#endif

#ifdef CompareOfTime
    cout<<"CompareOfTime"<<endl;
    cout<<"Don't use TS!"<<endl;
    outfile<<"CompareOfTime"<<endl;
#endif

    switch (*argv[2]){

        case '1':
            cout<<"ZF"<<endl;
            outfile<<"ZF"<<endl;
            break;

        case '2':
            cout<<"MMSE"<<endl;
            outfile<<"MMSE"<<endl;
            break;
        
        /* Ordered Successive Interference Cancellation */
        case '3':
            cout<<"OSIC"<<endl;
            outfile<<"OSIC"<<endl;
            break;

        default:
            break;
    }
}


void InitMapMatrix(MapMatrix& map, SymAfterMapMatrix* masterConstell, 
                   char* argv[]){
    
    double weight;

#ifdef CompareOfStructure
    /* initialize the mapping matrix's codewords */
    switch(*argv[1]){

        case '1':
            /* J = 8, NOS */
            cout<<"Non-overlapping section(NOS)"<<endl;
            outfile<<"Non-overlapping section(NOS)"<<endl;

            weight = sqrt(2);
            p1 = ComplexD(1/weight, 0/weight);
            p2 = ComplexD(0/weight, 1/weight);

            /*   Nt= 0 Nt= 1  */
            map<<p1,p0,p0,p0,
                 p2,p0,p0,p0,
                 p0,p1,p0,p0,
                 p0,p2,p0,p0,
                 p0,p0,p1,p0,
                 p0,p0,p2,p0,
                 p0,p0,p0,p1,
                 p0,p0,p0,p2;
            break;
        
        case '2':
            /* J = 8, OS */
            cout<<"Overlapping section(OS)"<<endl;
            outfile<<"Overlapping section(OS)"<<endl;

            weight = sqrt(26);
            p1 = ComplexD(1/weight, 0/weight);
            p2 = ComplexD(0/weight, 1/weight);
            p3 = ComplexD(2/weight, 0/weight);
            p4 = ComplexD(0/weight, 2/weight);
            p5 = ComplexD(4/weight, 0/weight);

            /*   Nt= 0 Nt= 1  */
            map<<p1,p0,p0,p0,
                 p2,p1,p0,p0,
                 p3,p2,p1,p0,
                 p4,p3,p2,p1,
                 p5,p4,p3,p2,
                 p0,p5,p4,p3,
                 p0,p0,p5,p4,
                 p0,p0,p0,p5;
            break;

        case '3':
            /* J = 8, TS */
            cout<<"Tailbiting section(TS)"<<endl;
            outfile<<"Tailbiting section(TS)"<<endl;

            weight = sqrt(10);
            p1 = ComplexD(1/weight, 0/weight);
            p2 = ComplexD(0/weight, 1/weight);
            p3 = ComplexD(2/weight, 0/weight);
            p4 = ComplexD(0/weight, 2/weight);

            /*   Nt= 0 Nt= 1  */
            map<<p1,p0,p0,p3,
                 p2,p0,p0,p4,
                 p3,p1,p0,p0,
                 p4,p2,p0,p0,
                 p0,p3,p1,p0,
                 p0,p4,p2,p0,
                 p0,p0,p3,p1,
                 p0,p0,p4,p2;
            break;
        
        default:
            break;
    }
#endif

#ifdef CompareOfBlock
    /* initialize the mapping matrix's codewords */
    switch(*argv[1]){

        case '1':{
            /* J = 6, NOS */
            cout<<"Useless!"<<endl;
            exit(0);
        }
        
        case '2':
            /* J = 6, OS */
            cout<<"Overlapping section(OS)"<<endl;
            outfile<<"Overlapping section(OS)"<<endl;

            weight = sqrt(6);
            p1 = ComplexD(1/weight, 0/weight);
            p2 = ComplexD(0/weight, 1/weight);
            p3 = ComplexD(2/weight, 0/weight);

            /*   Nt= 0 Nt= 1  */
            map<<p1,p0,p0,p0,
                 p2,p1,p0,p0,
                 p3,p2,p1,p0,
                 p0,p3,p2,p1,
                 p0,p0,p3,p2,
                 p0,p0,p0,p3;
            break;

        case '3':
            /* J = 6, TS */
            cout<<"Tailbiting section(TS)"<<endl;
            outfile<<"Tailbiting section(TS)"<<endl;

            weight = sqrt(6);
            p1 = ComplexD(1/weight, 0/weight);
            p2 = ComplexD(0/weight, 1/weight);
            p3 = ComplexD(2/weight, 1/weight);

            /*   Nt= 0 Nt= 1  */
            map<<p1,p0,p1,p0,
                 p2,p0,p2,p0,
                 p3,p0,p3,p0,
                 p0,p1,p0,p1,
                 p0,p2,p0,p2,
                 p0,p3,p0,p3;
            break;
        
        default:
            break;
    }
#endif

#ifdef CompareOfTime
    /* initialize the mapping matrix's codewords */
    switch(*argv[1]){

        case '1':{
            /* J = 8, NOS */
            cout<<"Non-overlapping section(NOS)"<<endl;
            outfile<<"Non-overlapping section(NOS)"<<endl;

            weight = sqrt(10);
            p1 = ComplexD(1/weight, 0/weight);
            p2 = ComplexD(0/weight, 2/weight);
            p3 = ComplexD(1/weight, 1/weight);
            p4 = ComplexD(2/weight, 0/weight);
            
            /*   Nt= 0 Nt= 1  */
            map<<p1,p0,p0,p0,
                 p2,p0,p0,p0,
                 p3,p0,p0,p0,
                 p4,p0,p0,p0,
                 p0,p0,p1,p0,
                 p0,p0,p2,p0,
                 p0,p0,p3,p0,
                 p0,p0,p4,p0;
            break;
        }
        
        case '2':
            /* J = 8, OS */
            cout<<"Overlapping section(OS)"<<endl;
            outfile<<"Overlapping section(OS)"<<endl;

            weight = sqrt(26);
            p1 = ComplexD(1/weight, 0/weight);
            p2 = ComplexD(0/weight, 1/weight);
            p3 = ComplexD(2/weight, 0/weight);
            p4 = ComplexD(0/weight, 2/weight);
            p5 = ComplexD(4/weight, 0/weight);
            
            /*   Nt= 0 Nt= 1  */
            map<<p1,p0,p0,p0,
                 p2,p0,p0,p0,
                 p3,p0,p0,p0,
                 p4,p0,p1,p0,
                 p5,p0,p2,p0,
                 p0,p0,p3,p0,
                 p0,p0,p4,p0,
                 p0,p0,p5,p0;
            break;

        case '3':
            /* J = 8, TS */
            cout<<"Useless!"<<endl;
            exit(0);
        
        default:
            break;
    }
#endif

    cout<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"M: "<<M<<setw(10)<<"J: "<<J<<setw(10)<<"NLoop: "<<NLoop<<endl;   
    outfile<<"Nt: "<<Nt<<setw(10)<<"Nr: "<<Nr<<setw(10)<<"M: "<<M<<setw(10)<<"J: "<<J<<setw(10)<<"NLoop: "<<NLoop<<endl;   
    
    cout<<"EbN0dB"<<setw(15)<<"BER"<<endl;
    outfile<<"EbN0dB"<<setw(15)<<"BER"<<endl;

    /* initialize the Masterconstellation */
#ifdef DebugMode
    cout<<"MasterConstell"<<endl;
#endif
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
    #ifdef DebugMode
        cout<<"...... i: "<<i<<" ........."<<endl<<masterConstell[i]<<endl;
    #endif
    }
}


/* key point */
void ChannelInitialize(int ebN0dB){

    double ebN0 = pow(10, (double)ebN0dB/10);
    double snr = BitperSymbol * ebN0;
    double N0 = Nt * power / snr;
    N_Var = N0/2;
    BER_TOTAL = 0;
    for(int nt = 0; nt < Nt; ++nt){
        WeightedIdentityMatrix(nt, nt) = ComplexD(N_Var, 0);
    }
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


void Modulation(SourceMatrix& source, ModuMatrix& modu){

    /* BPSK: 0->-1, 1->1 */
    for(int nj = 0; nj < Nj; nj++){
        for(int j = 0; j < J; j++){
            modu(nj, j) = 2 * source(nj * J + j) - 1;
        }
    }
#ifdef DebugMode
    cout<<".................bpsk signal................."<<endl;
    cout<<modu<<endl;
#endif
}


void Mapper(ModuMatrix& modu, SymAfterMapMatrix* symAfterMap, MapMatrix& map){

    for(int nj = 0; nj < Nj; nj++){
        for(int m = 0; m < M; m++){
            for(int nt = 0; nt < Nt; nt++){
                symAfterMap[nj](nt, m) = p0;
                for(int j = 0; j < J; j++){
                    symAfterMap[nj](nt, m) +=  map(j, m + nt * Nt) * modu(nj, j);
                }
            }
        }
#ifdef DebugMode
    if(nj == 0){
        cout<<"................symAfterMap[0]................ "<<endl;
        cout<<symAfterMap[0]<<endl; 
    }
#endif
    }
}


void FadingChannel(CSIMatrix* h){

    /* channel parameters */
    for(int u = 0;u<U;u++){
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
    cout<<"................h[0].................. "<<endl;
    cout<<h[0]<<endl; 
#endif
}


ComplexD AWGN(double nvar){

    double GG = sqrt(-2*log(randN() + 0.000000001));
    double B = randN();
    ComplexD GuassN(sqrt(nvar)* GG * cos(2*PI*B), sqrt(nvar)* GG * sin(2*PI*B));
    return GuassN;
}


void Receiver(SymAfterFCMatrix symAfterFC[Nj][U], SymAfterMapMatrix* symAfterMap, 
              CSIMatrix* h, SymAfterPPMatrix symAfterPP[Nj][U], 
              PPMatrix* v, char* argv[]){

    CSIMatrix ConjTrans_H;
    CSIMatrix Denomiator_H;

    for(int u = 0; u < U; u++){
        ConjTrans_H = h[u].conjugate().transpose();
        
        switch (*argv[2]){
            case '1':
                Denomiator_H = ConjTrans_H * h[u];
                break;
            
            case '2':
                Denomiator_H = ConjTrans_H * h[u] + WeightedIdentityMatrix;
                break;

            default:
                break;
        }
        v[u] = Denomiator_H.inverse() * ConjTrans_H;
        
        for(int nj = 0; nj < Nj; nj++){
            /* AWGN */
            SymAfterFCMatrix tmp;
            for(int nr = 0; nr < Nr; nr++){
                for(int m = 0; m < M; m++){
                    tmp(nr, m) = AWGN(N_Var);
                }
            }  
            symAfterFC[nj][u] = h[u] * symAfterMap[nj] + tmp;//
            symAfterPP[nj][u] = v[u] * symAfterFC[nj][u];
            #ifdef DebugMode
                if(u==1&&nj==0) cout<<symAfterPP[0][0]<<endl;
            #endif
        }
    }
}


void DirectDecoder(SymAfterMPAMatrix* symAfterMPA, SourceMatrix* decode, 
                   SourceMatrix& source){
    
    for(int u = 0; u < U; u++){        
        for(int nj = 0; nj < NJ; nj++){
            decode[u](nj) = symAfterMPA[1](u,nj) > symAfterMPA[0](u,nj) ? 1:0;
            BER_TOTAL += (decode[u](nj) != source(nj));
        }
        #ifdef DebugMode
            cout<<"decode of user: "<<u<<endl;
            for(int nj = 0; nj < NJ; nj++){
                cout<<decode[u](nj)<<" ";
            }
            cout<<endl;
        #endif
    }
}


void Compare(SourceMatrix* decode, SourceMatrix& source){

    for(int u = 0; u < U; u++){
        for(int nj = 0; nj < NJ; nj++){
            if(decode[u](nj) != source(nj)) BER_TOTAL++;
        }
    }
}
