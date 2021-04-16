/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: include the basic functions of transmitter & recevier & channels
***********************************************************/
#include "header.h"

double N_Var;                                       /* variance of Noise */
double BER_TOTAL = 0;                               /* total number of error symbols */
double BER = 0;

SubConstellMatrix SubConstell;                      /* sub-constelltion point */
SymAfterMapMatrix MasterConstell[Mpoint];           /* Master-constelltion point */
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

ComplexD p1;
ComplexD p2;
ComplexD p3;


void Initialize(char* argv[]){

    cout<<"SCMA_ST_Fadding"<<endl;
    outfile.open("SCMA_ST_Fadding.txt", ios::out|ios::app);
    outfile<<"....................................................."<<endl;
    outfile<<"SCMA_ST_Fadding"<<endl;

    switch (*argv[3]){
        case '1':
            cout<<"ML"<<endl;
            outfile<<"ML"<<endl;
            break;
        
        case '2':
            cout<<"MPA"<<endl;
            outfile<<"MPA"<<endl;
            break;

        default:
            break;
    }

    switch (*argv[4]){
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
            /* 178, 3, 4, 8 */
            cout<<"UDM-1,2i,1+1i"<<endl;
            outfile<<"UDM-1,2i,1+1i"<<endl;

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
    
    cout<<"EbN0dB"<<setw(15)<<"BER"<<endl;
    outfile<<"EbN0dB"<<setw(15)<<"BER"<<endl;

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
}


/* key point */
void ChannelInitialize(int ebN0dB){

    double ebN0 = pow(10, (double)ebN0dB/10);
    double snr = 1.5 * ebN0;
    double N0 = power*Nt / snr;
    N_Var = N0/2;
    for(int nt = 0; nt < Nt; ++nt){
        WeightedIdentityMatrix(nt, nt) = ComplexD(N_Var, 0);
    }
    BER_TOTAL = 0;
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

    for(int j = 0; j < ReciRate; j++){
        for(int i = 0; i < LenBit; i++){
            code(LenBit*j + i) = source(i);
        }
    }
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
            cout<<symaftermap[0]<<endl; 
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


void Receiver(SymAfterFCMatrix symAfterFC[Nj][U], SymAfterBFMatrix* symAfterBF, CSIMatrix* h,
              SymAfterPPMatrix symAfterPP[Nj][U], PPMatrix* v, char* argv[]){

    CSIMatrix ConjTrans_H;
    CSIMatrix Denomiator_H;

    for(int nj = 0; nj < Nj; nj++){
        for(int u = 0; u < U; u++){
            for(int u = 0; u < U; u++){
                /* AWGN */
                SymAfterFCMatrix tmp;
                for(int nr = 0; nr < Nr; nr++){
                    for(int m = 0; m < M; m++){
                        tmp(nr, m) = AWGN(N_Var);
                    }
                }  

                ConjTrans_H = h[u].conjugate().transpose();
                
                switch (*argv[4]){
                    case '1':
                        Denomiator_H = ConjTrans_H * h[u];
                        break;
                    
                    case '2':
                        Denomiator_H = ConjTrans_H * h[u] + WeightedIdentityMatrix;
                        break;

                    case '3':
                        break;

                    default:
                        break;
                }

                v[u] = Denomiator_H.inverse() * ConjTrans_H;
                symAfterFC[nj][u] = h[u] * symAfterBF[nj] + tmp;
                symAfterPP[nj][u] = v[u] * symAfterFC[nj][u];
            }
        }
    }
}


void DirectDecoder(SymAfterMPAMatrix* symAfterMPA, CodeMatrix* decodebuff, 
                   CodeMatrix& code){
    
    for(int u = 0; u < U; u++){        
        for(int nj = 0; nj < NJ; nj++){
            decodebuff[u](nj) = symAfterMPA[1](u,nj) > symAfterMPA[0](u,nj) ? 1:0;
        }
    }
}


void Compare(CodeMatrix* decodebuff, CodeMatrix& code){

    for(int u = 0; u < U; u++){
        for(int nj = 0; nj < NJ; nj++){
            if(decodebuff[u](nj) != code(nj)) BER_TOTAL++;
        }
    }
}
