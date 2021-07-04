/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: include the basic functions of transmitter & recevier & channels
***********************************************************/
#include "header.h"

double N_Var;                                  /* variance of Noise */
double capacity_Total;                         /* total number of error bits without code */
double capacity;                               /* BER without code */
CSIMatrix H[U];                                /* channel parameters , U*Nt*Nr */
CSIMatrix WeightedIdentityMatrix;              /* MMSE assistance matrix */
ModuMatrix vectorS;

void Initialize(){

    cout<<"SCMA_ST_Fading_Latin_Capacity"<<endl;
    outfile.open("SCMA_ST_Fading_Latin_Capacity.txt", ios::out|ios::app);
    outfile<<"....................................................."<<endl;
    outfile<<"SCMA_ST_Fading_Latin_Capacity"<<endl;   
}


/* key point */
void ChannelInitialize(int ebN0dB){

    double ebN0 = pow(10, (double)ebN0dB/10);
    double snr = (static_cast<double>(U)/(Nt*M)) * ebN0;
    double N0 = power*Nt / snr;
    N_Var = N0;
    #ifdef DebugMode
        cout<<"N_Var: "<<N_Var<<endl;
    #endif
}


void FadingChannel(CSIMatrix* h){

    /* channel parameters */
    for(int u = 0; u < U; ++u){
        for(int nr = 0; nr < Nr; nr++){
            for(int nt = 0; nt < Nt;nt++){
                double GG = sqrt(-2*log(randN() + 0.000000001));
                double B = randN();
                h[u](nr, nt) = ComplexD(sqrt(0.5) * GG * cos(2*PI*B),
                                        sqrt(0.5) * GG * sin(2*PI*B));
            }
        }
    }
}


ComplexD AWGN(double nvar){

    double GG = sqrt(-2*log(randN() + 0.000000001));
    double B = randN();
    ComplexD GuassN(sqrt(nvar)* GG * cos(2*PI*B), sqrt(nvar)* GG * sin(2*PI*B));
    return GuassN;
}


void Capacity(CSIMatrix* h, ModuMatrix& vectorS){
    for(int u = 0; u < U; ++u){
        SVD(h[u], vectorS);
        double tmp = 0;
        for(int nt = 0; nt < Nt; ++nt)
            tmp += pow(vectorS(nt),2);
        capacity_Total += log(1 + 0.125*tmp/N_Var);
    }
    
}

void SVD(CSIMatrix& h, ModuMatrix& vectorS){

    JacobiSVD<CSIMatrix> svd(h, ComputeFullU | ComputeFullV);
    // CSIMatrix u = svd.matrixU();
    // CSIMatrix v = svd.matrixV();
    vectorS = svd.singularValues();

    #ifdef DebugMode
        CSIMatrix S =  u.inverse() * h * (v.conjugate().transpose().inverse()); // S = U^-1 * H * VT^-1
        cout << "---------------------JacobiSVD------------------------"<< std::endl;
        cout << "H"<< endl << h << endl;
        cout << "S"<< endl << S << endl;
        cout << "U*S*V^H"<< endl << u*S*(v.conjugate().transpose()) << endl;
    #endif
}