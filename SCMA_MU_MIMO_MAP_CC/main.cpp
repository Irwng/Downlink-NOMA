/************************************************************
Author: Wangyi     Version: C++11  Date: 2021/2/21
Theme: UDM Downlink Multi-user MIMO
Channel: fading channel
Moduulaiton: BPSK
User number(U): 6 
Antenna: Nr=1, Nt=2
Decoding: Joint-MPA
Time slot(M): 2
Note: MAP
***********************************************************/
#include "header.h"

fstream outfile;

int main(){

    srand((unsigned)time(NULL));

#ifndef DebugMode    

    time_t start = time(NULL);
    
    constexpr long NLoop = pow(10, 4);          /* number of simulation loops  */
    constexpr int minEbN0dB = 0;
    
    #ifdef FaddingChannel
        constexpr int maxEbN0dB = 30;
    #else
        constexpr int maxEbN0dB = 10;
    #endif
    
    cout<<"SCMA_MU_MIMO_MAP_CC"<<endl;
    outfile.open("SCMA_MU_MIMO_MAP_CC.txt", ios::out|ios::app);
    outfile<<"....................................................."<<endl;
    outfile<<"SCMA_MU_MIMO_MAP_CC"<<endl;
    
    cout<<"Length of bits: "<<LenBit<<setw(20)<<"Rate: "<<ReciRate<<setw(20)<<"NLoop: "<<NLoop<<endl;   
    outfile<<"Length of bits: "<<LenBit<<setw(20)<<"Rate: "<<ReciRate<<setw(20)<<"NLoop: "<<NLoop<<endl; 

    NormalIO();
#else
    
    constexpr long NLoop = 1;                /* number of simulation loops  */
    constexpr int minEbN0dB = 30;
    constexpr int maxEbN0dB = minEbN0dB;
#endif //DebugMode

    InitMapMatrix(F, SubConstell, MasterConstell);
    for(int EbN0dB = minEbN0dB; EbN0dB <= maxEbN0dB; EbN0dB = EbN0dB + 2){
        ChannelInitialize(EbN0dB);
        
        for(int loop = 1; loop <= NLoop; loop++){
            BitSource(Source);
            ConvEncoder(Source, Code);
            Modulation(Code, Modu);
            Mapper(Modu, SymAfterMap, F);
            FadingChannel(H);
            Receiver(SymAfterFC, SymAfterMap, H);
            #ifdef MPAMode
                MPA(SymAfterPP, SymAfterMPA, SubConstell);
            #else
                ML(SymAfterFC, SymAfterMPA, MasterConstell, H, F);
            #endif
            
            for(int iter = 0; iter <= Iter; iter++){

                #ifdef NoCode
                    DirectDecoder(SymAfterMPA, DecodeBuff, Code);
                #else
                    #ifdef Soft_Prob
                        ViterbiSoftDecoder_Prob(SymAfterMPA, Source, Decode, Code);
                    #else
                        ViterbiHardDecoder(SymAfterMPA, Source, Decode, Code);
                    #endif
                #endif

                if(iter < Iter){
                    #ifdef DebugMode
                        cout<<"................iter: "<<iter<<"............"<<endl;
                    #endif
                    IterationHelper(DecodeBuff, Decode, SymAfterFC, H, SymAfterPP);
                    /* ML like in 2x2-MU-MIMO */
                    ML_Iter(SymAfterPP, SymAfterMPA, MasterConstell);
                }
            }
            Compare(DecodeBuff, Code, Source, Decode);
            /* process bar, one # means 5% */
            if(loop*20 % NLoop == 0) cout<<"#"<<flush;
            if(loop == NLoop) cout<<endl;
        }

        #ifndef DebugMode
            #ifdef NoCode
                BER = static_cast<double>(BER_TOTAL/(NLoop*U*NJ));
            #else
                BER = static_cast<double>(BER_TOTAL/(NLoop*U*LenBit));
            #endif
            cout<<EbN0dB<<setw(23)<<BER<<endl;
            outfile<<EbN0dB<<setw(23)<<BER<<endl;
            cout<<"time(s): "<<time(NULL) - start<<endl;
        #else
            BER = static_cast<double>(BER_TOTAL/(NLoop*U*LenBit));
            cout<<EbN0dB<<setw(23)<<BER;
            cout<<endl;
        #endif //DebugMode
    }

    #ifndef DebugMode
        cout<<"time(s): "<<time(NULL) - start<<endl;
        outfile<<"time(s): "<<time(NULL) - start<<endl;
        outfile.close();
    #endif //DebugMode

    return 0;
}
