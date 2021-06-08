/************************************************************
Author: Wangyi     Version: C++11  Date: 2021/3/15
Theme:  SCMA in AWGN channel
Channel: fading channel
Moduulaiton: BPSK
User number(U): 6 
Antenna: Nr=Nt=2
Time slot(M): 2
Subcarrier: equal to Nt*M 
Decoding: Joint-MPA
Comparison:J, codebook structure, M
***********************************************************/
#include "header.h"

fstream outfile;

int main(int argc, char * argv[]){

    srand((unsigned)time(NULL));

    time_t start = time(NULL);
    
    if(argc < 5){
        cout<<"1st argument: codebook codewords"<<endl;
        cout<<"1: UDM-202,2,4,9"<<endl;
        cout<<"2: UDM-42,1,2,4"<<endl;
        cout<<"3: UDM-178,3,4,8"<<endl;
        cout<<"4: UDM-1,1+i, 2i"<<endl;
        cout<<"6: PhaseRotation-pi/4"<<endl;
        cout<<"7: PhaseRotation-pi/5"<<endl;
        cout<<"8: PhaseRotation-pi/6"<<endl;

        cout<<"2nd argument: codebook structures"<<endl;
        cout<<"1: Latin"<<endl;
        cout<<"2: CLST"<<endl;
        cout<<"3: Orthogonal"<<endl;
        cout<<"4: half-Latin"<<endl;
        
        cout<<"3rd argument: Process"<<endl;
        cout<<"1: ZF"<<endl;
        cout<<"2: MMSE"<<endl;
        
        cout<<"4rd argument: OSIC"<<endl;
        cout<<"1: SINR"<<endl;
        cout<<"2: SNR"<<endl;
          
        return 0;
    }
    
    Initialize(argv);
    InitMapMatrix(F, SubConstell, MasterConstell, argv);
    
    for(int EbN0dB = MinEbN0dB; EbN0dB <= MaxEbN0dB; EbN0dB = EbN0dB + Step){
        ChannelInitialize(EbN0dB);
        
        for(int loop = 1; loop <= NLoop; loop++){
            BitSource(Source);
            ConvEncoder(Source, Code);
            Modulation(Code, Modu);
            Mapper(Modu, SymAfterMap, F);
            FadingChannel(H);
            #ifdef OSIC
                Receiver_OSIC(SymAfterFC, SymAfterMap, H, SymAfterPP, V, SubConstell, argv);
            #else
                Receiver(SymAfterFC, SymAfterMap, H, SymAfterPP, V, argv);
            #endif
            
            #ifdef LL
                MLLL(SymAfterPP, SymAfterMPA, MasterConstell);
                ViterbiSoftDecoderLL(SymAfterMPA, SourceEsti, CodeEsti);
            #else
                ML(SymAfterPP, SymAfterMPA, MasterConstell);
                ViterbiSoftDecoder(SymAfterMPA, SourceEsti, CodeEsti);
            #endif

            Compare(CodeEsti, Code, SourceEsti, Source);
            /* process bar, one # means 5% */
            if(loop*20 % NLoop == 0) cout<<"#"<<flush;
            if(loop == NLoop) cout<<endl;
        }

        BER_Nocode = static_cast<double>(BER_TOTAL_Nocode/NUM_Nocode);
        BER_Code = static_cast<double>(BER_TOTAL_Code/NUM_Code);
        cout<<EbN0dB<<setw(20)<<BER_Nocode<<setw(20)<<BER_Code<<endl;
        outfile<<EbN0dB<<setw(20)<<BER_Nocode<<setw(20)<<BER_Code<<endl;
        cout<<"time(s): "<<time(NULL) - start<<endl;
    }

    cout<<"time(s): "<<time(NULL) - start<<endl;
    outfile<<"time(s): "<<time(NULL) - start<<endl;
    outfile.close();

    return 0;
}
