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
    
    if(argc < 3){

        cout<<"1st argument: codebook structures"<<endl;
        cout<<"1: NOS"<<endl;
        cout<<"2: OS"<<endl;
        cout<<"3: TS"<<endl;
        
        cout<<"2nd argument: Preprocess"<<endl;
        cout<<"1: ZF"<<endl;
        cout<<"2: MMSE"<<endl;
        cout<<"3: OSIC"<<endl;
          
        return 0;
    }
    
    Initialize(argv);
    InitMapMatrix(F, MasterConstell, argv);
    
    for(int EbN0dB = MinEbN0dB; EbN0dB <= MaxEbN0dB; EbN0dB = EbN0dB + Step){
        ChannelInitialize(EbN0dB);
        
        for(int loop = 1; loop <= NLoop; loop++){
            BitSource(Source);
            Modulation(Source, Modu);
            Mapper(Modu, SymAfterMap, F);
            FadingChannel(H);
            Receiver(SymAfterFC, SymAfterMap, H, SymAfterPP, V, argv);
            ML(SymAfterPP, SymAfterMPA, MasterConstell);
            DirectDecoder(SymAfterMPA, Decode, Source);
            // Compare(Decode, Source);
        }

        BER = static_cast<double>(BER_TOTAL/(NLoop*U*NJ));
        cout<<EbN0dB<<setw(20)<<BER<<endl;
        outfile<<EbN0dB<<setw(20)<<BER<<endl;
        cout<<"time(s): "<<time(NULL) - start<<endl;
    }

    cout<<"time(s): "<<time(NULL) - start<<endl;
    outfile<<"time(s): "<<time(NULL) - start<<endl;
    outfile.close();

    return 0;
}
