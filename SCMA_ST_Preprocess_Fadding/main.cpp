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
Precoding: None, do preprocess in the recevier
***********************************************************/
#include "header.h"

fstream outfile;

int main(int argc, char * argv[]){

    srand((unsigned)time(NULL));

    time_t start = time(NULL);
    
    if(argc < 4){
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
        
        cout<<"3rd argument: MPA or ML"<<endl;
        cout<<"1: ML"<<endl;
        cout<<"2: MPA"<<endl;

        cout<<"4th argument: Preprocess"<<endl;
        cout<<"1: ZF"<<endl;
        cout<<"2: MMSE"<<endl;
        cout<<"3: OSIC"<<endl;
          
        return 0;
    }
    
    Initialize(argv);
    InitMapMatrix(F, SubConstell, MasterConstell, argv);
    
    for(int EbN0dB = minEbN0dB; EbN0dB <= maxEbN0dB; EbN0dB = EbN0dB + step){
        ChannelInitialize(EbN0dB);
        
        for(int loop = 1; loop <= NLoop; loop++){
            BitSource(Source);
            ConvEncoder(Source, Code);
            Modulation(Code, Modu);
            Mapper(Modu, SymAfterMap, F);
            FadingChannel(H);
            Receiver(SymAfterFC, SymAfterMap, H, SymAfterPP, V, argv);
            
            switch (*argv[3]){
                case '1':
                    ML(SymAfterPP, SymAfterMPA, MasterConstell);
                    break;
                case '2':
                    MPA(SymAfterPP, SymAfterMPA, SubConstell, argv);
                    break;
                default:
                    break;
            }

            DirectDecoder(SymAfterMPA, DecodeBuff, Code);
            Compare(DecodeBuff, Code);
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
