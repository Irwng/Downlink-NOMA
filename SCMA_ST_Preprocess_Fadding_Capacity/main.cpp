/************************************************************
Author: Wangyi     Version: C++11  Date: 2021/6/18
Theme:  Capacity of SCMA in Space-time strcuture
Channel: fading channel
Moduulaiton: BPSK
User number(U): 6 
Antenna: Nr=Nt=2
Time slot(M): 2
Decoding: ML
***********************************************************/
#include "header.h"

fstream outfile;

int main(int argc, char * argv[]){

    srand((unsigned)time(NULL));

    time_t start = time(NULL);
    
    Initialize();
    
    for(int EbN0dB = MinEbN0dB; EbN0dB <= MaxEbN0dB; EbN0dB = EbN0dB + Step){
        ChannelInitialize(EbN0dB);
        
        for(int loop = 1; loop <= NLoop; loop++){
            FadingChannel(H);
            Capacity(H,vectorS);
            
            /* process bar, one # means 5% */
            if(loop*20 % NLoop == 0) cout<<"#"<<flush;
            if(loop == NLoop) cout<<endl;
        }

        capacity = static_cast<double>(capacity_Total/NUM);
        cout<<EbN0dB<<setw(20)<<capacity<<endl;
        outfile<<EbN0dB<<setw(20)<<capacity<<endl;
        cout<<"time(s): "<<time(NULL) - start<<endl;
    }

    cout<<"time(s): "<<time(NULL) - start<<endl;
    outfile<<"time(s): "<<time(NULL) - start<<endl;
    outfile.close();

    return 0;
}
