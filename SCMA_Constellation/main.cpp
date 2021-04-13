/************************************************************
Author: Wangyi     Version: C++11  Date: 2021/2/20
Theme: Calculate the Euclidean distance of constellation
Moduulaiton: BPSK
Antenna: Nr=1, Nt=2
Time slot(M): 2
***********************************************************/
#include "header.h"

int main(int argc, char * argv[]){

    srand((unsigned)time(NULL));

    if(argc < 3){
        cout<<"1st argument: codewords"<<endl;
        cout<<"1: UDM-202,2,4,9"<<endl;
        cout<<"2: UDM-42,1,2,4"<<endl;
        cout<<"3: UDM-178,3,4,8"<<endl;
        cout<<"4: UDM-1,1+i, 2i"<<endl;
        cout<<"6: PhaseRotation-pi/4"<<endl;
        cout<<"7: PhaseRotation-pi/5"<<endl;
        cout<<"8: PhaseRotation-pi/6"<<endl;

        cout<<"2nd argument: codebook structure"<<endl;
        cout<<"1: Latin"<<endl;
        cout<<"2: CLST"<<endl;
        cout<<"3: Orthogonal"<<endl;
        cout<<"4: half-Latin"<<endl;

        return 0;
    }

    Initialize();
    InitMapMatrix(F, SubConstell, MasterConstell, argv);

    return 0;
}
