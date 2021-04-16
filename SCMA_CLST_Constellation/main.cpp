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
    Initialize();
    InitMapMatrix(F, SubConstell, MasterConstell);

    return 0;
}
