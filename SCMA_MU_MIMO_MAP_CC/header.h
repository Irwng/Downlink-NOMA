/***********************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: include the declaration of constants、variables、functions，
             and the definition of new type
************************************************************************/
#ifndef SCMA_MU_MIMO_MAP_CC_HEADER_H
#define SCMA_MU_MIMO_MAP_CC_HEADER_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <vector>
#include <iomanip>
#include <complex>
#include <bitset>
#include <Eigen/SVD>
#include <Eigen/Dense> 
using namespace std;
using namespace Eigen;

#define randN() (rand()/(double)RAND_MAX)         /* Random value in [0,30000] */
typedef complex<double> ComplexD;

/**********************************************************************************
 * Control Center: must choose one and only one of the choices in different blocks
 **********************************************************************************/

/* running mode*/
// #define DebugMode
#define MonteCarlo

/* iterate or not */
// #define IterMode
#define NoiterMode

/* code of not */
// #define NoCode
#define ConvCode

/* MPA or ML */
// #define MPAMode
#define MLMode

/* mapping matrix structure: Latin or CLST */
#define Latin
// #define CLST

/* mapping matrix codewords: UDM or PhaseRotation */
// #define UDM
#define PhaseRotation

/* channel */
#define FaddingChannel
// #define AWGNMode

/* viterbi decoder */
#define Soft_Prob
// #define Hard

/**********************************
 * basic constant model parameters
 **********************************/

constexpr int U = 6;                              /* Number of users */
#ifdef NoCode
    constexpr int K = 1;                          /* Length of per user's bit stream,length of vector(b)=K*U */
#else
    constexpr int K = 50;
#endif
constexpr int LenBit = U * K;                     /* number of bits of all users */
constexpr int ReciRate = 3;                       /* reciprocal of the rate */
constexpr int NODE = 4;
constexpr int Route = NODE;                       /* number of the routes to keep in viterbi decoder */
constexpr int NJ = LenBit * ReciRate;             /* number of coded bits, J*Nj = K*U*ReciRate one block per transmissi  on  */
constexpr int J = 6;                              /* bits per block, length of vector(c) = J*Nj */
constexpr int Nj = NJ / J;                        /* number of blocks,J*Nj = K*U one block per transmission  */

constexpr int M = 2;                              /* Number of time slots resources */
constexpr int Nt = 2;                             /* number of antennas at transmitter */
constexpr int Nr = 1;                             /* number of antennas at recevier */
constexpr int Mod = 2;				              /* BPSK modulation order */
constexpr int Mpoint = 64;                        /* Master-constellation Points */
constexpr int Spoint = 8;                         /* Sub-constellation Points */
constexpr double PI = 3.141592653589793;
constexpr int Du = 2;                             /* Number of resources connected to user */
constexpr int Dr = 3;                             /* Number of users connected to resources */
/* -745.13~709.78 */
constexpr double Alpha = 745;                     /* Max exponent keep the exp(-alpha)!=0 */
constexpr double Beta = 709;                      /* Max exponent keep the exp(alpha)!=INF */
// constexpr double varepsilon = pow(10, -10);

#ifdef IterMode
    constexpr int Iter = 1;
#else
    constexpr int Iter = 0;
#endif//IterMode

/*UDM or PhaseRotation*/
#ifdef UDM
    constexpr double weight = sqrt(42);
    constexpr ComplexD p0(0, 0);
    constexpr ComplexD p1(1/weight, 4/weight);
    constexpr ComplexD p2(4/weight, 1/weight);
    constexpr ComplexD p3(2/weight, 2/weight);
#else
    constexpr double weight = sqrt(3);
    constexpr double RF = PI/5;
    constexpr ComplexD p0(0, 0);
    constexpr ComplexD p1(cos(0*RF)/weight, sin(0*RF)/weight);
    constexpr ComplexD p2(cos(1*RF)/weight, sin(1*RF)/weight);
    constexpr ComplexD p3(cos(2*RF)/weight, sin(2*RF)/weight);
#endif

/***********************************************************
 * basic type defination and golbal variables in matrix type
 ***********************************************************/

/* sub-constelltion point */
typedef Matrix<ComplexD, Nr*M, Spoint> SubConstellMatrix;
extern SubConstellMatrix SubConstell;

/* source codewords */
typedef Matrix<int, 1, LenBit> SourceMatrix;
extern SourceMatrix Source;

/* codewords after coding */
typedef Matrix<int, 1, NJ> CodeMatrix;
extern CodeMatrix Code;

/* symbols after modulation */
typedef Matrix<ComplexD, Nj, J> ModuMatrix;
extern ModuMatrix Modu;

/* Mapping matrix*/
typedef Matrix<ComplexD, J, M*Nt> MapMatrix;
extern MapMatrix F;

/* signals after mapping, Nj*Nt*M */
typedef Matrix<ComplexD, Nt, M> SymAfterMapMatrix;
extern SymAfterMapMatrix SymAfterMap[Nj];
extern SymAfterMapMatrix MasterConstell[Mpoint];


/* channel parameters , U*Nr*Nt */
typedef Matrix<ComplexD, Nr, Nt> CSIMatrix;
extern CSIMatrix H[U];

/* MMSE assistance matrix */
extern CSIMatrix WeightedIdentityMatrix;

/* postprocessing matrix, Nj*Nr*M */
typedef Matrix<ComplexD, Nt, Nr> PPMatrix;
extern PPMatrix V[U];

/* signals after Beamforming, Nj*Nt*M */
typedef SymAfterMapMatrix SymAfterBFMatrix;
extern SymAfterBFMatrix SymAfterBF[Nj];

/* signals after fadding channal, Nj*U*Nr*M */
typedef Matrix<ComplexD, Nr, M> SymAfterFCMatrix;
extern SymAfterFCMatrix SymAfterFC[Nj][U];

/* signals after post processing, Nj*U*Nt*M */
/* After the post processing, the dimensionality of the signal matrix turns to Nt*M */
typedef Matrix<ComplexD, Nt, M> SymAfterPPMatrix; 
extern SymAfterPPMatrix SymAfterPP[Nj][U];

/* signals after MPA, Mod*U*NJ */
typedef Matrix<double, U, NJ> SymAfterMPAMatrix;
extern SymAfterMPAMatrix SymAfterMPA[Mod];

/* codewords after decode */
extern CodeMatrix DecodeBuff[U];

/* final decoded results */
typedef SourceMatrix DecodeMatrix;
extern DecodeMatrix Decode[U];

/*************************************
 * basic global variable declaration
 *************************************/

extern double N_Var;                      /* variance of Noise*/
extern double MaxExponent;                /* Max exponent keep the exp(-alpha/N_var)!=0 */
extern double MaxExponentUpper;           /* Max exponent keep the exp(alpha)!=INF */
extern double BER_TOTAL;                  /* total number of error bits*/
extern double BER;                        /* total number of error bits*/
extern fstream outfile;

/***********************************************
 * basic functions declarations
 ***********************************************/

/**************************************
 * description: normalize the output
 * date: 2020/12/16
 ***************************************/
void NormalIO();

/**************************************
 * description: add AAWGN noise
 * date: 2020/8/20
 * input parameters: variance of noise
 * output parameters: AWGN noise in complex
 ***************************************/
ComplexD AWGN(double nvar);


/**************************************
 * description: overload the operator - between ComplexD & int
 * date: 2020/11/20
 * input parameters: ComplexD comp, int b
 * output parameters: ComplexD tmp
 ***************************************/
ComplexD operator-(ComplexD comp, int b);


/**************************************
 * description: initialize the codebooks & constellations
 * date: 2020/8/13
 * input parameters: MapperType
 * output parameters: MapMatrix& map
 ***************************************/
void InitMapMatrix(MapMatrix& map, 
                   SubConstellMatrix& subConstell,
                   SymAfterMapMatrix* masterConstell);


/**************************************
 * description: calculate the Eb/N0 & initialize the error counter
 * date: 2020/8/18
 * input parameters: ebN0dB
 * output parameters: N_Var
 ***************************************/
void ChannelInitialize(int ebN0dB);


/**************************************
 * description: generater the source bits
 * date: 2020/9/24
 * input parameters: number of the users(U)
 * output parameters: the source bits(Source[U])
 ***************************************/
void BitSource(SourceMatrix& source);


/**************************************
 * description: convolutional code
 * date: 2020/9/24
 ***************************************/
void ConvEncoder(SourceMatrix& source, CodeMatrix& code);


/**************************************
 * description: Modulation
 * date: 2020/10/18
 ***************************************/
void Modulation(CodeMatrix& code, ModuMatrix& modu);


/**************************************
 * description: mapping the codewords based on the codebooks
 * date: 2020/8/12
 ***************************************/
void Mapper(ModuMatrix& modu, 
            SymAfterMapMatrix* symaftermap, 
            MapMatrix& f);


/**************************************
 * description: flat-fading channel
 * date: 2020/8/16
 ***************************************/
void FadingChannel(CSIMatrix* h);


/**************************************
 * description: Beamforming
 * date: 2020/9/20
 ***************************************/
void Beamforming(CSIMatrix* h, 
                 PPMatrix* v, 
                 SymAfterMapMatrix* symaftermap);


/**************************************
 * description: receiving signals
 * date: 2020/8/22
 ***************************************/
void Receiver(SymAfterFCMatrix symAfterFC[Nj][U], 
              SymAfterBFMatrix* symafterbf, 
              CSIMatrix* h);


/**************************************
 * description: decode the receiving signals
 * date: 2020/11/8
 ***************************************/
void MPA(SymAfterPPMatrix symAfterPP[Nj][U], 
         SymAfterMPAMatrix& symAfterMPA, 
         SubConstellMatrix& subConstell);


/**************************************
 * description: decode the receiving signals
 * date: 2020/12/8
 ***************************************/
void ML(SymAfterFCMatrix symAfterFC[Nj][U], 
        SymAfterMPAMatrix* symAfterMPA, 
        SymAfterMapMatrix* masterConstell, 
        CSIMatrix* h, 
        MapMatrix& map);


/**************************************
 * description: soft convolutional decoder
 * date: 2020/11/24
 ***************************************/
void ViterbiSoftDecoder_Prob(SymAfterMPAMatrix* symAfterMPA, 
                             SourceMatrix& source, 
                             SourceMatrix* decode,
                             CodeMatrix& code);

void ViterbiSoftDecoder(SymAfterMPAMatrix* symAfterMPA, 
                             SourceMatrix& source, 
                             SourceMatrix* decode,
                             CodeMatrix& code);

/**************************************
 * description: hard convolutional decoder
 * date: 2020/12/23
 ***************************************/                        
void ViterbiHardDecoder(SymAfterMPAMatrix* symAfterMPA, 
                        SourceMatrix& source, 
                        SourceMatrix* decode, 
                        CodeMatrix& code);


/**************************************
 * description: no convolutional decoder
 * date: 2020/12/23
 ***************************************/                        
void DirectDecoder(SymAfterMPAMatrix* symAfterMPA, 
                   CodeMatrix* decodebuff,
                   CodeMatrix& code);


/**************************************
 * description: Compare
 * date: 2020/12/01
 ***************************************/
void Compare(CodeMatrix* decodebuff, 
             CodeMatrix& code, 
             SourceMatrix& source, 
             SourceMatrix* decode);


/**************************************
 * description: IterationHelper
 * date: 2021/03/01
 ***************************************/
void IterationHelper(CodeMatrix* decodebuff, SourceMatrix* decode, 
                     SymAfterFCMatrix symAfterFC[Nj][U], CSIMatrix* h,
                     SymAfterPPMatrix symAfterPP[Nj][U]);


/**************************************
 * description: IterationHelper-ML
 * date: 2021/03/07
 ***************************************/
void ML_Iter(SymAfterPPMatrix symAfterPP[Nj][U], 
             SymAfterMPAMatrix* symAfterMPA, 
             SymAfterMapMatrix* masterConstell);


/**************************************
 * description: do Svd transform for the CSIMatrix
 * date: 2020/10/21
 ***************************************/
void Svd();

#endif //SCMA_MU_MIMO_MAP_CC_HEADER_H
