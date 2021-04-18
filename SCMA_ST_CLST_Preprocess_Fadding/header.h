/***********************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: include the declaration of constants、variables、functions，
             and the definition of new type
************************************************************************/
#ifndef SCMA_ST_PREPROCESS_FADDING_HEADER_H
#define SCMA_ST_PREPROCESS_FADDING_HEADER_H

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

/* Comparison */
// #define CompareOfStructure
#define CompareOfBlock
// #define CompareOfTime

/**********************************
 * basic constant model parameters
 **********************************/
#ifdef CompareOfStructure
    constexpr int J = 8;                              /* Length of per user's bit stream,length of vector(b)=K*U */
    constexpr int M = 2;                              /* Number of time slots resources */
    constexpr double BitperSymbol = 2;				            
#endif

#ifdef CompareOfBlock
    constexpr int J = 6;                              /* Length of per user's bit stream,length of vector(b)=K*U */
    constexpr int M = 2;                              /* Number of time slots resources */
    constexpr double BitperSymbol = 1.5;				            
#endif

#ifdef CompareOfTime
    constexpr int J = 8;                              /* Length of per user's bit stream,length of vector(b)=K*U */
    constexpr int M = 2;                              /* Number of time slots resources */
    constexpr double BitperSymbol = 4;				            
#endif

constexpr int U = 6;                              /* Number of users */
constexpr int K = 8;                              /* Length of per user's bit stream,length of vector(b)=K*U */
constexpr int LenBit = U * K;                     /* number of bits of all users */
constexpr double power = 1;
constexpr int NJ = LenBit;                        /* number of coded bits, J*Nj = K*U*ReciRate one block per transmission  */
constexpr int Nj = NJ / J;                        /* number of blocks,J*Nj = K*U one block per transmission  */

constexpr int Nt = 2;                             /* number of antennas at transmitter */
constexpr int Nr = 2;                             /* number of antennas at recevier */
constexpr int Mod = 2;				  /* BPSK modulation order */
constexpr int Mpoint = pow(2, J);                        /* Master-constellation Points */
// constexpr int Spoint = pow(2, J);                         /* Sub-constellation Points */

constexpr double PI = 3.141592653589793;
constexpr int MinEbN0dB = 0;
#ifdef DebugMode
    constexpr long NLoop = pow(10, 0);            /* number of simulation loops  */
    constexpr int MaxEbN0dB = MinEbN0dB;
#else
    constexpr long NLoop = pow(10, 5);            /* number of simulation loops  */
    constexpr int MaxEbN0dB = 30;           
#endif
constexpr int Step = 3;    

/*************************************
 * basic global variable declaration
 *************************************/

/* codewords */
constexpr ComplexD p0(0, 0);
extern ComplexD p1;
extern ComplexD p2;
extern ComplexD p3;
extern ComplexD p4;
extern ComplexD p5;

extern double N_Var;                      /* variance of Noise*/
extern double BER_TOTAL;                  /* total number of error bits*/
extern double BER;                        /* total number of error bits*/
extern fstream outfile;

/***********************************************************
 * basic type defination and golbal variables in matrix type
 ***********************************************************/

/* source codewords */
typedef Matrix<int, 1, LenBit> SourceMatrix;
extern SourceMatrix Source;

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

/* postprocessing matrix, Nj*Nt*Nr */
typedef Matrix<ComplexD, Nt, Nr> PPMatrix;
extern PPMatrix V[U];

/* signals after fadding channal, Nj*U*Nr*M */
typedef Matrix<ComplexD, Nr, M> SymAfterFCMatrix;
extern SymAfterFCMatrix SymAfterFC[Nj][U];

/* signals after post processing, Nj*U*Nt*M */
/* After the post processing, the dimensionality of the signal matrix turns to Nt*M */
typedef Matrix<ComplexD, Nt, M> SymAfterPPMatrix; 
extern SymAfterPPMatrix SymAfterPP[Nj][U];

/* signals after MPA, Nj*U*Nr*M */
typedef Matrix<double, U, NJ> SymAfterMPAMatrix;
extern SymAfterMPAMatrix SymAfterMPA[Mod];

/* final decoded results */
extern SourceMatrix Decode[U];

/***********************************************
 * basic functions declaration
 ***********************************************/

/**************************************
 * description: normalize the output
 * date: 2020/12/16
 ***************************************/
void Initialize(char* argv[]);

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
                   SymAfterMapMatrix* masterConstell,
                   char* argv[]);


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
 * description: Modulation
 * date: 2020/10/18
 * input parameters: CodeMatrix& code
 * output parameters: ModuMatrix& modu
 ***************************************/
void Modulation(SourceMatrix& source, ModuMatrix& modu);


/**************************************
 * description: mapping the codewords based on the codebooks
 * date: 2020/8/12
 * input parameters: codewords array
 * output parameters: map_symbols array
 ***************************************/
void Mapper(ModuMatrix& modu, 
            SymAfterMapMatrix* symAfterMap, 
            MapMatrix& f);


/**************************************
 * description: flat-fading channel
 * date: 2020/8/16
 * input parameters: transmitting signals
 * output parameters: receiving signals pointer
 ***************************************/
void FadingChannel(CSIMatrix* h);


/**************************************
 * description: receiving signals
 * date: 2020/8/22
 * input parameters: SymAfterBFMatrix* symafterbf
 * output parameters: SymAfterPPMatrix symAfterPP[Nj][U]
 ***************************************/
void Receiver(SymAfterFCMatrix symAfterFC[Nj][U], 
              SymAfterMapMatrix* symAfterMap, 
              CSIMatrix* h,
              SymAfterPPMatrix symAfterPP[Nj][U], 
              PPMatrix* v, 
              char* argv[]);


/**************************************
 * description: receiving signals
 * date: 2021/04/14
 * input parameters: SymAfterMapMatrix* symAfterMap
 * output parameters: SymAfterPPMatrix symAfterPP[Nj][U]
 ***************************************/
void Receiver_OSIC(SymAfterFCMatrix symAfterFC[Nj][U], 
              SymAfterMapMatrix* symAfterMap, 
              CSIMatrix* h,
              SymAfterPPMatrix symAfterPP[Nj][U], 
              PPMatrix* v, 
              char* argv[]);


/**************************************
 * description: decode the receiving signals
 * date: 2020/12/8
 * input parameters: SymAfterPPMatrix symAfterPP[Nj][U]
 * output parameters: SymAfterMPAMatrix& symAfterMPA
 ***************************************/
void ML(SymAfterPPMatrix symAfterPP[Nj][U], 
        SymAfterMPAMatrix* symAfterMPA, 
        SymAfterMapMatrix* masterConstell);


/**************************************
 * description: no convolutional decoder
 * date: 2020/12/23
 * input parameters: SymAfterPPMatrix symAfterMPA
 * output parameters: DecodeMatrix decode
 ***************************************/                        
void DirectDecoder(SymAfterMPAMatrix* symAfterMPA, 
                   SourceMatrix* decode, 
                   SourceMatrix& source);


/**************************************
 * description: Compare
 * date: 2020/12/01
 * input parameters: CodeMatrix* decodebuff, 
 * output parameters: CodeMatrix& code
 ***************************************/
void Compare(SourceMatrix* decode, 
             SourceMatrix& source);

#endif //SCMA_ST_PREPROCESS_FADDING_HEADER_H
