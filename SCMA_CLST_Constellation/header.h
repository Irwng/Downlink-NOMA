/***********************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: include the declaration of constants、variables、functions，
             and the definition of new type
************************************************************************/
#ifndef UDM_DL_CLST_TOVERR_MPA_VITERBI_V3_HEADER_H
#define UDM_DL_CLST_TOVERR_MPA_VITERBI_V3_HEADER_H

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

/* codebook structure mode*/
// #define NOS
// #define NOS2
#define OS
// #define OS2
// #define OS3
// #define TS

/**********************************
 * basic constant model parameters
 **********************************/

constexpr int M = 2;                              /* Number of time slots resources */
constexpr int Nt = 2;                             /* number of antennas at transmitter */
constexpr int J = 8;                              /* Number of user bits each block */
constexpr int Mod = 2;				              /* BPSK modulation order */
constexpr int Mpoint = pow(2, J);                        /* Master-constellation Points */
constexpr double PI = 3.141592653589793;
#ifdef NOS
    constexpr int Dr = 2;
#endif
#ifdef NOS2
    constexpr int Dr = 4;
#endif
#ifdef OS
    constexpr int Dr = 5;
#endif
#ifdef TS
    constexpr int Dr = 4;
#endif
#ifdef OS2
    constexpr int Dr = 5;
#endif
#ifdef OS3
    constexpr int Dr = 3;
#endif

constexpr int Spoint = pow(2, Dr);

/* UDM or PhaseRotation */
constexpr ComplexD p0(0, 0);
extern ComplexD p1;
extern ComplexD p2;
extern ComplexD p3;
extern ComplexD p4;
extern ComplexD p5;

/*************************************
 * basic global variable declaration
 *************************************/

extern fstream outfile1;
extern fstream outfile2;

/***********************************************************
 * basic type defination and golbal variables in matrix type
 ***********************************************************/

/* sub-constelltion point */
typedef Matrix<ComplexD, Nt*M, Spoint> SubConstellMatrix;
extern SubConstellMatrix SubConstell;

/* Mapping matrix*/
typedef Matrix<ComplexD, J, M*Nt> MapMatrix;
extern MapMatrix F;

/* signals after mapping, Nj*Nt*M */
typedef Matrix<ComplexD, Nt, M> SymAfterMapMatrix;
extern SymAfterMapMatrix MasterConstell[Mpoint];


/***********************************************
 * basic functions declaration
 ***********************************************/

/**************************************
 * description: normalize the output
 * date: 2020/12/16
 ***************************************/
void Initialize();


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
 * description: calculate the weight of codewords&minimum Euclidean distance
 * date: 2021/03/06
 * input parameters: subConstell
 * output parameters: MapMatrix& map
 ***************************************/
void Codewords(SubConstellMatrix subConstell, 
               SymAfterMapMatrix* masterConstell);

#endif //UDM_DL_CLST_TOVERR_V3_HEADER_H
