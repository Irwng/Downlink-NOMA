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

/**********************************
 * basic constant model parameters
 **********************************/

constexpr int U = 6;                              /* Number of users */
constexpr int M = 2;                              /* Number of time slots resources */
constexpr int Nt = 2;                             /* number of antennas at transmitter */
constexpr int Nr = 2;                             /* number of antennas at recevier */
constexpr int Mod = 2;				  /* BPSK modulation order */
constexpr double PI = 3.141592653589793;
constexpr int Du = 2;                             /* Number of resources connected to user */
constexpr int Dr = 3;                             /* Number of users connected to resources */
constexpr int Mpoint = 64;                        /* Master-constellation Points */
constexpr int Spoint = pow(2, Dr);                         /* Sub-constellation Points */

constexpr double power = 1;

#ifdef DebugMode
    constexpr int MinEbN0dB = 30;
    constexpr long NLoop = pow(10, 0);            /* number of simulation loops  */
    constexpr int MaxEbN0dB = MinEbN0dB;
#else
    constexpr int MinEbN0dB = 0;
    constexpr long NLoop = pow(10, 4);            /* number of simulation loops  */
    constexpr int MaxEbN0dB = 20;           
#endif//DebugMode
constexpr int Step = 2;

constexpr double NUM = NLoop*U;

/*************************************
 * basic global variable declaration
 *************************************/

extern double N_Var;                    /* variance of Noise*/
extern double capacity_Total;           /* total number of error bits without code */
extern double capacity;                 /* BER without code */
extern fstream outfile;

/***********************************************************
 * basic type defination and golbal variables in matrix type
 ***********************************************************/

/* symbols after modulation, Nt*1 */
typedef Matrix<double, Nt, 1> ModuMatrix;

/* channel parameters , U*Nt*Nr */
typedef Matrix<ComplexD, Nr, Nt> CSIMatrix;
extern CSIMatrix H[U];
extern ModuMatrix vectorS;

/***********************************************
 * basic functions declaration
 ***********************************************/

/**************************************
 * description: normalize the output
 * date: 2020/12/16
 ***************************************/
void Initialize();

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
 * description: calculate the Eb/N0 & initialize the error counter
 * date: 2020/8/18
 * input parameters: ebN0dB
 * output parameters: N_Var
 ***************************************/
void ChannelInitialize(int ebN0dB);


/**************************************
 * description: flat-fading channel
 * date: 2020/8/16
 * input parameters: transmitting signals
 * output parameters: receiving signals pointer
 ***************************************/
void FadingChannel(CSIMatrix* h);


/**************************************
 * description: Capacity
 * date: 2021/06/18
 ***************************************/
void Capacity(CSIMatrix* h, ModuMatrix& vectorS);


/**************************************
 * description: do Svd transform for the CSIMatrix
 * date: 2021/03/26
 * input parameters: CSIMatrix h
 * output parameters: CSIMatrix V, U
 ***************************************/
void SVD(CSIMatrix& h, ModuMatrix& vectorS);

#endif //SCMA_ST_PREPROCESS_FADDING_HEADER_H
