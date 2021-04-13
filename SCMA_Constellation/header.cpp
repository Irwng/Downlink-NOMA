/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: include the basic functions of transmitter & recevier & channels
***********************************************************/
#include "header.h"

fstream outfile1;
fstream outfile2;

ComplexD p1;
ComplexD p2;
ComplexD p3;

SubConstellMatrix SubConstell;                      /* sub-constelltion point */
SymAfterMapMatrix MasterConstell[Mpoint];           /* Master-constelltion point */
MapMatrix F;                                        /* Mapping matrix */


void Initialize(){

    outfile1.open("Euclidean distance of subConstellation.txt", ios::out|ios::app);
    outfile1<<"....................................................."<<endl;

    outfile2.open("Euclidean distance of masterConstellation.txt", ios::out|ios::app);
    outfile2<<"....................................................."<<endl;
}


void InitMapMatrix(MapMatrix& map, SubConstellMatrix& subConstell, 
                   SymAfterMapMatrix* masterConstell, 
                   char* argv[]){

    double weight;
    double RF;

    /* initialize the mapping matrix's codewords */
    switch(*argv[1]){
        /* UDM */
        case '1':
            /* 202,2,4,9 */
            cout<<"UDM-202,2,4,9"<<endl;
            outfile1<<"UDM-202,2,4,9"<<endl;
            outfile2<<"UDM-202,2,4,9"<<endl;

            weight = sqrt(202);
            p1 = ComplexD(2/weight, 9/weight);
            p2 = ComplexD(4/weight, 4/weight);
            p3 = ComplexD(9/weight, 2/weight);
            break;

        case '2':
            /* 42, 1, 2, 4 */
            cout<<"UDM-42,1,2,4"<<endl;
            outfile1<<"UDM-42,1,2,4"<<endl;
            outfile2<<"UDM-42,1,2,4"<<endl;

            weight = sqrt(42);
            p1 = ComplexD(1/weight, 4/weight);
            p2 = ComplexD(2/weight, 2/weight);
            p3 = ComplexD(4/weight, 1/weight);
            break;

        case '3':
            /* 178, 3, 4, 8 */
            cout<<"UDM-178,3,4,8"<<endl;
            outfile1<<"UDM-178,3,4,8"<<endl;
            outfile2<<"UDM-178,3,4,8"<<endl;

            weight = sqrt(178);
            p1 = ComplexD(3/weight, 8/weight);
            p2 = ComplexD(4/weight, 4/weight);
            p3 = ComplexD(8/weight, 3/weight);
            break;

        case '4':
            /* 178, 3, 4, 8 */
            cout<<"UDM-1,2i,1+1i"<<endl;
            outfile1<<"UDM-1,2i,1+1i"<<endl;
            outfile2<<"UDM-1,2i,1+1i"<<endl;

            weight = sqrt(7);
            p1 = ComplexD(1/weight, 0/weight);
            p2 = ComplexD(1/weight, 1/weight);
            p3 = ComplexD(0/weight, 2/weight);
            break;


        /* PhaseRotation */
        case '6':
            /* RF=pi/4 */
            cout<<"PhaseRotation-pi/4"<<endl;
            outfile1<<"PhaseRotation-pi/4"<<endl;
            outfile2<<"PhaseRotation-pi/4"<<endl;

            weight = sqrt(3);
            RF = PI/4;
            p1 = ComplexD(cos(0*RF)/weight, sin(0*RF)/weight);
            p2 = ComplexD(cos(1*RF)/weight, sin(1*RF)/weight);
            p3 = ComplexD(cos(2*RF)/weight, sin(2*RF)/weight);    
            break;

        case '7':
            /* RF=pi/5 */
            cout<<"PhaseRotation-pi/5"<<endl;
            outfile1<<"PhaseRotation-pi/5"<<endl;
            outfile2<<"PhaseRotation-pi/5"<<endl;

            weight = sqrt(3);
            RF = PI/5;
            p1 = ComplexD(cos(0*RF)/weight, sin(0*RF)/weight);
            p2 = ComplexD(cos(1*RF)/weight, sin(1*RF)/weight);
            p3 = ComplexD(cos(2*RF)/weight, sin(2*RF)/weight);    
            break;

        case '8':
            /* RF=pi/6 */
            cout<<"PhaseRotation-pi/6"<<endl;
            outfile1<<"PhaseRotation-pi/6"<<endl;
            outfile2<<"PhaseRotation-pi/6"<<endl;

            weight = sqrt(3);
            RF = PI/6;
            p1 = ComplexD(cos(0*RF)/weight, sin(0*RF)/weight);
            p2 = ComplexD(cos(1*RF)/weight, sin(1*RF)/weight);
            p3 = ComplexD(cos(2*RF)/weight, sin(2*RF)/weight);    
            break;

        default:
            break;
    }
    
    outfile1<<"p1: "<<p1<<endl<<"p2: "<<p2<<endl<<"p3: "<<p3<<endl;
    outfile2<<"p1: "<<p1<<endl<<"p2: "<<p2<<endl<<"p3: "<<p3<<endl;
    
    /* initialize the mapping matrix's structure */
    switch(*argv[2]){
        
        /* Latin */
        case '1':
            cout<<"Latin"<<endl;
            outfile1<<"Latin"<<endl;
            outfile2<<"Latin"<<endl;

              /* Nt= 0 Nt= 1 */
            map<<p1,p2,p0,p0,
                 p2,p0,p3,p0,
                 p3,p0,p0,p1,
                 p0,p3,p1,p0,
                 p0,p3,p0,p2,
                 p0,p0,p2,p3;
            break;

        /* CLST */
        case '2':
            cout<<"CLST"<<endl;
            outfile1<<"CLST"<<endl;
            outfile2<<"CLST"<<endl;

              /* Nt= 0 Nt= 1 */
            map<<p1,p0,p0,p0,
                 p2,p1,p0,p0,
                 p3,p2,p1,p0,
                 p0,p3,p2,p1,
                 p0,p0,p3,p2,
                 p0,p0,p0,p3;
            break;

        /* Orthogonal */
        case '3':
            cout<<"Orthogonal"<<endl;
            outfile1<<"Orthogonal"<<endl;
            outfile2<<"Orthogonal"<<endl;

            /* Nt= 0 Nt= 1 */
            map<<p1,p0,p0,p1,
                 p2,p0,p0,p2,
                 p3,p0,p0,p3,
                 p0,p1,p1,p0,
                 p0,p2,p2,p0,
                 p0,p3,p3,p0;
            break;

        /* half-Latin */
        case '4':
            cout<<"half-Latin"<<endl;
            outfile1<<"half-Latin"<<endl;
            outfile2<<"half-Latin"<<endl;

            /* Nt= 0 Nt= 1 */
            map<<p1,p1,p0,p0,
                 p2,p0,p1,p0,
                 p3,p0,p0,p1,
                 p0,p2,p2,p0,
                 p0,p3,p0,p2,
                 p0,p0,p3,p3;
            break;
        default:
            break;
    }

    /* initialize the subconstellation */
    int a,b;
    double subMsg[Spoint][Dr] = {0.0};

    for(int i = 0; i < Spoint; i++){
        for(int j = Dr; j > 0; j--){
            a = (int)pow(2,-j+Dr+1);
            b = (int)pow(2,-j+Dr);
            if(i%a >= b)
                subMsg[i][j-1] = 1;                             /* 1 - -1,0 - 1 */
            else
                subMsg[i][j-1] = -1;
        }
    } 

    int flag = 0;
    for(int k = 0; k < Nt*M; k++){
		for(int i = 0; i < Spoint; i++){
			flag = 0;
			for(int index = 0; index < J; index++){
				if(map(index,k) != p0){
					subConstell(k,i) += map(index, k) * subMsg[i][flag];
					flag++;
				}
			}
		}
	}

    /* initialize the Masterconstellation */
    for(int i = 0; i < Mpoint; ++i){                                  
        bitset<J> bit(i);
        for(int nt = 0; nt < Nt; ++nt){
            for(int m = 0; m < M; ++m){
                masterConstell[i](nt, m) = p0;
                for(int j = 0; j < J; ++j){
                    masterConstell[i](nt, m) += map(j, nt * Nt + m) * static_cast<double>(2 * bit[j] - 1);
                }
            }
        }
    }  

    Codewords(subConstell, masterConstell);
}


void Codewords(SubConstellMatrix subConstell, SymAfterMapMatrix* masterConstell){

    /* show the subconstellation */
    cout<<"subConstell"<<endl;
    outfile1<<"subConstell"<<endl;

    for(int i  = 0; i < Spoint; ++i){
        cout<<subConstell(0, i)<<endl;
        outfile1<<subConstell(0, i).real()<<" "<<subConstell(0, i).imag()<<endl;
    }

    /* calculate the weight of codewords */
    double tmpweight = 0.0;
    for(int i = 0; i < Spoint; i++){
        tmpweight += pow(abs(subConstell(0,i)), 2);
    }
    cout<<"subconstell weight: "<<tmpweight/static_cast<double>(Spoint)<<endl;
    outfile1<<"subconstell weight: "<<tmpweight/static_cast<double>(Spoint)<<endl;
    
    /* calculate the minimum Euclidean distance */
    double distance[Spoint][Spoint];
    for(auto &a : distance){
        for(auto &b: a){
            b = 0.0;
        }
    }

    double mindist = INT16_MAX;
    int x = 0, y = 0;
    for(int i = 0; i < Spoint; ++i){
        for(int j = i+1; j < Spoint; ++j){
            distance[i][j] = abs(subConstell(0,i) - subConstell(0,j)); 
            if(distance[i][j] < mindist){
                mindist = distance[i][j];
                x = i;
                y = j;
            }
        }
    }

    // cout<<"distance"<<endl;
    // for(auto &a : distance){
    //     for(auto &b: a){
    //         outfile1<<b<<" ";
    //     }
    //     outfile1<<endl;
    // }
    cout<<"The minimum distance is "<<mindist<<" between "<<subConstell(0, x)<<" and "<<subConstell(0,y)<<endl;
    outfile1<<"The minimum distance is "<<mindist<<" between "<<subConstell(0, x)<<" and "<<subConstell(0,y)<<"."<<endl;
    outfile1.close();

    /* show the masterconstellation */
    // outfile2<<"MasterConstell: "<<endl;
    double tmpweight2 = 0.0;
    for(int mpoint = 0; mpoint < Mpoint; ++mpoint){
        for(int nt = 0; nt < Nt; ++nt){
            for(int m = 0; m < M; ++m){
                tmpweight2 += pow(abs(masterConstell[mpoint](nt, m)),2); 
            }
        }
    }

    cout<<"MasterConstell weight: "<<tmpweight2/static_cast<double>(Mpoint)<<endl;
    outfile2<<"MasterConstell weight: "<<tmpweight2/static_cast<double>(Mpoint)<<endl;
    
    /* calculate the minimum Euclidean distance */
    double distance2[Mpoint][Mpoint];
    for(auto &a : distance2){
        for(auto &b: a){
            b = 0.0;
        }
    }
    double mindist2 = INT16_MAX;
    int x2 = 0, y2 = 0;
    long count = 0;
    for(int i = 0; i < Mpoint; ++i){
        for(int j = i+1; j < Mpoint; ++j){
            for(int nt = 0; nt < Nt; ++nt){
                for(int m = 0; m < M; ++m){
                    ++count;            
                    distance2[i][j] += abs(masterConstell[i](nt, m) - masterConstell[j](nt, m)); 
                }
            }
            if(distance2[i][j] < mindist2){
                mindist2 = distance2[i][j];
                x2 = i;
                y2 = j;
            }
        }
    }

    // for(auto &a : distance2){
    //     for(auto &b: a){
    //         outfile2<<b<<" ";
    //     }
    //     outfile2<<endl;
    // }
    cout<<"The minimum distance is "<<mindist2<<" between "<<endl;
    cout<<masterConstell[x2]<<endl<<" and "<<masterConstell[y2]<<"."<<endl;

    outfile2<<"The minimum distance is "<<mindist2<<" between "<<endl;
    outfile2<<masterConstell[x2]<<endl<<" and "<<endl<<masterConstell[y2]<<"."<<endl;
    outfile2.close();

}
