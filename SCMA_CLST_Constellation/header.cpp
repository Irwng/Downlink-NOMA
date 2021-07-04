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
ComplexD p4;
ComplexD p5;

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
                   SymAfterMapMatrix* masterConstell){

    double weight;

    /* initialize the mapping matrix's codewords */
#ifdef NOS
    /* J = 8, NOS */
    cout<<"Non-overlapping section(NOS)"<<endl;
    outfile1<<"Non-overlapping section(NOS)"<<endl;
    outfile2<<"Non-overlapping section(NOS)"<<endl;
    
    weight = sqrt(2);
    p1 = ComplexD(1/weight, 0/weight);
    p2 = ComplexD(0/weight, 1/weight);

    /*   Nt= 0 Nt= 1  */
    map<<p1,p0,p0,p0,
            p2,p0,p0,p0,
            p0,p1,p0,p0,
            p0,p2,p0,p0,
            p0,p0,p1,p0,
            p0,p0,p2,p0,
            p0,p0,p0,p1,
            p0,p0,p0,p2;
#endif

#ifdef NOS2
    /* J = 8, NOS */
    cout<<"Non-overlapping section(NOS):M=1"<<endl;
    outfile1<<"Non-overlapping section(NOS):M=1"<<endl;
    outfile2<<"Non-overlapping section(NOS):M=1"<<endl;
    
    weight = sqrt(10);
    p1 = ComplexD(1/weight, 0/weight);
    p2 = ComplexD(0/weight, 1/weight);
    p3 = ComplexD(0/weight, 2/weight);
    p4 = ComplexD(2/weight, 0/weight);

    /*   Nt= 0 Nt= 1  */
    map<<p1,p0,p0,p0,
         p2,p0,p0,p0,
         p3,p0,p0,p0,
         p4,p0,p0,p0,
         p0,p0,p1,p0,
         p0,p0,p2,p0,
         p0,p0,p3,p0,
         p0,p0,p4,p0;
#endif


#ifdef OS
    /* J = 8, OS */
    cout<<"Overlapping section(OS)"<<endl;
    outfile1<<"Overlapping section(OS)"<<endl;
    outfile2<<"Overlapping section(OS)"<<endl;
    
    weight = sqrt(26);
    p1 = ComplexD(1/weight, 0/weight);
    p2 = ComplexD(0/weight, 1/weight);
    p3 = ComplexD(0/weight, 2/weight);
    p4 = ComplexD(2/weight, 0/weight);
    p5 = ComplexD(4/weight, 0/weight);

    /*   Nt= 0 Nt= 1  */
    map<<p1,p0,p0,p0,
         p2,p1,p0,p0,
         p3,p2,p1,p0,
         p4,p3,p2,p1,
         p5,p4,p3,p2,
         p0,p5,p4,p3,
         p0,p0,p5,p4,
         p0,p0,p0,p5;
#endif

#ifdef TS      
    /* J = 8, TS */
    cout<<"Tailbiting section(TS)"<<endl;
    outfile1<<"Tailbiting section(TS)"<<endl;
    outfile2<<"Tailbiting section(TS)"<<endl;
    
    weight = sqrt(10);
    p1 = ComplexD(1/weight, 0/weight);
    p2 = ComplexD(0/weight, 1/weight);
    p3 = ComplexD(0/weight, 2/weight);
    p4 = ComplexD(2/weight, 0/weight);

    /*   Nt= 0 Nt= 1  */
    map<<p1,p0,p0,p3,
            p2,p0,p0,p4,
            p3,p1,p0,p0,
            p4,p2,p0,p0,
            p0,p3,p1,p0,
            p0,p4,p2,p0,
            p0,p0,p3,p1,
            p0,p0,p4,p2;
#endif
    
#ifdef OS2
    /* J = 8, OS */
    cout<<"Overlapping section(OS)"<<endl;
    outfile1<<"Overlapping section(OS)"<<endl;
    outfile2<<"Overlapping section(OS)"<<endl;
    
    weight = sqrt(26);
    p1 = ComplexD(1/weight, 0/weight);
    p2 = ComplexD(0/weight, 1/weight);
    p3 = ComplexD(2/weight, 0/weight);
    p4 = ComplexD(0/weight, 2/weight);
    p5 = ComplexD(4/weight, 0/weight);

    /*   Nt= 0 Nt= 1  */
    map<<p1,p0,p0,p0,
         p2,p0,p0,p0,
         p3,p0,p0,p0,
         p4,p0,p1,p0,
         p5,p0,p2,p0,
         p0,p0,p3,p0,
         p0,p0,p4,p0,
         p0,p0,p5,p0;
#endif

#ifdef OS3
    /* J = 8, OS */
    cout<<"Overlapping section(OS)"<<endl;
    outfile1<<"Overlapping section(OS)"<<endl;
    outfile2<<"Overlapping section(OS)"<<endl;
    
    weight = sqrt(6);
    p1 = ComplexD(1/weight, 0/weight);
    p2 = ComplexD(0/weight, 1/weight);
    p3 = ComplexD(0/weight, 2/weight);

    /*   Nt= 0 Nt= 1  */
    map<<p1,p0,p0,p0,
            p2,p1,p0,p0,
            p3,p2,p1,p0,
            p0,p3,p2,p1,
            p0,p0,p3,p2,
            p0,p0,p0,p3,
            p0,p0,p0,p0,
            p0,p0,p0,p0;
#endif


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
                    masterConstell[i](nt, m) += map(j, nt * M + m) * static_cast<double>(2 * bit[j] - 1);
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
