/************************************************************
Author:   wangyi     Version : C++11  Date: 2020/10/27
Description: include the functions of Viterbi Soft Decoder in recevier
***********************************************************/
#include "header.h"

/************************************
 * soft decoder based on probability 
 ************************************/

void ViterbiSoftDecoder(SymAfterMPAMatrix* symAfterMPA, SourceMatrix* sourceEsti, CodeMatrix* codeEsti){

#ifndef NoCode

    for(int u = 0; u < U; ++u){
        for(int nj = 0; nj < NJ; ++nj){
            codeEsti[u](nj) = symAfterMPA[1](u,nj) > symAfterMPA[0](u,nj)? 1:0;
        }
    }

    int i = 0;              /* index of bit */
    int step = 0;           /* index of step */
    int k = 0;              /* index of the route */
    int node = 0;           /* index of node */
    int prenode = 0;        /* index of prenode */
    double sum = 0.0;       /* temporary value for normalization */

    for(int u = 0; u < U; ++u){
        
        double Decodebuff[NJ] = {0.0};
        for(i = 0; i < LenBit; i++) sourceEsti[u](i) = 0;
        for(int nj = 0; nj < NJ; nj++){
            Decodebuff[nj] = symAfterMPA[1](u,nj);
        }

    #ifdef DebugMode 
        if(u==0){
            cout<<".................Decodebuff................."<<endl;
            for(int nj = 0; nj < NJ; nj++)    cout<<Decodebuff[nj]<<" ";
            cout<<endl;
        }
    #endif

        /* choose the route */ 
        int route[NODE][LenBit + 1];
        double routeProb[NODE];
        for(auto &dist: routeProb){
            dist = 1;
        }

        /* 
         * initialize the routes of the forward 3 steps
         * " route[A][B] = C " means that the node in the (B+1)th step 
         * the node linked with node A is the Cth node 
         */
        route[0][0] = 0;
        route[1][0] = 0;
        route[2][0] = 0;
        route[3][0] = 0;

        route[0][1] = 0;
        route[1][1] = 2;
        route[2][1] = 0;
        route[3][1] = 2;

        route[0][2] = 0;
        route[1][2] = 1;
        route[2][2] = 2;
        route[3][2] = 3;

        /* 
         * stateChange[i][j][k]: the codeword between the ith node & the jth node, 
         * 000 except (0,0) means no link 
         */
        /* (3,1,3) */
        // double stateChange[4][4][3] = {{{-1,-1,-1},{-1, 1, 1},{ 0, 0, 0},{ 0, 0, 0}},
        //                                {{ 0, 0, 0},{ 0, 0, 0},{-1, 1,-1},{-1,-1, 1}},
        //                                {{ 1, 1, 1},{ 1,-1,-1},{ 0, 0, 0},{ 0, 0, 0}},
        //                                {{ 0, 0, 0},{ 0, 0, 0},{ 1,-1, 1},{ 1, 1,-1}}};
        /* 5,7,7 */
        double stateChange[4][4][3] = {{{ -1,-1,-1},{ 1, 1, 1},{ 0, 0, 0},{ 0, 0, 0}},
                                       {{ 0, 0, 0},{ 0, 0, 0},{-1, 1, 1},{ 1,-1,-1}},
                                       {{ 1, 1, 1},{-1,-1,-1},{ 0, 0, 0},{ 0, 0, 0}},
                                       {{ 0, 0, 0},{ 0, 0, 0},{ 1,-1,-1},{-1, 1, 1}}};

        /* initialize the forward 1th & 2th step */
        for(step = 1; step <= 2; ++step){
            for(i = 0; i < ReciRate; ++i){
                for(node = 0; node < NODE; ++node){
                    routeProb[node] *= stateChange[ route[node][step] ][ route[node][step-1] ][i] == 1 ? 
                                       Decodebuff[i + (step-1) * ReciRate] : 
                                       (1 - Decodebuff[i + (step-1) * ReciRate]);
                }
            }
        }
        
        /* normalization */
        sum = routeProb[0] + routeProb[1] + routeProb[2] + routeProb[3];
        for(auto & tmpprob: routeProb) tmpprob /= sum;
        
    #ifdef DebugMode
        cout<<"..............routeProb............ "<<endl;
        for(auto &tmpprob1:routeProb) cout<<tmpprob1<<" ";
        cout<<endl;
    #endif

        /* index[node][k]: the prenodes linked with this node */
        int index[4][2]={{0, 1},
                         {2, 3},
                         {0, 1},
                         {2, 3}};

        for(step = 3; step <= LenBit; ++step){

            /* find each node's linked prenodes */
            int nodebuff[NODE];              // each node's linked prenode
            double valuebuff[NODE];          // prob between each node's and it's linked prenode

            /* choose each step's route per node */ 
            for(node = 0; node < NODE; ++node){
                double tmp[2] = {0.0};
                // initialize the prob
                for(k = 0; k < Route; ++k){
                    if(route[k][step-1] == index[node][0])
                        tmp[0] = routeProb[k];
                    if(route[k][step-1] == index[node][1])
                        tmp[1] = routeProb[k];
                }

                for(prenode = 0; prenode < 2; ++prenode){//k : forward node
                    for(i = 0; i < ReciRate; ++i){
                        tmp[prenode] *= stateChange[node][index[node][prenode]][i] == 1 ? 
                                        Decodebuff[ (step-1) * ReciRate + i] : 
                                        (1 - Decodebuff[ (step-1) * ReciRate + i] );
                    }
                }
                // cout<<"Soft: "<<"node : "<<node<<" "<<tmp[0]<<" "<<tmp[1]<<endl;
                nodebuff[node] = tmp[0] >= tmp[1] ? index[node][0] : index[node][1];
                valuebuff[node] = tmp[0] >= tmp[1] ? tmp[0] : tmp[1];
            }

            /* cover the useless branches from the tree */
            /* only two situation: 1,1 or 0,2 according to the half nodes */ 
            int prenodebuff[NODE] = {0};//the number of links with prenode
            for(node = 0; node < NODE; node++){
                prenodebuff[ nodebuff[node] ]++;
            }

            if(prenodebuff[0] == 0 || prenodebuff[0] == 2){
                // find the according route of the top2 node
                int tmproute0 = 0, tmproute1 = 0;
                for(k = 0; k < Route; k++){
                    if(route[k][step-1] == 0) tmproute0 = k;
                    if(route[k][step-1] == 1) tmproute1 = k;
                }

                // copy the route to cover the abandoned route
                if(prenodebuff[0] == 0){
                    for(int j = 0; j < step; j++){
                        route[tmproute0][j] = route[tmproute1][j];
                    }
                    routeProb[tmproute0] = routeProb[tmproute1];
                }else{
                    for(int j = 0; j < step; j++){
                        route[tmproute1][j] = route[tmproute0][j];
                    }
                    routeProb[tmproute1] = routeProb[tmproute0];
                }
            }

            if(prenodebuff[2] == 0 || prenodebuff[2] == 2){
                // find the according route of the top2 node
                int tmproute0 = 0, tmproute1 = 0;
                for(k = 0; k < Route; k++){
                    if(route[k][step-1] == 2) tmproute0 = k;
                    if(route[k][step-1] == 3) tmproute1 = k;
                }

                // copy the route to cover the abandoned route
                if(prenodebuff[2] == 0){
                    for(int j = 0; j < step; j++){
                        route[tmproute0][j] = route[tmproute1][j];
                    }
                    routeProb[tmproute0]   = routeProb[tmproute1];
                }else{
                    for(int j = 0; j < step; j++){
                        route[tmproute1][j] = route[tmproute0][j];
                    }
                    routeProb[tmproute1] = routeProb[tmproute0];
                }
            }

            /* here still existing Route(4) routes, begin the normal route connection */
            int routeOccupy[NODE] = {0}; // each route extend once
            for(node = 0; node < NODE; ++node){
                for(k = 0; k < Route; ++k){
                    if(route[k][step-1] == nodebuff[node] && routeOccupy[k] == 0){
                        routeOccupy[k] = 1;
                        route[k][step] = node;
                        routeProb[k] = valuebuff[node];
                        break;
                    }
                }
            }
            
            /* normalization */
            sum = routeProb[0] + routeProb[1] + routeProb[2] + routeProb[3];
            for(auto & tmpprob: routeProb) tmpprob /= sum;
        }

    #ifdef DebugMode 
        cout<<"...................route............ "<<endl;
        for(auto &routetmp: route){
            for(auto &tmp: routetmp) cout<<tmp<<" ";
            cout<<endl;
        }
        cout<<endl;

        cout<<"Prob: "<<endl;
        for(auto &tmp: routeProb) cout<<tmp<<" ";
        cout<<endl;
    #endif

        /* decoder */
        int routePtr = 0; //the choosen route pointer
        double resRoute = routeProb[0];
        for(k = 1;k < Route; ++k){
            if(routeProb[k] > resRoute){
                resRoute = routeProb[k];
                routePtr = k;
            }
        }
        for(step = 1; step <= LenBit; ++step){
            if(route[routePtr][step] >= 2) sourceEsti[u](step-1) = 1;
        }
    #ifdef DebugMode
        cout<<"routePtr: "<<routePtr<<endl;
        cout<<"SourceEsti of "<<u<<": "<<endl<<sourceEsti[u]<<endl;
    #endif
    }

#else
    for(int u = 0; u < U; ++u){
        for(int nj = 0; nj < NJ; ++nj){
            codeEsti[u](nj) = symAfterMPA[1](u,nj) > symAfterMPA[0](u,nj)? 1:0;
        }
    }
    #ifdef DebugMode
        cout<<"codeEsti"<<endl;
        for(int nj = 0; nj < J ; ++nj){
            cout<<codeEsti[0](nj)<<" "; 
        }
        cout<<endl;
    #endif 
#endif
}

/****************************************************************
 * soft decoder based on Euclidean distance after Log-likelihood
 ****************************************************************/

void ViterbiSoftDecoderLL(SymAfterMPAMatrix* symAfterMPA, SourceMatrix* sourceEsti, CodeMatrix* codeEsti){

#ifndef NoCode

    for(int u = 0; u < U; ++u){
        for(int nj = 0; nj < NJ; ++nj){
            codeEsti[u](nj) = symAfterMPA[1](u,nj) > symAfterMPA[0](u,nj)? 1:0;
        }
    }

    int i = 0;              /* index of bit */
    int step = 0;           /* index of step */
    int k = 0;              /* index of the route */
    int node = 0;           /* index of node */
    int prenode = 0;        /* index of prenode */
    double sum = 0.0;       /* temporary value for normalization */

    for(int u = 0; u < U; ++u){
        
        double Decodebuff[NJ][Mod] = {0.0};
        for(i = 0; i < LenBit; i++) sourceEsti[u](i) = 0;
        for(int nj = 0; nj < NJ; ++nj){
            for(int mod = 0; mod < Mod; ++mod){
                Decodebuff[nj][mod] = symAfterMPA[mod](u,nj);
            }
        }

    #ifdef DebugMode 
        if(u==0){
            cout<<".................Decodebuff................."<<endl;
            for(int mod = 0; mod < Mod; ++mod){
                for(int j = 0; j < J; ++j){
                    cout<<Decodebuff[j][mod]<<" ";
                }
                cout<<endl;
            }
        }
    #endif

        /* choose the route */ 
        int route[NODE][LenBit + 1];
        double routeProb[NODE];
        for(auto &dist: routeProb){
            dist = 0;
        }

        /* 
         * initialize the routes of the forward 3 steps
         * " route[A][B] = C " means that the node in the (B+1)th step 
         * the node linked with node A is the Cth node 
         */
        route[0][0] = 0;
        route[1][0] = 0;
        route[2][0] = 0;
        route[3][0] = 0;

        route[0][1] = 0;
        route[1][1] = 2;
        route[2][1] = 0;
        route[3][1] = 2;

        route[0][2] = 0;
        route[1][2] = 1;
        route[2][2] = 2;
        route[3][2] = 3;

        /* 
         * stateChange[i][j][k]: the codeword between the ith node & the jth node, 
         * 000 except (0,0) means no link 
         */
        /* (3,1,3) */
        // double stateChange[4][4][3] = {{{-1,-1,-1},{-1, 1, 1},{ 0, 0, 0},{ 0, 0, 0}},
        //                                {{ 0, 0, 0},{ 0, 0, 0},{-1, 1,-1},{-1,-1, 1}},
        //                                {{ 1, 1, 1},{ 1,-1,-1},{ 0, 0, 0},{ 0, 0, 0}},
        //                                {{ 0, 0, 0},{ 0, 0, 0},{ 1,-1, 1},{ 1, 1,-1}}};
        /* 5,7,7 */
        double stateChange[4][4][3] = {{{-1,-1,-1},{ 1, 1, 1},{ 0, 0, 0},{ 0, 0, 0}},
                                       {{ 0, 0, 0},{ 0, 0, 0},{-1, 1, 1},{ 1,-1,-1}},
                                       {{ 1, 1, 1},{-1,-1,-1},{ 0, 0, 0},{ 0, 0, 0}},
                                       {{ 0, 0, 0},{ 0, 0, 0},{ 1,-1,-1},{-1, 1, 1}}};

        /* initialize the forward 1th & 2th step */
        for(step = 1; step <= 2; ++step){
            for(i = 0; i < ReciRate; ++i){
                for(node = 0; node < NODE; ++node){
                    routeProb[node] += stateChange[ route[node][step] ][ route[node][step-1] ][i] == 1 ? 
                                       Decodebuff[ (step-1) * ReciRate + i][1]: Decodebuff[ (step-1) * ReciRate + i][0];
                }
            }
        }
        
    // #ifdef DebugMode
    //     cout<<"..............routeProb............ "<<endl;
    //     for(auto &tmpprob1:routeProb) cout<<tmpprob1<<" ";
    //     cout<<endl;
    // #endif

        /* index[node][k]: the prenodes linked with this node */
        int index[4][2]={{0, 1},
                         {2, 3},
                         {0, 1},
                         {2, 3}};

        for(step = 3; step <= LenBit; ++step){

            /* find each node's linked prenodes */
            int nodebuff[NODE];              // each node's linked prenode
            double valuebuff[NODE];          // prob between each node's and it's linked prenode

            /* choose each step's route per node */ 
            for(node = 0; node < NODE; ++node){
                double tmp[2] = {0.0};
                // initialize the prob
                for(k = 0; k < Route; ++k){
                    if(route[k][step-1] == index[node][0])
                        tmp[0] = routeProb[k];
                    if(route[k][step-1] == index[node][1])
                        tmp[1] = routeProb[k];
                }

                for(prenode = 0; prenode < 2; ++prenode){//k : forward node
                    for(i = 0; i < ReciRate; ++i){
                        tmp[prenode] += stateChange[node][index[node][prenode]][i] == 1 ? 
                                        Decodebuff[ (step-1) * ReciRate + i][1]: Decodebuff[ (step-1) * ReciRate + i][0];
                    }
                }
                // cout<<"Soft: "<<"node : "<<node<<" "<<tmp[0]<<" "<<tmp[1]<<endl;
                nodebuff[node] = tmp[0] >= tmp[1] ? index[node][0] : index[node][1];
                valuebuff[node] = tmp[0] >= tmp[1] ? tmp[0] : tmp[1];
            }

            /* cover the useless branches from the tree */
            /* only two situation: 1,1 or 0,2 according to the half nodes */ 
            int prenodebuff[NODE] = {0};//the number of links with prenode
            for(node = 0; node < NODE; node++){
                prenodebuff[ nodebuff[node] ]++;
            }

            if(prenodebuff[0] == 0 || prenodebuff[0] == 2){
                // find the according route of the top2 node
                int tmproute0 = 0, tmproute1 = 0;
                for(k = 0; k < Route; k++){
                    if(route[k][step-1] == 0) tmproute0 = k;
                    if(route[k][step-1] == 1) tmproute1 = k;
                }

                // copy the route to cover the abandoned route
                if(prenodebuff[0] == 0){
                    for(int j = 0; j < step; j++){
                        route[tmproute0][j] = route[tmproute1][j];
                    }
                    routeProb[tmproute0] = routeProb[tmproute1];
                }else{
                    for(int j = 0; j < step; j++){
                        route[tmproute1][j] = route[tmproute0][j];
                    }
                    routeProb[tmproute1] = routeProb[tmproute0];
                }
            }

            if(prenodebuff[2] == 0 || prenodebuff[2] == 2){
                // find the according route of the top2 node
                int tmproute0 = 0, tmproute1 = 0;
                for(k = 0; k < Route; k++){
                    if(route[k][step-1] == 2) tmproute0 = k;
                    if(route[k][step-1] == 3) tmproute1 = k;
                }

                // copy the route to cover the abandoned route
                if(prenodebuff[2] == 0){
                    for(int j = 0; j < step; j++){
                        route[tmproute0][j] = route[tmproute1][j];
                    }
                    routeProb[tmproute0]   = routeProb[tmproute1];
                }else{
                    for(int j = 0; j < step; j++){
                        route[tmproute1][j] = route[tmproute0][j];
                    }
                    routeProb[tmproute1] = routeProb[tmproute0];
                }
            }

            /* here still existing Route(4) routes, begin the normal route connection */
            int routeOccupy[NODE] = {0}; // each route extend once
            for(node = 0; node < NODE; ++node){
                for(k = 0; k < Route; ++k){
                    if(route[k][step-1] == nodebuff[node] && routeOccupy[k] == 0){
                        routeOccupy[k] = 1;
                        route[k][step] = node;
                        routeProb[k] = valuebuff[node];
                        break;
                    }
                }
            }
        }

    // #ifdef DebugMode 
    //     cout<<"...................route............ "<<endl;
    //     for(auto &routetmp: route){
    //         for(auto &tmp: routetmp) cout<<tmp<<" ";
    //         cout<<endl;
    //     }
    //     cout<<endl;

    //     cout<<"Prob: "<<endl;
    //     for(auto &tmp: routeProb) cout<<tmp<<" ";
    //     cout<<endl;
    // #endif

        /* decoder */
        int routePtr = 0; //the choosen route pointer
        double resRoute = routeProb[0];
        for(k = 1;k < Route; ++k){
            if(routeProb[k] > resRoute){
                resRoute = routeProb[k];
                routePtr = k;
            }
        }
        for(step = 1; step <= LenBit; ++step){
            if(route[routePtr][step] >= 2) sourceEsti[u](step-1) = 1;
        }
    #ifdef DebugMode
        if(u == 0){
            cout<<"routePtr: "<<routePtr<<endl;
            cout<<"SourceEsti"<<endl;
            cout<<sourceEsti[u]<<endl;
        }
    #endif
    }

#else
    for(int u = 0; u < U; ++u){
        for(int nj = 0; nj < NJ; ++nj){
            codeEsti[u](nj) = symAfterMPA[1](u,nj) > symAfterMPA[0](u,nj)? 1:0;
        }
    }
    #ifdef DebugMode
        cout<<"codeEsti"<<endl;
        for(int nj = 0; nj < J ; ++nj){
            cout<<codeEsti[0](nj)<<" "; 
        }
        cout<<endl;
    #endif 
#endif
}

/*****************************************
 * hard decoder based on hamming distance
 *****************************************/

void ViterbiHardDecoder(SymAfterMPAMatrix* symAfterMPA, SourceMatrix* sourceEsti, CodeMatrix* codeEsti){

#ifndef NoCode
    int i = 0;          /* index of bit */
    int step = 0;       /* index of step */
    int k = 0;          /* index of the route */
    int node = 0;       /* index of node */
    int prenode = 0;    /* index of prenode */

    for(int u = 0; u < U; u++){

        /* Demodulation */
        int Decodebuff[NJ] = {0};
        for(i = 0; i < LenBit; ++i) sourceEsti[u](i) = 0;
        for(int nj = 0; nj < NJ; ++nj){
            Decodebuff[nj] = symAfterMPA[1](u,nj) > symAfterMPA[0](u,nj)? 1:0;
        }

    #ifdef DebugMode 
        if(u==0){
            cout<<".................Decodebuff................."<<endl;
            for(auto tmp:Decodebuff) cout<<tmp<<" ";
            cout<<endl;
        }  
    #endif

        /* choose the route */
        int route[NODE][LenBit + 1];
        int routeHMDistance[NODE];
        for(int &dist:routeHMDistance){
            dist = 0;
        }
        
        /* 
         * initialize the routes of the forward 3 steps
         * " route[a][b] = c " means that the node in the (B-1)th step 
         * linked with node A in the (B-1)th step is the Cth node 
         */
        route[0][0] = 0;
        route[1][0] = 0;
        route[2][0] = 0;
        route[3][0] = 0;

        route[0][1] = 0;
        route[1][1] = 2;
        route[2][1] = 0;
        route[3][1] = 2;

        route[0][2] = 0;
        route[1][2] = 1;
        route[2][2] = 2;
        route[3][2] = 3;

        /* 
         * stateChange[i][j][k]: the codeword between the ith node & the jth node, 
         * 000 except (0,0) means no link 
         */
        int stateChange[4][4][3] = {{{ 0, 0, 0},{ 1, 1, 1},{ 0, 0, 0},{ 0, 0, 0}},
                                    {{ 0, 0, 0},{ 0, 0, 0},{ 0, 1, 1},{ 1, 0, 0}},
                                    {{ 1, 1, 1},{ 0, 0, 0},{ 0, 0, 0},{ 0, 0, 0}},
                                    {{ 0, 0, 0},{ 0, 0, 0},{ 1, 0, 0},{ 0, 1, 1}}};

        for(step = 1; step <= 2; ++step){
            for(i = 0; i < ReciRate; ++i){
                for(node = 0; node < NODE; ++node){
                    routeHMDistance[node] += abs(Decodebuff[i + (step-1)*ReciRate] - \
                                             stateChange[route[node][step]][route[node][step-1]][i]);
                }
            }
        }

        int index[4][2]={{0, 1},
                         {2, 3},
                         {0, 1},
                         {2, 3}};

        for(step = 3; step <= LenBit; ++step){

            /* find each node's linked prenodes */
            int nodebuff[NODE];// each node's linked prenode
            int valuebuff[NODE];// hamming distance between each node's  and it's linked prenode
            
            // choose each step's route per node
            for(node = 0; node < NODE; ++node){
                int tmp[2] = {0};
                // initialize the hamming distance
                for(k = 0; k < Route; ++k){
                    if(route[k][step-1] == index[node][0]){
                        tmp[0] = routeHMDistance[k];
                    }
                    if(route[k][step-1] == index[node][1]){
                        tmp[1] = routeHMDistance[k];
                    }
                }

                for(prenode = 0; prenode < 2; ++prenode){//k : forward node
                    for(i = 0; i < ReciRate; ++i){
                        tmp[prenode] += abs(Decodebuff[ (step-1) * ReciRate + i] - \
                                        stateChange[node][index[node][prenode]][i]);
                    }
                }
                nodebuff[node] = tmp[0] <= tmp[1] ? index[node][0] : index[node][1];
                valuebuff[node] = tmp[0] <= tmp[1] ? tmp[0] : tmp[1];
            }

            /* cover the useless branches from the tree */
            /* only two situation: 1,1 or 0,2 according to the half nodes */ 
            int prenodebuff[NODE] = {0};//the number of links with prenode
            for(node = 0; node < NODE; ++node){
                prenodebuff[nodebuff[node]] ++;
            }

            if(prenodebuff[0] == 0 || prenodebuff[0] == 2){
                // find the according route of the top2 node
                int tmproute0 = 0, tmproute1 = 0;
                for(k = 0; k < Route; ++k){
                    if(route[k][step-1] == 0) tmproute0 = k;
                    if(route[k][step-1] == 1) tmproute1 = k;
                }

                // copy the route to cover the abandoned route
                if(prenodebuff[0] == 0){
                    for(int j = 0; j < step; ++j){
                        route[tmproute0][j] = route[tmproute1][j];
                    }
                    routeHMDistance[tmproute0] = routeHMDistance[tmproute1];
                }else{
                    for(int j = 0; j < step; ++j){
                        route[tmproute1][j] = route[tmproute0][j];
                    }
                    routeHMDistance[tmproute1] = routeHMDistance[tmproute0];
                }
            }
            if(prenodebuff[2] == 0 || prenodebuff[2] == 2){
                // find the according route of the top2 node
                int tmproute0 = 0, tmproute1 = 0;
                for(k = 0; k < Route; ++k){
                    if(route[k][step-1] == 2) tmproute0 = k;
                    if(route[k][step-1] == 3) tmproute1 = k;
                }

                // copy the route to cover the abandoned route
                if(prenodebuff[2] == 0){
                    for(int j = 0; j < step; ++j){
                        route[tmproute0][j] = route[tmproute1][j];
                    }
                    routeHMDistance[tmproute0] = routeHMDistance[tmproute1];
                }else{
                    for(int j = 0; j < step; ++j){
                        route[tmproute1][j] = route[tmproute0][j];
                    }
                    routeHMDistance[tmproute1] = routeHMDistance[tmproute0];
                }
            }

            /* here still existing Route(4) routes, begin the normal route connection */
            int routeOccupy[NODE] = {0}; // each route extend once
            for(node = 0; node < NODE; ++node){
                for(k = 0; k < Route; k++){
                    if(route[k][step-1] == nodebuff[node] && routeOccupy[k] == 0){
                        routeOccupy[k] = 1;
                        route[k][step] = node;
                        routeHMDistance[k] = valuebuff[node];
                        break;
                    }
                }
            }
        }

    #ifdef DebugMode 
        if(u == 0){
            cout<<"...................route..................."<<endl;
            for(auto &routetmp: route){
                for(auto &tmp: routetmp) cout<<tmp<<" ";
                cout<<endl;
            }
            cout<<endl;

            cout<<"...................HMDistance..................."<<endl;
            for(auto &tmp: routeHMDistance) cout<<tmp<<" ";
            cout<<endl;
        }
    #endif

        /* decoder */
        int routePtr = 0; //the choosen route pointer
        int resRoute = routeHMDistance[0];
        for(k = 1; k < Route; ++k){
            if(routeHMDistance[k] < resRoute){
                resRoute = routeHMDistance[k];
                routePtr = k;
            }
        }

        for(step = 1; step <= LenBit; ++step){
            if(route[routePtr][step] >= 2 ) sourceEsti[u](step-1) = 1;
        }

    #ifdef DebugMode
        if(u == 0){
            cout<<"routePtr: "<<routePtr<<endl;
            cout<<"SourceEsti of user"<<u<<": "<<endl<<sourceEsti[u]<<endl;
        }
    #endif
    }
#else
    for(int u = 0; u < U; ++u){
        for(int nj = 0; nj < NJ; ++nj){
            codeEsti[u](nj) = symAfterMPA[1](u,nj) > symAfterMPA[0](u,nj)? 1:0;
        }
    }
#endif
}

