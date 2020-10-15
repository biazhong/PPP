//
//  _batchingcpp.cpp
//  
//
// 
//

#include <iostream>
#include <vector>
#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fstream>


#include <chrono>
#include <thread>
using namespace std;
using namespace std::this_thread;
using namespace std::chrono; 

double exprand(double lambda){
    double u;
    
    u = (double)rand() / (RAND_MAX * 1.0);
    
    return -log(1- u) / lambda;
}



typedef struct{
    int label;
    int sampleSize;
    float sampleXSum;
    float sampleX2Sum;
} alt;



void master_send_out_tasks(vector<alt>* outgoing_alts, int outgoing_rank){
    MPI_Send((void*)outgoing_alts->data(),outgoing_alts->size()*sizeof(alt), MPI_BYTE, outgoing_rank, 0, MPI_COMM_WORLD);
}



void master_receive_tasks(vector<alt>* incoming_alts){
    MPI_Status status;
    int incoming_rank;
    int incoming_size;
    
    MPI_Probe(MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
    MPI_Get_count(&status, MPI_BYTE,&incoming_size);
    incoming_rank = status.MPI_SOURCE;
    
    incoming_alts->resize(incoming_size / sizeof(alt));
    
    MPI_Recv((void*)incoming_alts->data(),incoming_size, MPI_BYTE,incoming_rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}


void worker_send_out_tasks(vector<alt>* outgoing_alts){
    MPI_Send((void*)outgoing_alts->data(), outgoing_alts->size()*sizeof(alt), MPI_BYTE, 0, 0 ,MPI_COMM_WORLD);
}


void worker_receive_tasks(vector<alt>* incoming_alts){
    MPI_Status status;
    int incoming_rank = 0;
    int incoming_size;
    
    MPI_Probe(incoming_rank,0,MPI_COMM_WORLD,&status);
    MPI_Get_count(&status,MPI_BYTE,&incoming_size);
    
    incoming_alts->resize(incoming_size / sizeof(alt));
    
    MPI_Recv((void*)incoming_alts->data(),incoming_size, MPI_BYTE,incoming_rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}



void generate_obv (double* _sim, int label){
    int RB = 50;
    vector<int> x_disc(5);
    int rr = RB * 2 - 3;
    int n = label/(RB-1)+1;
    double _t = MPI_Wtime();


    x_disc[3] =  label%(RB-1) + 1;
    x_disc[4] = RB - x_disc[3];
    x_disc[0] = (int)ceil(((float)rr - sqrt((float)rr*rr - 8. * n )) / 2 );
    x_disc[1] = RB - 1 + (n - (rr-x_disc[0]) * x_disc[0] / 2) - x_disc[0];
    x_disc[2] = RB - x_disc[0] - x_disc[1];
    
    int i, j; double t;
    
    
    int _njobs = 2050;// (int)(floor(rand()/(RAND_MAX * 1.0)*1000)+1550);
    int _nstages = 3;
    int _burnin = _njobs-50;
    int _diff = _njobs - _burnin;
    
    vector<int> r(3);
    r[0] = x_disc[0];
    r[1] = x_disc[1];
    r[2] = x_disc[2];
    vector<int> b(2);
    b[0] = x_disc[3];
    b[1] = x_disc[4];
    
    vector<double> st(_nstages);
    vector<vector<double> > sTime (_njobs,st);
    
    vector<double> et(_njobs+1);
    vector<vector<double> > eTime(_nstages,et);
    
    for(i = 0 ; i < _njobs; i++){
        for(j = 0; j < _nstages; j++){
            sTime[i][j] = exprand(r[j]);
        }
    }
    
    for(i = 1; i <= _njobs; ++i){
        t = sTime[i-1][0];
        if (eTime[1][max(0,i-b[0])] <= eTime[0][i-1]+t){
            eTime[0][i] = eTime[0][i-1]+t;
        }else {
            eTime[0][i] = eTime[1][max(0,i-b[0])];
        }
        
        for (j=2; j<=_nstages-1; j++) {
            t = sTime[i-1][j-1];
            if(eTime[j-1][i-1]>eTime[j-2][i]) {
                if(eTime[j][max(0,i-b[j-1])] <= eTime[j-1][i-1]+t) {
                    eTime[j-1][i] = eTime[j-1][i-1]+t;
                } else {
                    eTime[j-1][i] = eTime[j][max(0,i-b[j-1])];
                }
            } else {
                if (eTime[j][max(0,i-b[j-1])] <= eTime[j-2][i]+t) {
                    eTime[j-1][i] = eTime[j-2][i]+t;
                } else{
                    eTime[j-1][i] = eTime[j][max(0,i-b[j-1])];
                }
            }
        }
        
        t=sTime[i-1][_nstages-1];
        
        if (eTime[_nstages-1][i-1] <= eTime[_nstages-2][i]) {
            eTime[_nstages-1][i] = eTime[_nstages-2][i]+t;
        } else {
            eTime[_nstages-1][i] = eTime[_nstages-1][i-1]+t;
        }
    }

    _t = MPI_Wtime()-_t;
    //int sleepTime = (int)floor(_t*(rand()/(RAND_MAX * 1.0)*1.0+0.5)*1000000000);
    int sleepTime = (int)floor((rand()/(RAND_MAX * 1.0)*1.0+0.5)*1000000);
    //sleep_for(nanoseconds(sleepTime));
    double _t_ = MPI_Wtime();
    double stop = 0.0;
    while(stop<sleepTime){
    	stop = (MPI_Wtime()-_t_)*1000000000;		
    }
    
    *_sim = ((double)(_njobs-_burnin))/(eTime[_nstages-1][_njobs]-eTime[_nstages-1][_burnin]);
}


int main(int argc, char** argv){
    int k = 57624;
    double _starttime = 0.0, _endtime;
    MPI_Init(NULL,NULL);
    int repeat = 100;
    vector<double> recordTime;    

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    srand(time(NULL)*world_rank);
    vector<alt> outgoing_alts;
    vector<alt> incoming_alts;
    double obv;

    for(int count = 0; count < repeat; count++){
    
    if(world_rank == 0){
        _starttime = MPI_Wtime();
        vector<alt> I;
        alt tempAlt;
        for (int i = 0; i <k; i++){
            tempAlt.label = i+1;
            tempAlt.sampleSize = 0;
            tempAlt.sampleXSum = 0.0;
            tempAlt.sampleX2Sum = 0.0;
            I.push_back(tempAlt);
        }
        for(int i = 1; i <= world_size - 1; i++){
            for(int j = 0; j < k; j++){
                if((j%(world_size-1)+1)==i){
                    outgoing_alts.push_back(I[j]);
                }
            }
            
            master_send_out_tasks(&outgoing_alts, i);
            outgoing_alts.clear();
        }
        
        for(int i = 0; i < world_size-1; i++){
            master_receive_tasks(&incoming_alts);
            for(int j = 0; j < incoming_alts.size(); j++){
                I[incoming_alts[j].label-1].sampleSize = incoming_alts[j].sampleSize;
                I[incoming_alts[j].label-1].sampleXSum = incoming_alts[j].sampleXSum;
                I[incoming_alts[j].label-1].sampleX2Sum = incoming_alts[j].sampleX2Sum;
            }
            incoming_alts.clear();
        }
        _endtime = MPI_Wtime();
        recordTime.push_back(-_starttime+_endtime);
        printf("Batching: The wall clock time is %f\n", -_starttime+_endtime);
        
    }else{
        worker_receive_tasks(&incoming_alts);
        for(int i = 0; i < incoming_alts.size(); i++){
            
            generate_obv(&obv, incoming_alts[i].label);
	    
            incoming_alts[i].sampleSize++;
            incoming_alts[i].sampleXSum = incoming_alts[i].sampleXSum + obv;
            incoming_alts[i].sampleX2Sum = incoming_alts[i].sampleX2Sum + obv * obv;
        }
        worker_send_out_tasks(&incoming_alts);
        incoming_alts.clear();
    }

    }
    if(world_rank==0){
	ofstream myfile ("Batchingoutput.txt");
	for(int i = 0; i < recordTime.size();i++){
        myfile<<recordTime[i]<<"\n";
        }
        //    myfile<<"The wall clock time is: "<<-_starttime+_endtime<<"\n"<<"The best alternative is: "<<bestLabel<<"\n"<<"The total sample size is: "<<totalSampleSize<<"\n"<<"Total comparison time is: "<<total_comparison_time<<"\n";
        myfile.close();
    }
    MPI_Finalize();
    return 0;
    
}

