//
//  GSP.cpp
// This code implements the GSP of Ni et al.（2017） in MPI.
// The simulation optimization problem considered in this code is the three-stage buffer allocation problem.
// The parameters follows by comments "Input Parameter:..." should be adjusted from one problem instance to another.
// The parameter eta could be calculated by EtaFunc.java and Rinott's constant can be calculated by Rinott.java
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
    
    u = (double)rand() / (RAND_MAX * 1.0+1);
    
    return -log(1- u) / lambda;
}


typedef struct{
    int label;
    int batchsize;
    int sampleSize;
    double sampleXSum;
    double S2;
    
    
} alt;


void master_send_out_tasks(vector<alt>* outgoing_alts, int outgoing_rank){
    MPI_Send((void*)outgoing_alts->data(),outgoing_alts->size()*sizeof(alt), MPI_BYTE, outgoing_rank, 0, MPI_COMM_WORLD);
}


void master_receive_tasks(vector<alt>* incoming_alts, int* incoming_rank){
    MPI_Status status;
    int incoming_size;
    
    MPI_Probe(MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
    MPI_Get_count(&status, MPI_BYTE,&incoming_size);
    *incoming_rank = status.MPI_SOURCE;
    
    incoming_alts->resize(incoming_size / sizeof(alt));
    
    MPI_Recv((void*)incoming_alts->data(),incoming_size, MPI_BYTE,*incoming_rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
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
    int RB = 128;   //Input Parameter: Problem parameter
    vector<int> x_disc(5);
    int rr = RB * 2 - 3;
    int n = label/(RB-1)+1;
    //double _t = MPI_Wtime();    
    
    x_disc[3] =  label%(RB-1) + 1;
    x_disc[4] = RB - x_disc[3];
    x_disc[0] = (int)ceil(((double)rr - sqrt((double)rr*rr - 8. * n )) / 2 );
    x_disc[1] = RB - 1 + (n - (rr-x_disc[0]) * x_disc[0] / 2) - x_disc[0];
    x_disc[2] = RB - x_disc[0] - x_disc[1];
    
    int i, j; double t;
    
    
    
    int _njobs = 2050;
    int _nstages = 3;
    int _burnin = 2000;
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
	
    *_sim = ((double)(_njobs-_burnin))/(eTime[_nstages-1][_njobs]-eTime[_nstages-1][_burnin]);
}


int main(int argc, char** argv){

    MPI_Init(NULL,NULL);
    int repeat = 10;
    
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    for(int count1 = 0; count1 < repeat; count1++){
    int k = 1016127; //Input Parameter: Total number of alternatives
    int n0 = 50; //Input Parameter: First stage sample size
    double alpha = 0.05; //Input Parameter: Desired PAC
    double delta = 0.1; //Input Parameter: The IZ parameter delta
    int beta = 200;     //Input Parameter: Parameter beta for GSP
    int maxR = 5;       //Input Parameter: Parameter rbar for GSP
    int bestsys;
    int r = 0;
    double eta = 0.9744071960449219;    //Input Parameter: Parameter eta used for constructing the continuation region. NOTE: For different problems one should use EtaFunc.java to calculate eta values.
    double rinotth = 8.527191162109375;     //Input Parameter: Rinott's constant used for determine the maximum sample sizes of different alternativfes. NOTE: For different problems one should use Rinott.java to calculate Rinott's constants.
    
    

    int incoming_rank;
    double total_comparison_time = 0.0;
    double total_sim_time = 0.0;
    double _tt;
    
    double totalSampleSize = 0;
    
    double _starttime = 0.0, _endtime;
    
    srand(time(NULL)*world_rank);
    vector<alt> outgoing_alts;
    vector<alt> incoming_alts;
    
    
    
    double obv=0.0;
    vector<alt> I;
    vector<int> _record_surviving;
    int num_of_surviving;
    
    if(world_rank == 0){
        _starttime = MPI_Wtime();
        alt tempAlt;
        for (int i = 0; i < k; i++){
            tempAlt.label = i;
            tempAlt.batchsize = beta;
            tempAlt.sampleSize = 0;
            tempAlt.sampleXSum = 0.0;
            tempAlt.S2 = 0.0;
            I.push_back(tempAlt);
            _record_surviving.push_back(i);
        }
        for(int i = 1; i <= world_size - 1; i++){
            for(int j = 0; j < _record_surviving.size(); j++){
                if((j%(world_size-1)+1)==i){
                    outgoing_alts.push_back(I[_record_surviving[j]]);
                }
            }
            
            master_send_out_tasks(&outgoing_alts, i);
            outgoing_alts.clear();
        }
        for(int i = 0; i < world_size-1; i++){
            master_receive_tasks(&incoming_alts,&incoming_rank);
            for(int j = 0; j < incoming_alts.size(); j++){
                I[incoming_alts[j].label].sampleSize = incoming_alts[j].sampleSize;
                I[incoming_alts[j].label].sampleXSum = incoming_alts[j].sampleXSum;
                I[incoming_alts[j].label].S2 = incoming_alts[j].S2;
            }
            incoming_alts.clear();
        }
        double sumS = 0.0;
        for(int i = 0; i < _record_surviving.size(); i++){
            sumS = sumS+sqrt (I[_record_surviving[i]].S2);
        }
        for(int i = 0; i < _record_surviving.size(); i++){
            I[_record_surviving[i]].batchsize = beta * sqrt(I[_record_surviving[i]].S2) * k / sumS;
        }
        for(int i = 1; i <= world_size - 1; i++){
            for(int j = 0; j < _record_surviving.size(); j++){
                if((j%(world_size-1)+1)==i){
                    outgoing_alts.push_back(I[_record_surviving[j]]);
                }
            }
            
            master_send_out_tasks(&outgoing_alts, i);
            outgoing_alts.clear();
        }
        _record_surviving.clear();
        
        for(int i = 0; i < world_size-1; i++){
            master_receive_tasks(&incoming_alts,&incoming_rank);
            for(int j = 0; j < incoming_alts.size(); j++){
                _record_surviving.push_back(incoming_alts[j].label);
            }
            incoming_alts.clear();
        }
        
        num_of_surviving = (int)_record_surviving.size();
        for(int i = 0; i < world_size-1; i++){
            MPI_Send(&num_of_surviving,1, MPI_INT, i+1, 0, MPI_COMM_WORLD);
        }
        //printf("Number of surviving alt is:%lu\n",_record_surviving.size());
    }else{
       
        worker_receive_tasks(&incoming_alts);
        for(int i = 0; i < incoming_alts.size(); i++){
            _tt=MPI_Wtime();
            for(int count = 0; count < n0; count++){
                generate_obv(&obv, incoming_alts[i].label);
                incoming_alts[i].sampleSize++;
                incoming_alts[i].sampleXSum = incoming_alts[i].sampleXSum + obv;
                incoming_alts[i].S2 = incoming_alts[i].S2 + obv * obv;
            }
            incoming_alts[i].S2 = (incoming_alts[i].S2-incoming_alts[i].sampleXSum*incoming_alts[i].sampleXSum/n0)/(n0-1);
            total_sim_time = total_sim_time + MPI_Wtime() - _tt;
        }
        worker_send_out_tasks(&incoming_alts);
        incoming_alts.clear();
        
        
        worker_receive_tasks(&incoming_alts);
        _tt = MPI_Wtime();
        for(int i = 0; i < incoming_alts.size();i++){
            bool eli = false;
            for(int j = 0; j < incoming_alts.size();j++){
                double tij = 1/(incoming_alts[i].S2/incoming_alts[i].sampleSize+incoming_alts[j].S2/incoming_alts[j].sampleSize);
                double Y = tij*(incoming_alts[i].sampleXSum/incoming_alts[i].sampleSize-incoming_alts[j].sampleXSum/incoming_alts[j].sampleSize);
                double tauij=1/(incoming_alts[i].S2/(n0+incoming_alts[i].batchsize*maxR)+incoming_alts[j].S2/(n0+incoming_alts[i].batchsize*maxR));
                double aij = eta*sqrt((n0-1)*tauij);
                if(i!=j&& Y<-aij){
                    eli =true;
                    break;
                }
            }
            if(eli==false){
                outgoing_alts.push_back(incoming_alts[i]);
            }
        }
        total_comparison_time = total_comparison_time + MPI_Wtime() - _tt;
        worker_send_out_tasks(&outgoing_alts);
        incoming_alts.clear();
        outgoing_alts.clear();
        
        MPI_Recv(&num_of_surviving,1, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    
    while(num_of_surviving > 1 && r < maxR){
        
        
        //generate obv.
        if(world_rank == 0){
            if(num_of_surviving <= world_size - 1){
                for(int i = 0; i < world_size - 1; i++){
                    outgoing_alts.push_back(I[_record_surviving[i]]);
                    master_send_out_tasks(&outgoing_alts,i+1);
                    outgoing_alts.clear();
                }
                for(int i = 0; i < world_size - 1; i++){
                    master_receive_tasks(&incoming_alts,&incoming_rank);
                    I[incoming_alts[0].label].sampleSize = incoming_alts[0].sampleSize;
                    I[incoming_alts[0].label].sampleXSum = incoming_alts[0].sampleXSum;
                    incoming_alts.clear();
                }
            }else{
                for(int i = 0; i < world_size - 1; i++){
                    outgoing_alts.push_back(I[_record_surviving[i]]);
                    master_send_out_tasks(&outgoing_alts,i+1);
                    outgoing_alts.clear();
                }
                
                int count = world_size - 1;
                while(count < num_of_surviving + world_size - 1){
                    master_receive_tasks(&incoming_alts,&incoming_rank);
                    I[incoming_alts[0].label].sampleSize = incoming_alts[0].sampleSize;
                    I[incoming_alts[0].label].sampleXSum = incoming_alts[0].sampleXSum;
                    incoming_alts.clear();
                    if(count < num_of_surviving){
                        outgoing_alts.push_back(I[_record_surviving[count]]);
                        master_send_out_tasks(&outgoing_alts,incoming_rank);
                        outgoing_alts.clear();
                    }else{
                        alt pseudoalt;
                        pseudoalt.label = -1;
                        outgoing_alts.push_back(pseudoalt);
                        master_send_out_tasks(&outgoing_alts,incoming_rank);
                        outgoing_alts.clear();
                    }
                    count++;
                }
            }
        }else{
            if(num_of_surviving <= world_size - 1){
                worker_receive_tasks(&incoming_alts);
                _tt = MPI_Wtime();
                for(int i = 0; i < incoming_alts.size(); i++){
                    for(int count = 0; count < incoming_alts[i].batchsize; count++){
                        generate_obv(&obv, incoming_alts[i].label);
                        incoming_alts[i].sampleSize++;
                        incoming_alts[i].sampleXSum = incoming_alts[i].sampleXSum + obv;
                    }
                }
                total_sim_time = total_sim_time + MPI_Wtime()-_tt;
                worker_send_out_tasks(&incoming_alts);
                incoming_alts.clear();
            }else{
                int _label_verify = 1;
                while(_label_verify >= 0){
                    worker_receive_tasks(&incoming_alts);
                    _label_verify = incoming_alts[0].label;
                    if(_label_verify >= 0){
                        _tt =MPI_Wtime();
                        for(int i = 0; i < incoming_alts.size(); i++){
                            for(int count = 0; count < incoming_alts[i].batchsize; count++){
                                generate_obv(&obv, incoming_alts[i].label);
                                incoming_alts[i].sampleSize++;
                                incoming_alts[i].sampleXSum = incoming_alts[i].sampleXSum + obv;
                            }
                        }
                        total_sim_time = total_sim_time + MPI_Wtime() -_tt;
                        worker_send_out_tasks(&incoming_alts);
                        incoming_alts.clear();
                    }else{
                        incoming_alts.clear();
                        break;
                    }
                }
            }
        }
        
        
        //screening;
        if(world_rank == 0){
            if(num_of_surviving <= world_size){
                _tt = MPI_Wtime();
                for(int i = 0; i < _record_surviving.size();i++){
                    bool eli = false;
                    for(int j = 0; j < _record_surviving.size();j++){
                        double tij = 1/(I[_record_surviving[i]].S2/I[_record_surviving[i]].sampleSize+I[_record_surviving[j]].S2/I[_record_surviving[j]].sampleSize);
                        double Y = tij*(I[_record_surviving[i]].sampleXSum/I[_record_surviving[i]].sampleSize-I[_record_surviving[j]].sampleXSum/I[_record_surviving[j]].sampleSize);
                        double tauij=1/(I[_record_surviving[i]].S2/(n0+I[_record_surviving[i]].batchsize*maxR)+I[_record_surviving[j]].S2/(n0+I[_record_surviving[j]].batchsize*maxR));
                        double aij = eta*sqrt((n0-1)*tauij);
                        if(i!=j&& Y<-aij){
                            eli =true;
                            break;
                        }
                    }
                    if(eli==false){
                        outgoing_alts.push_back(I[_record_surviving[i]]);
                    }
                }
                total_comparison_time = total_comparison_time + MPI_Wtime() - _tt;
                _record_surviving.clear();
                for(int i = 0 ; i < outgoing_alts.size();i++){
                    _record_surviving.push_back(outgoing_alts[i].label);
                }
                outgoing_alts.clear();
            }else{
                int max_index = 0;
                double _save_max = I[_record_surviving[0]].sampleXSum/I[_record_surviving[0]].sampleSize;
		int endding_ = (int)(floor(_record_surviving.size()/48.0)+1);
                for(int i = 1; i < endding_;i++){
                    if(I[_record_surviving[i]].sampleXSum/I[_record_surviving[i]].sampleSize>_save_max){
                        max_index = i;
                        _save_max = I[_record_surviving[i]].sampleXSum/I[_record_surviving[i]].sampleSize;
                    }
                }
                int temp_save = _record_surviving[max_index];
                _record_surviving[max_index] = _record_surviving[0];
                _record_surviving[0] = temp_save;
                
                int groupSize = (int)floor((_record_surviving.size() - 1)/world_size);
                for(int i = 1; i < world_size; i++){
                    outgoing_alts.push_back(I[_record_surviving[0]]);
                    int stoppingpoint = (i + 1)* groupSize;
                    if(i == world_size - 1){
                        stoppingpoint = _record_surviving.size()-1;
                    }
                    for(int j = i * groupSize+1;j<=stoppingpoint;j++){
                        outgoing_alts.push_back(I[_record_surviving[j]]);
                    }
                    
                    
                    /**for(int j = 1; j < _record_surviving.size(); j++){
                        if(j% world_size==i){
                            outgoing_alts.push_back(I[_record_surviving[j]]);
                        }
                    }**/
                    
                    master_send_out_tasks(&outgoing_alts, i);
                    outgoing_alts.clear();
                }
                for(int j = 0; j <= groupSize; j ++){
                    outgoing_alts.push_back(I[_record_surviving[j]]);
                }
                
                /**for(int j = 0; j < _record_surviving.size(); j++){
                    if(j% world_size==0){
                        outgoing_alts.push_back(I[_record_surviving[j]]);
                    }
                }**/
                
                _record_surviving.clear();
                _tt = MPI_Wtime();
                for(int i = 0; i < outgoing_alts.size();i++){
                    bool eli = false;
                    for(int j = 0; j < outgoing_alts.size();j++){
                        double tij = 1/1/(outgoing_alts[i].S2/outgoing_alts[i].sampleSize+outgoing_alts[j].S2/outgoing_alts[j].sampleSize);
                        double Y = tij*(outgoing_alts[i].sampleXSum/outgoing_alts[i].sampleSize-outgoing_alts[j].sampleXSum/outgoing_alts[j].sampleSize);
                        double tauij=1/(outgoing_alts[i].S2/(n0+outgoing_alts[i].batchsize*maxR)+outgoing_alts[j].S2/(n0+outgoing_alts[j].batchsize*maxR));
                        double aij = eta*sqrt((n0-1)*tauij);
                        if(i!=j&& Y<-aij){
                            eli =true;
                            break;
                        }
                    }
                    if(eli==false){
                        _record_surviving.push_back(outgoing_alts[i].label);
                    }
                }
                total_comparison_time = total_comparison_time + MPI_Wtime() - _tt;
                outgoing_alts.clear();
                
                for(int i = 0; i < world_size-1; i++){
                    master_receive_tasks(&incoming_alts,&incoming_rank);
                    for(int j = 0; j < incoming_alts.size(); j++){
                        _record_surviving.push_back(incoming_alts[j].label);
                    }
                    incoming_alts.clear();
                }
                
                num_of_surviving = (int)_record_surviving.size();
                
                for(int i = 0; i < world_size-1; i++){
                    MPI_Send(&num_of_surviving,1, MPI_INT, i+1, 0, MPI_COMM_WORLD);
                }
		/**if(r==maxR-1){
			double firstStageSampleSize = 0.0;
			
			for(int u=0; u < num_of_surviving; u++){
				firstStageSampleSize = firstStageSampleSize + I[_record_surviving[u]].sampleSize;			
			}
			printf("r=%d, Number of surviving alt is:%d and FirstStageSampleSize:%f\n",r+1,num_of_surviving,firstStageSampleSize);
		}**/
            }
        }else{
            if(num_of_surviving > world_size){
                worker_receive_tasks(&incoming_alts);
                _tt = MPI_Wtime();
                for(int i = 1; i < incoming_alts.size();i++){
                    bool eli = false;
                    for(int j = 0; j < incoming_alts.size();j++){
                        double tij = 1/(incoming_alts[i].S2/incoming_alts[i].sampleSize+incoming_alts[j].S2/incoming_alts[j].sampleSize);
                        double Y = tij*(incoming_alts[i].sampleXSum/incoming_alts[i].sampleSize-incoming_alts[j].sampleXSum/incoming_alts[j].sampleSize);
                        double tauij=1/(incoming_alts[i].S2/(n0+incoming_alts[i].batchsize*maxR)+incoming_alts[j].S2/(n0+incoming_alts[i].batchsize*maxR));
                        double aij = eta*sqrt((n0-1)*tauij);
                        if(i!=j&& Y<-aij){
                            eli =true;
                            break;
                        }
                    }
                    if(eli==false){
                        outgoing_alts.push_back(incoming_alts[i]);
                    }
                }
                total_comparison_time = total_comparison_time + MPI_Wtime() - _tt;
                worker_send_out_tasks(&outgoing_alts);
                incoming_alts.clear();
                outgoing_alts.clear();
                MPI_Recv(&num_of_surviving,1, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
        }
        
        
        r++;
    }
    
    if(num_of_surviving == 1){
        if(world_rank == 0){
            bestsys = I[_record_surviving[0]].label;
        }
        for(int i = 0; i < I.size();i++){
            totalSampleSize = totalSampleSize+I[i].sampleSize;
        }
    }else{
        if(world_rank == 0){
            if(num_of_surviving <= world_size - 1){
                for(int i = 0; i < world_size - 1; i++){
                    outgoing_alts.push_back(I[_record_surviving[i]]);
                    master_send_out_tasks(&outgoing_alts,i+1);
                    outgoing_alts.clear();
                }
                for(int i = 0; i < world_size - 1; i++){
                    master_receive_tasks(&incoming_alts,&incoming_rank);
                    I[incoming_alts[0].label].sampleSize = incoming_alts[0].sampleSize;
                    I[incoming_alts[0].label].sampleXSum = incoming_alts[0].sampleXSum;
                    incoming_alts.clear();
                }
            }else{
                for(int i = 0; i < world_size - 1; i++){
                    outgoing_alts.push_back(I[_record_surviving[i]]);
                    master_send_out_tasks(&outgoing_alts,i+1);
                    outgoing_alts.clear();
                }
                
                int count = world_size - 1;
                while(count < num_of_surviving + world_size - 1){
                    master_receive_tasks(&incoming_alts,&incoming_rank);
                    I[incoming_alts[0].label].sampleSize = incoming_alts[0].sampleSize;
                    I[incoming_alts[0].label].sampleXSum = incoming_alts[0].sampleXSum;
                    incoming_alts.clear();
                    if(count < num_of_surviving){
                        outgoing_alts.push_back(I[_record_surviving[count]]);
                        master_send_out_tasks(&outgoing_alts,incoming_rank);
                        outgoing_alts.clear();
                    }else{
                        alt pseudoalt;
                        pseudoalt.label = -1;
                        outgoing_alts.push_back(pseudoalt);
                        master_send_out_tasks(&outgoing_alts,incoming_rank);
                        outgoing_alts.clear();
                    }
                    count++;
                }
            }
            double _maxavg = I[_record_surviving[0]].sampleXSum/I[_record_surviving[0]].sampleSize;
            bestsys=I[_record_surviving[0]].label;
            for(int i = 0; i < I.size();i++){
                totalSampleSize = totalSampleSize+I[i].sampleSize;
                
            }
            for(int i = 1; i < _record_surviving.size();i++){
                if(I[_record_surviving[0]].sampleXSum/I[_record_surviving[0]].sampleSize > _maxavg){
                    _maxavg = I[_record_surviving[0]].sampleXSum/I[_record_surviving[0]].sampleSize ;
                    bestsys = I[_record_surviving[0]].label;
                }
                
            }
            //printf("total sample size is:%i\n",totalSampleSize);
            //printf("best system is:%i\n",bestsys);
            //printf("Wall clock time is:%f\n",MPI_Wtime()-_starttime);
            
        }else{
            if(num_of_surviving <= world_size - 1){
                worker_receive_tasks(&incoming_alts);
                _tt = MPI_Wtime();
                for(int i = 0; i < incoming_alts.size(); i++){
                    int maxN = ceil(incoming_alts[i].S2 * rinotth * rinotth / (delta*delta));
                    for(int count = incoming_alts[i].sampleSize; count < maxN; count++){
                        generate_obv(&obv, incoming_alts[i].label);
                        incoming_alts[i].sampleSize++;
                        incoming_alts[i].sampleXSum = incoming_alts[i].sampleXSum + obv;
                    }
                }
                total_sim_time = total_sim_time + MPI_Wtime() - _tt;
                worker_send_out_tasks(&incoming_alts);
                incoming_alts.clear();
            }else{
                int _label_verify = 1;
                while(_label_verify >= 0){
                    worker_receive_tasks(&incoming_alts);
                    _label_verify = incoming_alts[0].label;
                    if(_label_verify >= 0){
                        _tt = MPI_Wtime();
                        for(int i = 0; i < incoming_alts.size(); i++){
                            int maxN = ceil(incoming_alts[i].S2 * rinotth * rinotth / (delta*delta));
                            for(int count = incoming_alts[i].sampleSize; count < maxN; count++){
                                generate_obv(&obv, incoming_alts[i].label);
                                incoming_alts[i].sampleSize++;
                                incoming_alts[i].sampleXSum = incoming_alts[i].sampleXSum + obv;
                            }
                        }
                        total_sim_time = total_sim_time + MPI_Wtime() - _tt;
                        worker_send_out_tasks(&incoming_alts);
                        incoming_alts.clear();
                    }else{
                        incoming_alts.clear();
                        break;
                    }
                }
            }
        }
    }
    
    if(world_rank > 0){
        MPI_Send(&total_sim_time,1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        MPI_Send(&total_comparison_time,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
    }
    
    if(world_rank == 0){
        double store_time;
        for(int i = 0; i < world_size-1; i++){
            int incoming_rank;
            MPI_Status status1;
            MPI_Probe(MPI_ANY_SOURCE, 1,MPI_COMM_WORLD,&status1);
            
            incoming_rank = status1.MPI_SOURCE;
            MPI_Recv(&store_time,1, MPI_DOUBLE, incoming_rank, 1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            total_sim_time = total_sim_time + store_time;
            
            
            MPI_Probe(MPI_ANY_SOURCE, 2,MPI_COMM_WORLD,&status1);
            
            incoming_rank = status1.MPI_SOURCE;
            MPI_Recv(&store_time,1, MPI_DOUBLE, incoming_rank, 2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            total_comparison_time = total_comparison_time + store_time;
        }
        printf("%d, %f, %f, %f, %f \n",count1+1,totalSampleSize,MPI_Wtime()-_starttime,total_sim_time,total_comparison_time);
        ofstream myfile ("GSPoutput.txt");
        myfile<<bestsys<<" "<<totalSampleSize<<" "<<MPI_Wtime()-_starttime<<" "<<total_sim_time<<" "<<total_comparison_time<<"\n";
        
        //    myfile<<"The wall clock time is: "<<-_starttime+_endtime<<"\n"<<"The best alternative is: "<<bestLabel<<"\n"<<"The total sample size is: "<<totalSampleSize<<"\n"<<"Total comparison time is: "<<total_comparison_time<<"\n";
        myfile.close();
    }
    
    }
    MPI_Finalize();
    return 0;

}
