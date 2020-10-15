//
//  VKN.cpp
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
#include <string>

using namespace std;

double exprand(double lambda){
    double u;
    
    u = (double)rand() / (RAND_MAX * 1.0+1);
    
    return -log(1- u) / lambda;
}


typedef struct{
    int label;
    int position;
    double sim_obv;
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
    
    
    int RB = 50;    //Input Parameter: Problem parameter
    vector<int> x_disc(5);
    int rr = RB * 2 - 3;
    int n = label/(RB-1)+1;
    
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
    
    
    int k = 57624;  //Input Parameter: Total number of alternatives
    int n0 = 50;    //Input Parameter: First-stage sample size
    double alpha = 0.05;    //Input Parameter: Desired PAC
    double delta = 0.1;     //Input Parameter: IZ parameter delta
    double h2 = (n0-1)*(pow(2.0*alpha/(k-1),-2.0/(n0-1))-1);
    int incoming_rank;
    int _survivingk = k;
    double maxN = 0.0;
    int t = 0;
    double total_comparison_time = 0.0;
    double totalSampleSize = 0.0;
    double total_simulation_time = 0.0;
        
    double _tt = 0.0;
    double _starttime = 0.0, _endtime;
    MPI_Init(NULL,NULL);
    
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    srand(time(NULL)*world_rank);
    vector<alt> outgoing_alts;
    vector<alt> incoming_alts;
    alt tempalt;
    
    double obv=0.0;
    vector<int> order_position;
    vector<int> _record_surviving;
    
    
    vector<double> S1(k);
    vector<vector<double> > S2(k,S1);
    vector<double> AVG;
    
    vector<double> Xl(k);
    vector<vector<double> > X(40000,Xl);
    int sim_index;
    
    if(world_rank == 0){
        _starttime = MPI_Wtime();
        for(int i = 0; i < k; i++){
            order_position.push_back(0);
            _record_surviving.push_back(i);
            AVG.push_back(0.0);
        }
        
        for(int i = 0; i < world_size - 1; i++){
            sim_index = i;
            tempalt.label = _record_surviving[sim_index];
            tempalt.position = order_position[tempalt.label];
            order_position[tempalt.label]++;
            outgoing_alts.push_back(tempalt);
            master_send_out_tasks(&outgoing_alts,i+1);
            outgoing_alts.clear();
        }
    }
    
    while(_survivingk > 1){
        if(world_rank == 0){
            master_receive_tasks(&incoming_alts,&incoming_rank);
            X[incoming_alts[0].position][incoming_alts[0].label]=incoming_alts[0].sim_obv;
            incoming_alts.clear();
            sim_index++;
            sim_index = sim_index % (int)_record_surviving.size();
            int record_simlabel = _record_surviving[sim_index];
            bool check = true;
            
            
            
            for(int i = 0; i < _record_surviving.size(); i++){
                if(X[t][_record_surviving[i]] == 0){
                    check = false;
                    break;
                }
            }
            
            if(check == true){
                t++;
                _tt = MPI_Wtime();
                if(t == n0){
                    for(int i = 0; i < k; i++){
                        for(int j = 0 ; j < n0; j++){
                            AVG[i] = AVG[i]+X[j][i];
                        }
                        AVG[i] = AVG[i] / n0;
                    }
                    for(int i = 0; i < k; i++){
                        for(int j = i+1; j < k; j++){
                            double _diff_sum = 0.0;
                            double _diff_sum2 = 0.0;
                            for(int count = 0; count < n0; count++){
                                _diff_sum = _diff_sum+ X[count][i]-X[count][j];
                                _diff_sum2 =_diff_sum2+ (X[count][i]-X[count][j])*(X[count][i]-X[count][j]);
                            }
                            S2[i][j] = (_diff_sum2 - _diff_sum * _diff_sum/n0)/(n0-1);
                            S2[j][i] = S2[i][j];
                        }
                    }
                }
                
                if(t >= n0){
                    if(t > n0 - 1){
                        for (int i = 0 ; i < _record_surviving.size();i++){
                            AVG[_record_surviving[i]]=(AVG[_record_surviving[i]] * (t-1) + X[t-1][_record_surviving[i]])/t;
                        }
                    }
                    for(int i = 0; i < _record_surviving.size();i++){
                        bool eli = false;
                        for(int j = 0 ; j < _record_surviving.size();j++){
                            if(i!=j && AVG[_record_surviving[i]]-AVG[_record_surviving[j]]<min(0.0,-h2*S2[_record_surviving[i]][_record_surviving[j]]/(2*t*delta)+delta/2)){
                                eli = true;
                                break;
                            }
                        }
                        if(eli == true){
                            _record_surviving.erase(_record_surviving.begin()+i);
                            i--;
                        }
                    }
                }
                total_comparison_time = total_comparison_time + MPI_Wtime() -_tt;
                _survivingk =  _record_surviving.size();
               // printf("current time is %i\n",t);
               // printf("current surviving is %lu\n",_record_surviving.size());
                for(int i = 0; i < _record_surviving.size();i++){
                    if(_record_surviving[i] >= record_simlabel){
                        sim_index = i;
                        record_simlabel = _record_surviving[i];
                        break;
                    }
                }
            }
            
            if(_survivingk == 1){
                tempalt.label = -1;
                outgoing_alts.push_back(tempalt);
                master_send_out_tasks(&outgoing_alts,incoming_rank);
                vector<int> recordincoming;
                recordincoming.push_back(incoming_rank);
                for(int i = 0; i < world_size - 2;i++){
                    master_receive_tasks(&incoming_alts,&incoming_rank);
                    bool checkEli = true;
                    for(int count = 0; count < recordincoming.size();count++){
                        if(incoming_rank == recordincoming[count]){
                            checkEli = false;
                            break;
                        }
                    }
                    if(checkEli == true){
                        incoming_alts.clear();
                        master_send_out_tasks(&outgoing_alts,incoming_rank);
                    }else{
                        i--;
                        incoming_alts.clear();
                    }
                }
                outgoing_alts.clear();
            }else{
                tempalt.label = record_simlabel;
                tempalt.position = order_position[record_simlabel];
                order_position[record_simlabel]++;
                outgoing_alts.push_back(tempalt);
                master_send_out_tasks(&outgoing_alts,incoming_rank);
                outgoing_alts.clear();
            }
        }else{
            worker_receive_tasks(&incoming_alts);
            _tt = MPI_Wtime();
            if(incoming_alts[0].label == -1){
                _survivingk = 1;
            }else{
                generate_obv(&obv, incoming_alts[0].label);
                incoming_alts[0].sim_obv = obv;
                worker_send_out_tasks(&incoming_alts);
                incoming_alts.clear();
                totalSampleSize++;
            }
            total_simulation_time = total_simulation_time + MPI_Wtime()-_tt;
        }
    }
    if(world_rank > 0){
        MPI_Send(&total_simulation_time,1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        
        MPI_Send(&totalSampleSize,1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    }
    
    if(world_rank == 0){
        double store_simtime;
        double store_samplesize;
        for(int i = 0; i < world_size-1; i++){
        
            MPI_Status status1;
            MPI_Probe(MPI_ANY_SOURCE, 1,MPI_COMM_WORLD,&status1);
        
            incoming_rank = status1.MPI_SOURCE;
            MPI_Recv(&store_simtime,1, MPI_DOUBLE, incoming_rank, 1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            total_simulation_time = total_simulation_time + store_simtime;
            
            MPI_Probe(MPI_ANY_SOURCE, 2,MPI_COMM_WORLD,&status1);
            
            incoming_rank = status1.MPI_SOURCE;
            MPI_Recv(&store_samplesize,1, MPI_DOUBLE, incoming_rank, 2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            totalSampleSize = totalSampleSize + store_samplesize;
            
            
        }
        ofstream myfile ("VKNoutput.txt");
        myfile<<_record_surviving[0]<<" "<<totalSampleSize<<" "<<MPI_Wtime()-_starttime<<" "<<total_simulation_time<<" "<<total_comparison_time<<"\n";
        //printf("Total sample size is: %i\n",totalSampleSize);
        //printf("Wall clock time is: %f\n",MPI_Wtime()-_starttime);
        //printf("Total comparison time is: %f\n",total_comparison_time);
        myfile.close();
    }
    MPI_Finalize();
    return 0;
}
