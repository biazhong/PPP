//
//  PPP-PACNew.cpp
// This code implements the PPP in MPI.
// The simulation optimization problem considered in this code is the three-stage buffer allocation problem.
// The parameters follows by comments "Input Parameter:..." should be adjusted from one problem instance to another.
// While generating observations, no pause time is considered in this code.
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
    int sampleSize;
    double sampleXSum;
    double S2;
    int batchsize;
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
    int RB = 128;   //Input Parameter: Problem parameter
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
    int k = 1016127;    //Input Parameter: Total number of alternatives
    int n0 = 50;    //Input Parameter: First-stage sample size
    double alpha = 0.05;    //Input Parameter: Desired PAC
    double delta = 0.1;     //Input Parameter: IZ parameter delta
    double lambda = delta / 2;
    double h2 = (n0-1)/(4*(delta-lambda))*(pow(alpha/(k-1), -2.0/(n0-1))-1);
    
    int _survivingK = k;
    double maxN =0.0 ;
    int t = n0;
    
    
    double total_comparison_time = 0.0;
    double totalSampleSize = 0.0;
    double total_sim_time = 0.0;
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
    
    
    double obv=0.0;
    vector<alt> I;
    vector<int> _record_surviving;
    if(world_rank == 0){
        _starttime = MPI_Wtime();
        for (int i = 0; i <k; i++){
	    alt tempAlt;
            tempAlt.label = i;
            tempAlt.sampleSize = 0;
            tempAlt.sampleXSum = 0.0;
            tempAlt.S2 = 0.0;
	    tempAlt.batchsize = 0;
            I.push_back(tempAlt);
            _record_surviving.push_back(i);
        }
    }
    
    if(world_rank == 0){
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
            master_receive_tasks(&incoming_alts);
            for(int j = 0; j < incoming_alts.size(); j++){
                I[incoming_alts[j].label].sampleSize = incoming_alts[j].sampleSize;
                I[incoming_alts[j].label].sampleXSum = incoming_alts[j].sampleXSum;
                I[incoming_alts[j].label].S2 = incoming_alts[j].S2;
            }
            incoming_alts.clear();
        }
    }else{
        worker_receive_tasks(&incoming_alts);
        _tt = MPI_Wtime();
        for(int i = 0; i < incoming_alts.size(); i++){
            for(int count = 0; count < n0; count++){
                generate_obv(&obv, incoming_alts[i].label);
                incoming_alts[i].sampleSize++;
                incoming_alts[i].sampleXSum = incoming_alts[i].sampleXSum + obv;
                incoming_alts[i].S2 = incoming_alts[i].S2 + obv * obv;
            }
            incoming_alts[i].S2 = (incoming_alts[i].S2 - incoming_alts[i].sampleXSum* incoming_alts[i].sampleXSum / n0)/(n0-1);
        }
        worker_send_out_tasks(&incoming_alts);
        incoming_alts.clear();
        total_sim_time = total_sim_time + MPI_Wtime() -  _tt;
    }
    
    if(world_rank == 0){
        double maxS2=0;
        int maxS2INDEX = -1;
        
        
        for(int i = 0; i < k; i++){
            if(I[i].S2 > maxS2){
                maxS2 = I[i].S2;
                maxS2INDEX = i;
            }
        }
        
        double sec_maxS2 = 0;
        int sec_maxS2INDEX = -1;
        for(int i = 0; i < k; i++){
            if(I[i].S2 > sec_maxS2 && i!=maxS2INDEX){
                sec_maxS2 = I[i].S2;
                sec_maxS2INDEX = i;
            }
        }
        
        maxN = floor((maxS2 + sec_maxS2) * h2/lambda) + 1;

	printf("max Index : %d\n", maxS2INDEX);
	printf("max S2: %f\n", maxS2);

	printf("second max Index: %d\n", sec_maxS2INDEX);
	printf("second max S2: %f\n", sec_maxS2);        

        printf("maxN is :%f\n",maxN);
        
        for(int i = 1; i <= world_size - 1; i++){
            MPI_Send(&maxN,1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
        /**
         printf("First Stage Finished,%f\n",MPI_Wtime()-_starttime);
         printf("maxS2 is %f\n",maxS2);
         printf("Corresponding average is: %f\n",I[maxS2INDEX].sampleXSum/n0);
         printf("h2 is %f\n",h2);
         printf("maxN is %f\n",maxN);**/
    }else{
        MPI_Recv(&maxN,1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    
    while(_survivingK > world_size - 1 && t <= maxN){
        if(world_rank == 0){
            double _tt = MPI_Wtime();
            int i_star = 0;
            double temp = I[_record_surviving[0]].sampleXSum - h2*I[_record_surviving[0]].S2;
            for(int i = 1; i < _record_surviving.size();i++){
                if(I[_record_surviving[i]].sampleXSum - h2 * I[_record_surviving[i]].S2>temp){
                    temp = I[_record_surviving[i]].sampleXSum - h2 * I[_record_surviving[i]].S2;
                    i_star = i;
                }
            }
            int swapTemp = _record_surviving[i_star];
            _record_surviving[i_star] = _record_surviving[0];
            _record_surviving[0] = swapTemp;
            
            for(int i = 1; i < _record_surviving.size();i++){
                if(I[_record_surviving[i]].sampleXSum + h2 * I[_record_surviving[i]].S2 < I[_record_surviving[0]].sampleXSum - h2 * I[_record_surviving[0]].S2-(delta-lambda) * I[_record_surviving[0]].sampleSize){
                    _record_surviving[i] = _record_surviving[_record_surviving.size()-1];
                    _record_surviving.erase(_record_surviving.begin()+_record_surviving.size()-1);
                    
                    i--;
                }
            }
            total_comparison_time = total_comparison_time + MPI_Wtime() - _tt;
            _survivingK = _record_surviving.size();
            
            
            
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
                master_receive_tasks(&incoming_alts);
                for(int j = 0; j < incoming_alts.size(); j++){
                    I[incoming_alts[j].label].sampleSize = incoming_alts[j].sampleSize;
                    I[incoming_alts[j].label].sampleXSum = incoming_alts[j].sampleXSum;
                }
                incoming_alts.clear();
            }
            
            for(int i = 1; i <= world_size - 1; i++){
                MPI_Send(&_survivingK,1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
            t++;
            //printf("Number of Surviving alt is: %lu\n",_record_surviving.size());
        }else{
            worker_receive_tasks(&incoming_alts);
            _tt = MPI_Wtime();
            for(int i = 0; i < incoming_alts.size(); i++){
                generate_obv(&obv, incoming_alts[i].label);
                incoming_alts[i].sampleSize++;
                incoming_alts[i].sampleXSum = incoming_alts[i].sampleXSum + obv;
            }
            worker_send_out_tasks(&incoming_alts);
            incoming_alts.clear();
            total_sim_time = total_sim_time + MPI_Wtime()- _tt;
            
            MPI_Recv(&_survivingK,1, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            t++;
        }
    }
    if(t <= maxN){
        
        if(world_rank == 0){
		double maxS2=0;
        	int maxS2INDEX = -1;
        
        
        	for(int i = 0; i < _record_surviving.size(); i++){
            		if(I[_record_surviving[i]].S2 > maxS2){
                		maxS2 = I[_record_surviving[i]].S2;
                		maxS2INDEX = i;
            		}
        	}
        
        	double sec_maxS2 = 0;
        	int sec_maxS2INDEX = -1;
        	for(int i = 0; i < _record_surviving.size(); i++){
            		if(I[_record_surviving[i]].S2 > sec_maxS2 && i!=maxS2INDEX){
                		sec_maxS2 = I[_record_surviving[i]].S2;
               			sec_maxS2INDEX = i;
            		}
        	}
	        

        	maxN = floor((maxS2 + sec_maxS2) * h2/lambda) + 1;

		printf("max Index : %d\n", _record_surviving[maxS2INDEX]);
		printf("max S2: %f\n", maxS2);

		printf("second max Index: %d\n",_record_surviving[sec_maxS2INDEX]);
		printf("second max S2: %f\n", sec_maxS2);   
		printf	("maxN is : %f\n",maxN);

            //printf("Final Stage\n");
	    int batch = maxN - t;
	    if(batch<0){
	    	batch = 0;
	    }
            for(int i = 0; i < _record_surviving.size(); i++){
		I[_record_surviving[i]].batchsize = batch;
                outgoing_alts.push_back(I[_record_surviving[i]]);
                master_send_out_tasks(&outgoing_alts, i+1);
                outgoing_alts.clear();
            }
            
            for(int i = 0; i < _record_surviving.size(); i++){
                master_receive_tasks(&incoming_alts);
                for(int j = 0; j < incoming_alts.size(); j++){
                    I[incoming_alts[j].label].sampleSize = incoming_alts[j].sampleSize;
                    I[incoming_alts[j].label].sampleXSum = incoming_alts[j].sampleXSum;
                }
                incoming_alts.clear();
            }
        }else if(world_rank <= _survivingK){
            worker_receive_tasks(&incoming_alts);
            _tt = MPI_Wtime();
            for(int i = 0; i < incoming_alts[0].batchsize; i++){
                generate_obv(&obv, incoming_alts[0].label);
                incoming_alts[0].sampleSize++;
                incoming_alts[0].sampleXSum = incoming_alts[0].sampleXSum + obv;
            }
            worker_send_out_tasks(&incoming_alts);
            incoming_alts.clear();
            total_sim_time = total_sim_time + MPI_Wtime() - _tt;
        }else{
            
        }
    }
    if(world_rank > 0){
        MPI_Send(&total_sim_time,1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }
    
    
    if(world_rank==0){
        double store_simtime;
        for(int i = 0; i < world_size-1; i++){
            int incoming_rank;
            MPI_Status status1;
            MPI_Probe(MPI_ANY_SOURCE, 1,MPI_COMM_WORLD,&status1);
            
            incoming_rank = status1.MPI_SOURCE;
            MPI_Recv(&store_simtime,1, MPI_DOUBLE, incoming_rank, 1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            total_sim_time = total_sim_time + store_simtime;
            
        }
        
        
        
        
        int bestLabel = I[_record_surviving[0]].label;
        double tempBest = I[_record_surviving[0]].sampleXSum;
        for(int i = 1; i < _record_surviving.size(); i++){
            if(I[_record_surviving[i]].sampleXSum > tempBest){
                tempBest = I[_record_surviving[i]].sampleXSum;
                bestLabel = I[_record_surviving[i]].label;
            }
        }
        for(int i = 0; i < I.size();i++){
            totalSampleSize = totalSampleSize + I[i].sampleSize;
        }
        
        _endtime = MPI_Wtime();
        
        
        ofstream myfile ("PPPPACoutput.txt");
        myfile<<bestLabel<<" "<<totalSampleSize<<" "<<-_starttime+_endtime<<" "<<total_sim_time<<" "<<total_comparison_time<<"\n";
        
        //    myfile<<"The wall clock time is: "<<-_starttime+_endtime<<"\n"<<"The best alternative is: "<<bestLabel<<"\n"<<"The total sample size is: "<<totalSampleSize<<"\n"<<"Total comparison time is: "<<total_comparison_time<<"\n";
        myfile.close();
    }
    MPI_Finalize();
    return 0;
}
