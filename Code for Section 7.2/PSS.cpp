//
//  PSS.cpp
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
    int sampleSize;
    int batch_size;
    double sampleXSum;
    double S2;
    int eli;
    double budget;
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

void master_receive_tasks_main(vector<alt>* incoming_alts, int* outgoing_rank){
    MPI_Status status;
    int incoming_rank;
    int incoming_size;
    
    MPI_Probe(MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
    MPI_Get_count(&status, MPI_BYTE,&incoming_size);
    incoming_rank = status.MPI_SOURCE;
    
    *outgoing_rank = incoming_rank;
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
    int RB = 20;    //Input Parameter: Problem parameter
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
	int k = 3249;
	double ref_mu=5.776121751463529; //Input Parameter: Highest mu
	int n0 = 50;    //Input Parameter: First-stage sample size
	double alpha = 0.05;    //Input Parameter: Desired PAC
	double delta = 0.1;     //Input Parameter: IZ parameter delta
	double budget = 410000; //Input Parameter: Total sample size
	double c=4.35;          //Input Parameter: Constant c used to construct continuation regions. NOTE: For different problems, the value of c could be different.
	double batch = 10;
    
	int _survivingK = k;
	int t = n0;
    
    	int outgoing_rank;
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
			tempAlt.batch_size = n0;
            		tempAlt.sampleXSum = 0.0;
            		tempAlt.S2 = 0.0;
			tempAlt.eli = 0;
			tempAlt.budget = budget;
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
        	budget = budget - n0*k;
		
        	for(int i = 0; i < world_size-1; i++){
			
            		master_receive_tasks(&incoming_alts);

			
			           		
			for(int j = 0; j < incoming_alts.size(); j++){
                		I[incoming_alts[j].label].sampleSize = incoming_alts[j].sampleSize;
                		I[incoming_alts[j].label].sampleXSum = incoming_alts[j].sampleXSum;
                		I[incoming_alts[j].label].S2 = incoming_alts[j].S2;
			}
			
            	}
		incoming_alts.clear();
			            	
		
		_tt = MPI_Wtime();
		
		for(int i = 0; i < _record_surviving.size();i++){
			double tau = t/I[_record_surviving[i]].S2;
			if(tau*(I[_record_surviving[i]].sampleXSum/I[_record_surviving[i]].sampleSize-ref_mu)< sqrt((c+log(tau+1))*(tau+1))*(-1)){
				_record_surviving[i] = _record_surviving[_record_surviving.size()-1];
        	            	_record_surviving.erase(_record_surviving.begin()+_record_surviving.size()-1);
        	            	i--;	
			}

        	}
		_survivingK = _record_surviving.size();

		
		total_comparison_time = total_comparison_time + MPI_Wtime() - _tt;

    	}else{

		worker_receive_tasks(&incoming_alts);
        	_tt = MPI_Wtime();

		
        	for(int i = 0; i < incoming_alts.size(); i++){
            		for(int count = 0; count < incoming_alts[i].batch_size; count++){
                		generate_obv(&obv, incoming_alts[i].label);
                		incoming_alts[i].sampleXSum = incoming_alts[i].sampleXSum + obv;
                		incoming_alts[i].S2 = incoming_alts[i].S2 + obv * obv;
            		}
			incoming_alts[i].sampleSize = incoming_alts[i].sampleSize + incoming_alts[i].batch_size;
			incoming_alts[i].batch_size = 0;
            		incoming_alts[i].S2 = (incoming_alts[i].S2 - incoming_alts[i].sampleXSum* incoming_alts[i].sampleXSum / n0)/(n0-1);
        	}
		
        	worker_send_out_tasks(&incoming_alts);
        	incoming_alts.clear();
        	total_sim_time = total_sim_time + MPI_Wtime() -  _tt;

        	
    	}

	if(world_rank == 0){
		
		for(int i = 1; i <= world_size - 1; i++){
			int round_batch = batch;
			if(budget < batch){
				round_batch = budget;		
			}

			I[_record_surviving[0]].batch_size = round_batch;
			
			budget = budget - round_batch;
			
	
			
            		outgoing_alts.push_back(I[_record_surviving[0]]);
			_record_surviving.erase(_record_surviving.begin());
            		
            		master_send_out_tasks(&outgoing_alts, i);
            		outgoing_alts.clear();
        	}

		int leftWorkers = 47;

		while(budget > 0){
			master_receive_tasks_main(&incoming_alts,&outgoing_rank);
			if(incoming_alts[0].label >= 0){
				for(int j = 0; j < incoming_alts.size(); j++){
                			I[incoming_alts[j].label].sampleSize = incoming_alts[j].sampleSize;
                			I[incoming_alts[j].label].sampleXSum = incoming_alts[j].sampleXSum;
                			I[incoming_alts[j].label].S2 = incoming_alts[j].S2;
					I[incoming_alts[j].label].eli = incoming_alts[j].eli;
					
					if(I[incoming_alts[j].label].eli == 0){
						_record_surviving.push_back(incoming_alts[j].label);
					}
				}			
			}
			incoming_alts.clear();

			if(_record_surviving.size()>=1){

			int round_batch = batch;
			if(budget < batch){
				round_batch = budget;		
			}
			
			I[_record_surviving[0]].batch_size = round_batch;
			
			budget = budget - round_batch;

			
			
			
			
			

			incoming_alts.push_back(I[_record_surviving[0]]);
			_record_surviving.erase(_record_surviving.begin());

			master_send_out_tasks(&incoming_alts, outgoing_rank);
            		incoming_alts.clear();
			double si = _record_surviving.size();
			printf("Q remaining: %f\n",si);
			}else{
			alt tempAlt;
			tempAlt.label = -1;
            		tempAlt.sampleSize = 0;
			tempAlt.batch_size = 0;
            		tempAlt.sampleXSum = 0.0;
            		tempAlt.S2 = 0.0;
			tempAlt.eli = 0;
			tempAlt.budget = 0;
			incoming_alts.push_back(tempAlt);
			master_send_out_tasks(&incoming_alts, outgoing_rank);
			leftWorkers--;
			incoming_alts.clear();
			}
			printf("Budget remaining: %f\n",budget);				
		}
		
	
		for(int i = 1; i <= leftWorkers ; i++){
			master_receive_tasks_main(&incoming_alts,&outgoing_rank);
			if(incoming_alts[0].label >= 0){
				for(int j = 0; j < incoming_alts.size(); j++){
                			I[incoming_alts[j].label].sampleSize = incoming_alts[j].sampleSize;
                			I[incoming_alts[j].label].sampleXSum = incoming_alts[j].sampleXSum;
                			I[incoming_alts[j].label].S2 = incoming_alts[j].S2;
					I[incoming_alts[j].label].eli = incoming_alts[j].eli;
					if(I[incoming_alts[j].label].eli == 0){
						_record_surviving.push_back(incoming_alts[j].label);
					}
				}			
			}

			incoming_alts.clear();	
			alt tempAlt;
			tempAlt.label = -1;
            		tempAlt.sampleSize = 0;
			tempAlt.batch_size = 0;
            		tempAlt.sampleXSum = 0.0;
            		tempAlt.S2 = 0.0;
			tempAlt.eli = 0;
			tempAlt.budget = 0;
			incoming_alts.push_back(tempAlt);
			master_send_out_tasks(&incoming_alts, outgoing_rank);
			incoming_alts.clear();
		}
			
	}else{
		while(budget > 0){
			worker_receive_tasks(&incoming_alts);
            		_tt = MPI_Wtime();
            		for(int i = 0; i < incoming_alts.size(); i++){
				budget = incoming_alts[i].budget;				
				for(int count = 0; count < incoming_alts[i].batch_size; count++){
					
                			generate_obv(&obv, incoming_alts[i].label);
                			incoming_alts[i].sampleXSum = incoming_alts[i].sampleXSum + obv;
				}
				incoming_alts[i].sampleSize = incoming_alts[i].sampleSize + incoming_alts[i].batch_size;
				incoming_alts[i].batch_size = 0;
            		}

		
			
			total_sim_time = total_sim_time + MPI_Wtime()- _tt;	


			_tt = MPI_Wtime();
			double tau = incoming_alts[0].sampleSize/incoming_alts[0].S2;
			
			

			if(tau*(incoming_alts[0].sampleXSum/incoming_alts[0].sampleSize-ref_mu)< sqrt((c+log(tau+1))*(tau+1))*(-1)){
				incoming_alts[0].eli = 1;	
			}
	
			total_comparison_time = total_comparison_time + MPI_Wtime() - _tt;
			if(budget > 0){
            		worker_send_out_tasks(&incoming_alts);
			}
            		incoming_alts.clear();
			
		}

	}

	if(world_rank > 0){
        	MPI_Send(&total_sim_time,1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		MPI_Send(&total_comparison_time,1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    	}
	
	if(world_rank==0){
        	double store_simtime;
        	for(int i = 0; i < 2*(world_size-1); i++){
        		int incoming_rank;
        		MPI_Status status1;
        		MPI_Probe(MPI_ANY_SOURCE, 1,MPI_COMM_WORLD,&status1);
            
        		incoming_rank = status1.MPI_SOURCE;
        		MPI_Recv(&store_simtime,1, MPI_DOUBLE, incoming_rank, 1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			if(store_simtime <0.01){
		    		total_comparison_time = total_comparison_time + store_simtime;	
			}else{
        	    		total_sim_time = total_sim_time + store_simtime;
			}
        	}
       		_endtime = MPI_Wtime();
        	
        	
        	ofstream myfile ("PSSoutput.txt");
		myfile<<totalSampleSize<<" "<<-_starttime+_endtime<<" "<<total_sim_time<<" "<<total_comparison_time<<"\n";
		for(int count =0; count < _record_surviving.size();count++){
			myfile<< I[_record_surviving[count]].label<<"\n";			
		}
        	myfile.close();
    	}
	MPI_Finalize();
	return 0;
}
