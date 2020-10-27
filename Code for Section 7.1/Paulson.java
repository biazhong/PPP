  
//This class implements the classical Paulson's procedure used to conduct the experiments in Section 7.1 of 
//the paper entitled Speeding Up Paulson's Procedure for Large-Scale Problem Using Parallel Computing
package comparePaulsonandModifiedPaulson;
import java.util.ArrayList;
import java.util.Random;


public class Paulson {
	private int k,n0;
	private double delta,lambda,alpha;
	private final int bestK = 1;
	private final double sigma = 0.5;
	private Random R = new Random();
	public Paulson() {
		this.k = 10;
		this.n0 = 10;
		this.delta = 0.01;
		this.lambda = 0.005;
		this.alpha = 0.05;
	}
	
	public Paulson(int k, int n0, double delta,double lambda,double alpha) {
		this.k = k;
		this.n0 = n0;
		this.delta = delta;
		this.lambda = lambda;
		this.alpha = alpha;
	}
	
	protected ArrayList<Integer> getResult() {
		ArrayList<Integer> RESULT = new ArrayList<Integer>();
		int totalSampleSize = 0;
		int selectedA = -1;
		
		
		final double h2 = (n0-1)/(4*(delta-lambda))*(Math.pow(alpha/(k-1), -2d/(n0-1))-1);
		ArrayList<Integer> I = new ArrayList<Integer>();
		ArrayList<Double> mu = new ArrayList<Double>();
		for(int i = 0 ; i < k;i++) {
			I.add(i);
			if(i==bestK) {
				mu.add(delta);
			}else {
				mu.add(0d);
			}
		}
		
		double[][] X = new double[n0][k];
		firstStageSampleGeneration(X,sigma,mu);
		
		
		double[] sampleMean = new double[k];
		meanCalculation(sampleMean,X);
		
		double[][] S = new double[k][k];
		
		//double[] S1 = new double[k];
		
		//double[][] COR = new double[k][k];
		
		varianceCalculation(S, X, sampleMean);
		//varianceCalculation1(S1,X,sampleMean);
		//corCalculation(COR,S1,X,sampleMean);
		
		int t = n0;
		while(I.size()>1) {
			for(int i  = 0 ; i < I.size(); i++) {
				for(int j = 0 ; j < I.size(); j++) {
				//	if(I.get(i)!=1&t*(sampleMean[I.get(i)]-sampleMean[1])<-1d*h2*S[I.get(i)][1]+lambda*t) {
				
					if(i!=j & t*(sampleMean[I.get(i)]-sampleMean[I.get(j)])<-1d*h2*S[I.get(i)][I.get(j)]+lambda*t) {
						
						
						I.set(i, I.get(I.size()-1));
						I.remove(I.size()-1);
						totalSampleSize = totalSampleSize + t;
						i--;
						break;
					}
				}
			}
			
			updateSampleMean(t, I, mu, sampleMean);
			t++;
			
		}
		totalSampleSize = totalSampleSize + t;
		selectedA =  I.get(0);
		
		if(selectedA==bestK) {
			RESULT.add(1);
		}else {
			RESULT.add(0);
		}
		RESULT.add(totalSampleSize);
		return RESULT;
	}
	
	protected void firstStageSampleGeneration(double[][] X, double sigma, ArrayList<Double> mu) {
		int size = X[0].length;
		int nSamples = X.length;
		for(int i = 0; i < nSamples; i++) {
			for (int j = 0; j < size; j++) {
				X[i][j]= R.nextGaussian()*sigma+mu.get(j);
			}
		}
	}
	
	protected void meanCalculation(double[] sampleMean, double[][] X) {
		int size = X[0].length;
		int nSamples = X.length;
		for(int i = 0 ; i < size;i++) {
			for(int j = 0 ; j < nSamples;j++) {
				sampleMean[i]=sampleMean[i] + X[j][i];
			}
			sampleMean[i] = sampleMean[i]/nSamples;
		}
		
	}
	
	protected void varianceCalculation(double[][] S, double[][] X, double[] sampleMean) {
		int k = S.length;
		int nSamples = X.length;
		
		for(int i = 0 ; i < k ; i++) {
			for(int j = i; j < k; j++) {
				for (int count = 0 ; count < nSamples; count++) {
					S[i][j] = S[i][j] + Math.pow((X[count][i]-sampleMean[i])-(X[count][j]-sampleMean[j]), 2);
				}
				S[i][j] =S[i][j]/(nSamples-1);
				S[j][i]=S[i][j];
			}
		}
	}
	

	protected void updateSampleMean(int t, ArrayList<Integer> I, ArrayList<Double> mu, double[] sampleMean) {
		for(int i = 0; i < I.size(); i++) {
			sampleMean[I.get(i)] = (sampleMean[I.get(i)] * t + R.nextGaussian()*sigma + mu.get(I.get(i)))/(t+1);
		}
	}
}
