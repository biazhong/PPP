  //This class implements the PPP used to conduct the experiments in Section 7.1 of 
//the paper entitled Speeding Up Paulson's Procedure for Large-Scale Problem Using Parallel Computing

package comparePaulsonandModifiedPaulson;

import java.util.ArrayList;
import java.util.Random;

public class ModifiedPaulson {
	private int k,n0;
	private double delta,lambda,alpha;
	private final int bestK = 1;
	private final double sigma = 0.5;
	private Random R = new Random();
	public ModifiedPaulson() {
		this.k = 10;
		this.n0 = 10;
		this.delta = 0.01;
		this.lambda = 0.005;
		this.alpha = 0.05;
	}
	
	public ModifiedPaulson(int k, int n0, double delta,double lambda,double alpha) {
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
		//System.out.println(h2);
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
		
		double[] S = new double[k];
		
		varianceCalculation(S, X, sampleMean);
		
		int t = n0;
		while(I.size()>1) {
			int i_star = 0;
			double temp = t*sampleMean[I.get(0)]-h2*S[I.get(0)];
			for(int i  = 1 ; i < I.size(); i++) {
				if(t*sampleMean[I.get(i)]-h2*S[I.get(i)]>temp) {
					temp = t*sampleMean[I.get(i)]-h2*S[I.get(i)];
					i_star = i;
				}
			}
			int tempIndex = I.get(i_star);
			I.set(i_star, I.get(0));
			I.set(0, tempIndex);
			
			
			for(int i  = 1 ; i < I.size(); i++) {
				if( t*sampleMean[I.get(i)]+h2*S[I.get(i)]<t*sampleMean[I.get(0)]-h2*S[I.get(0)]+lambda*t) {
					I.set(i, I.get(I.size()-1));
					I.remove(I.size()-1);
					totalSampleSize = totalSampleSize + t;
					i--;
				}
			}
			
			for(int i = 1; i < I.size();i++) {
				if( t*sampleMean[I.get(0)]+h2*S[I.get(0)]<t*sampleMean[I.get(i)]-h2*S[I.get(i)]+lambda*t) {
					I.set(0, I.get(I.size()-1));
					I.remove(I.size()-1);
					totalSampleSize = totalSampleSize + t;
					break;
				}
			}
			
			
			
			
			updateSampleMean(t, I, mu, sampleMean);
			t++;
		}
		//System.out.println(t);
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
	
	protected void varianceCalculation(double[] S, double[][] X, double[] sampleMean) {
		int k = S.length;
		int nSamples = X.length;
		
		for(int i = 0 ; i < k ; i++) {
			for (int count = 0 ; count < nSamples; count++) {
				S[i] = S[i] + Math.pow((X[count][i]-sampleMean[i]), 2.0);
			}
			S[i]=S[i]/(nSamples-1);
			//System.out.println(S[i]);
		}
	}
	
	protected void updateSampleMean(int t, ArrayList<Integer> I, ArrayList<Double> mu, double[] sampleMean) {
		for(int i = 0; i < I.size(); i++) {
			sampleMean[I.get(i)] = (sampleMean[I.get(i)] * t + R.nextGaussian()*sigma + mu.get(I.get(i)))/(t+1);
		}
	}
}
