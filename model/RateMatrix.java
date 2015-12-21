
package yeswecan.model;

import yeswecan.model.parameters.*;

import yeswecan.utils.MatrixPrinter;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import swmutsel.model.parameters.BaseFrequencies;

/**
 *
 * @author Christopher Monit (c.monit.12@ucl.ac.uk)
 * 
 * Constructs a rate matrix for the CAN model, from specified parameter values
 * 
 * 
 */
public class RateMatrix extends Array2DRowRealMatrix {

    private BaseFrequencies pi; 
    private TsTvRatioAdvanced kappa;
    
    public BaseFrequencies getBaseFrequencies(){ 
        return this.pi; 
    }
    
    public TsTvRatioAdvanced getKappa(){
        return this.kappa;
    }
    
    public RateMatrix( TsTvRatioAdvanced kappa, BaseFrequencies pi){ 
        //suitable for constructing a classic HKY85 Q matrix
        
        super(); // this will create a matrix with no data; populated at end of constructor
        this.kappa = kappa;
        this.pi = pi;
        
        
        int numStates = this.pi.get().length;
        
        double[][] matrixData = new double[numStates][numStates];
                
        //populate off-diagonal elements
        for (int i = 0; i < numStates; i++) {
            for (int j = 0; j < numStates; j++) {
                
                if (i == j) continue;
                matrixData[i][j] = this.pi.get()[j] * this.kappa.getKappaIfTransition(i, j);
            }// j
        }// i
        
        
        //populate diagonal elements
        
        for (int i = 0; i < numStates; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < numStates; j++) {
                rowSum += matrixData[i][j];
            }
            matrixData[i][i] = -rowSum;
        }
        
        //MatrixPrinter.PrintMatrix(matrixData, "matrix data before scaling");
        
        // scale matrix such that average substitution rate is 1
        
        // compute nu

        double nuDenominator = 0.0;
        
        for (int i = 0; i < numStates; i++) {
            nuDenominator += this.pi.get()[i] * matrixData[i][i];
        }
        double nu = 1.0 / -nuDenominator;
        
        // do the scaling
        for (int i = 0; i < numStates; i++) {
            for (int j = 0; j < numStates; j++) {
                matrixData[i][j] *= nu;
            }
        }
        
        //MatrixPrinter.PrintMatrix(matrixData, "matrix data after scaling");

        
        super.setSubMatrix(matrixData, 0, 0); //replace the (hitherto blank) matrix data
               
    }//constructor
    
    
    
    public static void main(String[] args){
        System.out.println("Testing RateMatrix");
        TsTvRatioAdvanced k = new TsTvRatioAdvanced(2.0);
        BaseFrequencies freqs = new BaseFrequencies(new double[]{ 0.25, 0.25, 0.25, 0.25 });
        
        RateMatrix Q = new RateMatrix( k, freqs );
        double[][] matrixData = Q.getData();
        MatrixPrinter.PrintMatrix(matrixData, "matrix data");

    }

}//class
