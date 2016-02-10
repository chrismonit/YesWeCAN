
package yeswecan.model;

import yeswecan.model.parameters.*;


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
    
    private int numStates;
    
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
        
        this.numStates = this.pi.get().length;
        
        super.setSubMatrix(new double[numStates][numStates], 0, 0); // makes a square zero matrix
        
        //populate off-diagonal elements
        for (int i = 0; i < numStates; i++) {
            for (int j = 0; j < numStates; j++) {
                
                if (i == j) continue;
                double q_ij = this.pi.get()[j] * this.kappa.getKappaIfTransition(i, j);
                this.setEntry(i,j, q_ij);
            }// j
        }// i
          
        this.populateDiagonals();
        this.scale();

    }//constructor
    
    
    public final void populateDiagonals(){
         for (int i = 0; i < numStates; i++) {
                double rowSum = 0.0;
                for (int j = 0; j < numStates; j++) {
                    if (i == j) continue;
                    rowSum += this.getEntry(i, j);
                }
                this.setEntry(i, i, -rowSum);
            }
    }
    
    
    public final void scale(){
        // scale matrix such that average substitution rate is 1
        // compute nu
        
        double nuDenominator = 0.0;
        
        for (int i = 0; i < this.numStates; i++) {
            nuDenominator += this.pi.get()[i] * this.getEntry(i, i);
        }
        double nu = 1.0 / -nuDenominator;
        
        // apply the scaling
        super.setSubMatrix(this.scalarMultiply(nu).getData(), 0, 0);
    }

    

}//class
