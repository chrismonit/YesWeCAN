/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.parameters;

import org.apache.commons.math3.util.MathArrays;
import swmutsel.model.parameters.Probabilities;
import swmutsel.utils.CoreUtils;
import yeswecan.Constants;

/**
 *
 * @author cmonit1
 */
public class ProbabilitiesFixLast extends Probabilities {
    
    private double[] probsFixLast;
    
    protected static double fixValue = 0.0;
    
   
    public ProbabilitiesFixLast(double[] probabilities){
        super(probabilities); 
        // this calls set method in superclass.
        // that in turn asigns normalised values to its probabilities field.
        // but that field is private.
        // Therefore the outcome of calling the super constructor has no effect
        // on this class
        
        set(probabilities);
    }// constructor
    
    
    private static void checkInput(double[] input, double fixValue){
        double lastProb = input[input.length-1];
        if (lastProb > fixValue + Constants.EPSILON || lastProb < fixValue - Constants.EPSILON) {
            throw new RuntimeException("ProbabilitiesFixLast: last element is not close to required value ("+fixValue+")");
        }
    }
    
    // Note that every method in Probabilities superclass is overriden
    
    @Override
    public double[] get() {
        return this.probsFixLast;
    }
    
    @Override
    public void set(double[] x) {
        checkInput(x, this.fixValue);
        this.probsFixLast = MathArrays.normalizeArray(x, 1); 
        // if last element is close to zero, normalising the whole array will still leave the last element close to zero
    }
    
    // NOTE mius 2, because last element is fixed and the remaining n-1 elements must sum to 1
    
    @Override
    public double[] getOptimisable() {
        double[] x = new double[probsFixLast.length - 2]; // see NOTE, above
        System.arraycopy(this.probsFixLast, 0, x, 0, x.length);
        return CoreUtils.alr(x); // return the optimisible params in opt space
    }

    @Override
    public void setOptimisable(double[] params) {
        double[] in = CoreUtils.alr_inv(params);
        double sum = 0;
        for (int i = 0; i < in.length; i++) {
            this.probsFixLast[i] = in[i];
            sum += in[i];
        }
        this.probsFixLast[this.probsFixLast.length - 2] = 1.0 - sum;
        this.probsFixLast[this.probsFixLast.length - 1] = this.fixValue;
    }

    @Override
    public int getOptimisableCount() {
        return this.probsFixLast.length - 2; // see NOTE, above
    }

    @Override
    public String toString() {
        return "ProbabilitiesFixLast{" +
                "p=" +
                CoreUtils.join("%.7f", ", ", this.probsFixLast) +
                '}';
    }
    
    

    
    
}
