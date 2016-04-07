/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model.ratioscaling;

import yeswecan.Constants;
import yeswecan.utils.ArrayPrinter;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class ProportionScaler implements RatioScaler{
    /*
        Which probability of bening nonsynonymous to use depends on both the site type {alpha, beta, gamma} and the frame {a,b,c}
        codonPositions gives you the nucleotide position within the codon
            e.g. if the frame is a and the site type is alpha, the codon position is 1 (0 in zero based)
        The index is then used to access probNonSyn for the probability of a change at that site being nonsynonymous
    */
    
    
    private int[][] codonPositions = new int[][]{
            // a frame, b frame, c frame
            { 0, 2, 1 }, // alpha site
            { 1, 0, 2 }, // beta site
            { 2, 1, 0 }  // gamma site
        };
    
    // probabilities of changes at the 1st, 2nd and 3rd codon positions resulting in nonsynonymous change
    // derived from counting 
    private double[] probNonsyn = new double[]{ 0.954022988506, 1.0, 0.284090909091 };
    
    
    
    public ProportionScaler(){
        // empty constructor
    };
    
    
    public double get(double ratio, int siteType, int frame){
        
        int codonPos = codonPositions[siteType][frame];
        double probNonSyn = probNonsyn[codonPos];
        double probSyn = 1.0 - probNonSyn;
        return (probSyn * 1.0) + (probNonSyn * ratio);
    }
    
    @Override
    public String toString(){
        return ArrayPrinter.toString(probNonsyn, Constants.DEL);
    }
    
}
