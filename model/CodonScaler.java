/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model;

import yeswecan.sim.CodonFrequencies;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class CodonScaler implements RatioScaler {
    
    
    public CodonScaler(CodonFrequencies codonFrequencies){
    
    
    }
    
    // ordered AAA to TTT
    // includes stop codons, with their values set to 0.0
    double[] nonSynFractions = new double[]{0.875, 0.8888888888888888, 0.875, 0.8888888888888888, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.75, 0.8888888888888888, 0.7777777777777778, 0.8888888888888888, 0.7777777777777778, 0.7777777777777778, 1.0, 0.7777777777777778, 0.875, 0.8888888888888888, 0.875, 0.8888888888888888, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.5, 0.6666666666666666, 0.5555555555555556, 0.6666666666666666, 0.5555555555555556, 0.6666666666666666, 0.5555555555555556, 0.6666666666666666, 0.875, 0.8888888888888888, 0.875, 0.8888888888888888, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.625, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.6666666666666666, 0.0, 0.8571428571428571, 0.0, 0.8571428571428571, 0.5714285714285714, 0.6666666666666666, 0.625, 0.0, 0.6666666666666666, 0.875, 1.0, 0.875, 0.7142857142857143, 0.8888888888888888, 0.75, 0.8888888888888888 };
    
    // copied from ProportionScaler
    //TODO: should put in RatioScaler interface as static thing to be used by all
    
    private int[][] codonPositions = new int[][]{
            // a frame, b frame, c frame
            { 0, 2, 1 }, // alpha site
            { 1, 0, 2 }, // beta site
            { 2, 1, 0 }  // gamma site
        };
    
    
    
    // TODO
//    public double get(double ratio, int siteType, int frame){
//        
//        int codonPos = codonPositions[siteType][frame];
//        double probNonSyn = probNonsyn[codonPos];
//        double probSyn = 1.0 - probNonSyn;
//        return (probSyn * 1.0) + (probNonSyn * ratio);
//    }
    
    
    
}
