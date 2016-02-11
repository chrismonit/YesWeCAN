/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class ProbMatrixFactory {
    
    /**
     * To change the algorithm used to compute P, change the ProbMatrixGenerator 
     * type returned by this method
     */
    public static ProbMatrixGenerator getPGenerator(RateMatrix Q){
        return new YangEigenDecomposition(Q);
        //return new PtSeriesExpansion(Q, 4);

    }

}//class
