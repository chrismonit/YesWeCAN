/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model;

import org.apache.commons.math3.linear.RealMatrix;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public interface ProbMatrixGenerator {
    
    
    public RealMatrix getP(double t);
    
    public RateMatrix getQ();
}
