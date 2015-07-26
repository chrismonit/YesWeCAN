/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class PtSeriesExpansion implements ProbMatrixGenerator {
    
    
    private RateMatrix Q;
    private int highestTerm;
    
    public PtSeriesExpansion(RateMatrix Q, int highestTerm){
        this.Q = Q;
        this.highestTerm = highestTerm;
    }
    
    private static double factorial(int n){
        double product = 1.0;
        for (int i = n; i > 0 ; i--) {
            product *= i;
        }
        return product;
    }
    
    public RealMatrix getP(double t){
        return seriesExpansion(this.Q, t, this.highestTerm);
    }
    
    
    private static RealMatrix seriesExpansion(RateMatrix Q, double t, int highestTerm){
        
        /*
            exp(Qt) = I + Qt + ((Qt)^2)/2! + ((Qt)^3)/3! + ... + ((Qt)^n)/n!
        */
        
        RealMatrix Qt = Q.scalarMultiply(t);
        
        // make identity matrix
        double[][] identity = new double[Q.getRowDimension()][Q.getColumnDimension()];
        for (int i = 0; i < Q.getRowDimension(); i++) {
            identity[i][i] = 1.0;
        }
                
        RealMatrix sum = new Array2DRowRealMatrix(identity); // sum == I
        
        sum = sum.add( Qt ); // sum == I + Qt
        
        for (int n = 2; n <= highestTerm; n++) { // sum = I + Qt + ((Qt)^2)/2! + ((Qt)^3)/3! + ... + ((Qt)^n)/n!
                        
            double[][] termData = Qt.power(n).scalarMultiply(1.0/factorial(n)).getData();
                        
            RealMatrix term = new Array2DRowRealMatrix(termData);
            sum = sum.add(term);
            
        }// for
        
        return sum;
        
    }// seriesExpansion
    
    
    
    
}
