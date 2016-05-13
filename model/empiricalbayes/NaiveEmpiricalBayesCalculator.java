/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.empiricalbayes;

/**
 *
 * @author cmonit1
 */
public class NaiveEmpiricalBayesCalculator extends EmpiricalBayesCalculator {
    
    protected double[][][] probValues;
    
    
    
    
    
    
    @Override
    public double[][][] getEBValues(){
        return this.probValues;
    }
    
}
