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
public abstract class NucEmpiricalBayesCalculator extends EmpiricalBayesCalculator {
     
    public abstract double[][][] getEBValues();
}
