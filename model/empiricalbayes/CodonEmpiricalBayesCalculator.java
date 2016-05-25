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
public abstract class CodonEmpiricalBayesCalculator extends EmpiricalBayesCalculator {
    
    public abstract double getNormalisationFactor(int[] codonSites);
    
    public abstract double getNumerator(int[] codonSites, int codonVSiteClass);

}
