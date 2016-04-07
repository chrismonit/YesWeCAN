/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model.submodels;

import java.util.List;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.Probabilities;
import yeswecan.model.submodels.CANModel;
import yeswecan.model.parameters.TsTvRatioAdvanced;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class CANModelMixture extends CANModel {

    private List<Probabilities> probabilities;
    int numSiteClasses; // used for accessing omegas
    
    public CANModelMixture(HKYModel hky, BranchScaling scaling, List<Omega> omegas, 
            List<Probabilities> probabilities, int numSiteClasses){
        this(hky.getKappa(), hky.getPi(), scaling, omegas, probabilities, numSiteClasses);
    }
     
    public CANModelMixture(TsTvRatioAdvanced kappa, BaseFrequencies pi, BranchScaling scaling, 
            List<Omega> omegas, List<Probabilities> probabilities, int numSiteClasses){
        // NB the 0th omega has to be an unoptimisible 1.0 value
        
        super(kappa, pi, scaling, omegas);               
        
        this.probabilities = probabilities;
        this.numSiteClasses = numSiteClasses;
        
        for (Probabilities p : this.probabilities){
            super.addParameters(p);
        }
    }
    // NB this does NOT override parent method getOmega(int gene)
    // calling that method from this class would give you the wrong omega instance
    public Omega getOmega(int gene, int siteClass){
        int index = gene * this.numSiteClasses + siteClass;
        return super.getOmegas().get(index);
    }
    
    public double getProbability(int gene, int siteClass){
        return this.probabilities.get(gene).get()[siteClass];
    }
    
    public Probabilities getGeneProbabilities(int gene){
        return this.probabilities.get(gene);
    }
    
}
