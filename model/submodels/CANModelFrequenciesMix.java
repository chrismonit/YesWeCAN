/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.submodels;

import java.util.List;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.Probabilities;
import yeswecan.model.parameters.TsTvRatioAdvanced;

/**
 *
 * @author cmonit1
 */
public class CANModelFrequenciesMix extends CANModelFrequencies {
    
    private List<Probabilities> probabilities;
    private int numSiteClasses; // used for accessing omegas
    
    public CANModelFrequenciesMix(
        TsTvRatioAdvanced kappa,
        BranchScaling scaling, List<Omega> omegas,
        List<Probabilities> probabilities, int numSiteClasses
    ){
    
        super(kappa, scaling, omegas);
        
        this.probabilities = probabilities;
        this.numSiteClasses = numSiteClasses;
        
        for (Probabilities p : this.probabilities){
            super.addParameters(p);
        }
    }
    

    public Omega getGeneAndSiteClassOmega(int gene, int siteClass){
        int index = gene * this.numSiteClasses + siteClass;
        return super.getOmegas().get(index);
    }
    
    public double getProbability(int gene, int siteClass){
        return this.probabilities.get(gene).get()[siteClass];
    }
    
    public Probabilities getGeneProbabilities(int gene){
        return this.probabilities.get(gene);
    }
    
    public List<Probabilities> getProbabilities(){
        return this.probabilities;
    }
}
