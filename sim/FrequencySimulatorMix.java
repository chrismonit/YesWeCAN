/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.sim;

import java.util.List;
import java.util.Random;
import pal.tree.Tree;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.Probabilities;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;
import yeswecan.phylo.States;

/**
 *
 * @author cmonit1
 */
public class FrequencySimulatorMix extends FrequencySimulator {
    
    private List<Probabilities> probabilities;
    private int numSiteClasses;
    
    private int[][] geneSiteClasses;
    
    public FrequencySimulatorMix(Tree tree, Random rand, GeneticStructure genStruct,
        TsTvRatioAdvanced kappa, List<Omega> omegas, List<CodonFrequencies> codonFrequencies,
        List<Probabilities> probabilities, int numSiteClasses
    
        ){

        super(tree, rand, genStruct, kappa, omegas, codonFrequencies);
        
        this.probabilities = probabilities;
        this.numSiteClasses = numSiteClasses; // do we need this?
        
        this.geneSiteClasses = assignGeneSiteClasses(
                genStruct.getTotalLength(), genStruct.getNumberOfGenes(), 
                rand, probabilities);
        
        // print the matrix so we know the true classes at each site
    }
    
//    @Override
//    public double computeRate(int[] quintStates, int j, int site, double nu){
//        int[] genes = super.genStruct.getGenes(site);
//        double product = 1.0; // multiplicative identity
//                
//        product *= this.kappa.getKappaIfTransition(quintStates[2], j);
//        
//        
//    }
    
    public static int[][] assignGeneSiteClasses(int numSites, int numGenes, 
            Random rand, List<Probabilities> probabilities){
        
        int[][] siteClasses = new int[numSites][numGenes];
        
        for (int iSite = 0; iSite < numSites; iSite++) {
            for (int iGene = 0; iGene < numGenes; iGene++) {
                
                double[] probs = probabilities.get(iGene).get();
                int siteClass = States.draw(probs, rand.nextDouble());
                siteClasses[iSite][iGene] = siteClass;
            }
        }
        
        return siteClasses;
    }

}
