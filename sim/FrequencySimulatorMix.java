/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.sim;

import java.util.List;
import java.util.Random;
import pal.datatype.Codons;
import pal.tree.Tree;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.Probabilities;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;
import yeswecan.phylo.ReorderFrequencies;
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
                genStruct.getTotalLength(), genStruct.getNumberOfGenes()+1, 
                rand, probabilities);
        
    }
    
    public int[][] getGenSiteClasses(){
        return this.geneSiteClasses;
    }
    
    
    @Override
    public double computeRate(int[] quintStates, int j, int site, double nu){
        int[] genes = this.genStruct.getGenes(site);
        double product = 1.0; // multiplicative identity
                
        product *= this.kappa.getKappaIfTransition(quintStates[2], j);
        for (int iFrame = 0; iFrame < 3; iFrame++) {
            int[] codonI_array = getCodon(quintStates, quintStates[2], iFrame, site%3);
            int codonI = Codons.getCodonIndexFromNucleotideStates(codonI_array);
            
            int[] codonJ_array = getCodon(quintStates, j, iFrame, site%3);

            int codonJ = Codons.getCodonIndexFromNucleotideStates(codonJ_array);
               
            if (!this.codonTable.isSynonymous(codonI, codonJ)) { 
                int gene = genes[iFrame];
                int siteClass = this.geneSiteClasses[site][gene];
                
                int omegaIndex = (gene * this.numSiteClasses) + siteClass;
                
                product *= this.omegas.get(omegaIndex).get();
            }

            CodonFrequencies geneCodonFreq = this.codonFrequencies.get(genes[iFrame]);
            int[] mappedToPaml = ReorderFrequencies.alphaToPaml(codonJ_array); // expecting codons will be ordered TCAG in codonFrequencies instances
            double pi_J = geneCodonFreq.getFrequency(mappedToPaml); 
            product *= pi_J;
            
        }// iFrame
        
        product *= nu;
        
        return product;
        
    }
    
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
