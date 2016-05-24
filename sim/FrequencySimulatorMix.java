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
    
    private static final int TERMINAL_SITES = 2; // number of sites we want to skip at end of sequence (2 because there are 3 sites in codon)
    private static final int NUM_FRAMES = 3; 
    private static final int CODON_LENGTH = 3; 
    
    
    
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
    
    public int[][] getGeneSiteClasses(){
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
    
    public static int[][] assignCodonSiteClasses( 
            Random rand, List<Probabilities> probabilities, GeneticStructure genStruct){
        
        int[][] siteClasses = new int[genStruct.getTotalLength()][NUM_FRAMES];
        
        // set initial values in array to -1 so we don't mistake an un-assigned element for site class 0
        for (int iSite = 0; iSite < siteClasses.length; iSite++) {
            for (int iFrame = 0; iFrame < siteClasses[0].length; iFrame++) {
                siteClasses[iSite][iFrame] = -1;
            }
        }
        
        for (int iFrame = 0; iFrame < NUM_FRAMES; iFrame++) {
            
            for (int iSite = iFrame; iSite < genStruct.getTotalLength()-TERMINAL_SITES; iSite += CODON_LENGTH) {
                
                int[] genes = genStruct.getGenes(iSite);
                int gene = genes[iFrame];
                double[] probs = probabilities.get(gene).get();
                int codonSiteClass = States.draw(probs, rand.nextDouble());
                for (int iCodonPosition = 0; iCodonPosition < 3; iCodonPosition++) {
                    siteClasses[iSite+iCodonPosition][iFrame] = codonSiteClass;
                }// iCodonPosition
            }// iSite
        }// iFrame
        return siteClasses;
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
