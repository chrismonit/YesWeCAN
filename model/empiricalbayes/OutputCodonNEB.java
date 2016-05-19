/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.empiricalbayes;

import java.util.ArrayList;
import java.util.List;
import pal.datatype.CodonTable;
import pal.tree.Tree;
import yeswecan.model.submodels.CANModelFrequenciesMix;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;

/**
 *
 * @author cmonit1
 */
public class OutputCodonNEB {
    
    protected AdvancedAlignment alignment;
    protected Tree tree;
    protected GeneticStructure genStruct;
    protected CANModelFrequenciesMix canModel;
    
    protected CodonFrequencies[] codonFrequenciesArray;
    protected CodonTable codonTable;
    
    protected int numSiteClasses;
    
    protected CodonNaiveEmpiricalBayesCalculator codonNEB;
    
    public OutputCodonNEB(
            AdvancedAlignment alignment, Tree tree, 
            GeneticStructure genStruct, CANModelFrequenciesMix canModel,
            CodonFrequencies[] codonFrequenciesArray, CodonTable codonTable,
            int numSiteClasses
    ){
        this.alignment = alignment;
        this.tree = tree;
        this.genStruct = genStruct;
        
        this.canModel = canModel;
        // NB 0th omega is fixed to 1.0 for neutral evolution
           
        this.codonFrequenciesArray = codonFrequenciesArray;
        this.codonTable = codonTable;
        
        this.numSiteClasses = numSiteClasses;
        
        this.codonNEB = 
                new CodonNaiveEmpiricalBayesCalculator(
                        this.alignment, this.tree, this.genStruct, 
                        this.canModel, this.numSiteClasses);
        
    }
    

    public int[] getCodonSites(int gene){
        
        List<Integer> codonStartSites = new ArrayList<Integer>();
        
        // -2 because you can't have a codon starting right at the end of the sequence
        for (int iSite = 0; iSite < this.genStruct.getTotalLength()-2; iSite++) {
            
            int partition = this.genStruct.getPartitionIndex(iSite);
            int siteType = iSite%3;
            
            if (this.genStruct.containsGene(partition, gene)) {
                
                int frame = this.genStruct.getFrame(partition, gene);
                
//                int plusTwoPartition = this.genStruct.getPartitionIndex(iSite+2);
//                boolean completeCodon = this.genStruct.
                
                // if this is leading site in codon. eg if frame==0 && siteType==0 then this is alpha site in frame A, 
                // and is a first codon position by definition
                if (frame == siteType) {
                    codonStartSites.add(iSite);
                }
                                
            }// if gene present in this partition
            
        }// iSite
    
        int[] codonStartSitesArray = new int[codonStartSites.size()];
        for (int i = 0; i < codonStartSites.size(); i++) {
            codonStartSitesArray[i] = codonStartSites.get(i);
        }
        return codonStartSitesArray;
        
    }// getCodonSites
    
    
}
