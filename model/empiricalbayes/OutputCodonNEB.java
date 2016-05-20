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
import yeswecan.utils.ArrayPrinter;

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
    
    // make new method which gets all the contiguous sites in a gene
    // i.e. splice it
    // then have a method which pulls out the codons sites
    
    public static int[] intArray(List<Integer> intList){
        int[] array = new int[intList.size()];
        for (int i = 0; i < intList.size(); i++) {
            array[i] = intList.get(i);
        }
        return array;
    }
    
    public int[] getGeneContiguousSites(int gene){
        List<Integer> geneSites = new ArrayList<Integer>();
        
        for (int iSite = 0; iSite < this.genStruct.getTotalLength(); iSite++) {
            int partition = this.genStruct.getPartitionIndex(iSite);
            
            if (this.genStruct.containsGene(partition, gene)) {
                geneSites.add(iSite);
            }
        }//for iSite
        
        return intArray(geneSites);
    }
    
    // if array length is not multiple of three, will just ignore the last 1 or 2 elements
    public int[][] getCodons(int[] sites){
        
        int nCodons = sites.length/3; // integer division. Only considers number of complete triplets
        
        int[][] codons = new int[nCodons][3];
        
        for (int iSite = 0; iSite < sites.length-2; iSite+=3) {
            int codonIndex = iSite/3;
            codons[codonIndex] = new int[]{ sites[iSite], sites[iSite+1], sites[iSite+2] };
        }
        return codons;
    }
    
    
    
    
    
    // this is redundant now I think
    public int[] getCodonSites(int gene){
        
        List<Integer> codonStartSites = new ArrayList<Integer>();
        
        // -2 because you can't have a codon starting right at the end of the sequence
        for (int iSite = 0; iSite < this.genStruct.getTotalLength()-2; iSite++) {
            
            int partition = this.genStruct.getPartitionIndex(iSite);
            int siteType = iSite%3;
            
            if (this.genStruct.containsGene(partition, gene)) {
                
                int frame = this.genStruct.getFrame(partition, gene);

                // two sites up from here, is the same gene present in the same frame?
                // if so, we have a complete codon
                int plusTwoPartition = this.genStruct.getPartitionIndex(iSite+2);
                boolean completeCodon = this.genStruct.genePresent(gene, plusTwoPartition, frame);
                
                // if this is leading site in codon. eg if frame==0 && siteType==0 then this is alpha site in frame A, 
                // and is a first codon position by definition
                if (frame == siteType && completeCodon) {
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
