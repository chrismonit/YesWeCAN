/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.io;

import java.util.ArrayList;
import java.util.List;
import yeswecan.Constants;
import yeswecan.phylo.GeneticStructure;
import yeswecan.utils.ArrayPrinter;

/**
 *
 * @author cmonit1
 */
public class TrueCodonSiteClassTable extends CodonSiteClass {
    
    protected List<String> rows;
    protected String header;
    protected int numSiteClasses;
    protected GeneticStructure genStruct;
    protected int[][] siteClasses; // i dimension is nt site, j dimension is frame
    
    public TrueCodonSiteClassTable(int[][] siteClasses, int numSiteClasses, GeneticStructure genStruct){
        // assign fields
        
        this.numSiteClasses = numSiteClasses;
        this.genStruct = genStruct;
        this.siteClasses = siteClasses;
    }
    
    public String getHeader(){
        List<String> row = new ArrayList<String>();
        row.add(Constants.HEADER); // blank to match the regular expression at the start of each row
        row.add(Constants.GENE);
        row.add(Constants.CODON_SITE);
        for (int iSiteClass = 0; iSiteClass < this.numSiteClasses; iSiteClass++) {
            row.add(Integer.toString(iSiteClass));
        }
        row.add(Constants.NT_SITES);
//        row.add(Constants.CODON);
//        row.add(Constants.AA);
        return String.join(Constants.DEL, row);
    }
    
    public List<String> getGeneCodonTrueRows(int gene){
        int[] geneNucSites = getGeneContiguousSites(gene, this.genStruct);
        int[][] codons = getCodons(geneNucSites);
        
        List<String> allRows = new ArrayList<String>();

        for (int iCodon = 0; iCodon < codons.length; iCodon++) {
            
            if (!sitesSequential(codons[iCodon])) {
                throw new RuntimeException("Non-sequential nt sites within codon in simulation"); // this would suggest there's a bug
            }
            
            int partition = this.genStruct.getPartitionIndex(codons[iCodon][0]);// this will vary across codons, depending on layout of other genes
            int frame = this.genStruct.getFrame(partition, gene);// should be constant for all gene's codons, because we do not permit discontiguous genes 
            
            // get site classes for each site within the codon
            int[] codonSiteClasses = new int[3];
            for (int i = 0; i < codonSiteClasses.length; i++) {
                codonSiteClasses[i] = this.siteClasses[ codons[iCodon][i] ][frame];
            }
            // check that they are all the same. Would suggest a bug if not
            if (!statesEqual(codonSiteClasses)) {
                throw new RuntimeException("Nucleotide sites within a codon have differing site classes");
            }
            
            int codonSiteClass = codonSiteClasses[0]; // we can use any of the three since we now know they are the same
            
            // make row
            List<String> row = new ArrayList<String>();
            row.add(Constants.CLASSES);
            row.add(Integer.toString(gene));
            row.add(Integer.toString(iCodon+1)); // correct for zero based
            
            for (int iSiteClass = 0; iSiteClass < this.numSiteClasses; iSiteClass++) {
                String siteClassString;
                if (iSiteClass == codonSiteClass) {
                    row.add(Constants.TRUE_STRING);
                }else{
                    row.add(Constants.FALSE_STRING);
                }
            }// for iSiteClass
            
            // show the 3 nt site indices for this codon 
            int[] codonSitesPlus1 = new int[codons[iCodon].length];
            for (int i = 0; i < codonSitesPlus1.length; i++) {
                codonSitesPlus1[i] = codons[iCodon][i]+1;
            }
            row.add( ArrayPrinter.toString(codonSitesPlus1, Constants.ARGS_DELIMITER));
            
        }// for iCodon
        
        return allRows;
    }
    
    protected static boolean statesEqual(int[] states){
        for (int i = 1; i < states.length; i++) {
            if (states[i] != states[i-1]) {
                return false;
            }
        }
        return true;
    }
    
    public static void main(String[] args){
        int[] states = new int[]{ 1,2,3 };
        boolean equal = statesEqual(states);
        System.out.println(ArrayPrinter.toString(states, ",")+"\t"+equal);
        
    }
    
}
