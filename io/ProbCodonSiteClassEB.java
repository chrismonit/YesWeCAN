/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.io;

import java.util.ArrayList;
import java.util.List;
import pal.datatype.CodonTable;
import yeswecan.Constants;
import yeswecan.model.empiricalbayes.BaseCodonEBCalculator;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.GeneticStructure;
import yeswecan.utils.ArrayPrinter;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author cmonit1
 */
public class ProbCodonSiteClassEB extends CodonSiteClass {
    
    protected AdvancedAlignment alignment;
    protected GeneticStructure genStruct;
    protected CodonTable codonTable;
    protected int numSiteClasses;
    protected boolean roundNEBValues;
    protected int representativeSequence = Constants.DISPLAY_SEQUENCE_INDEX;
    
    protected BaseCodonEBCalculator codonEB;
    
    public ProbCodonSiteClassEB(
            BaseCodonEBCalculator codonEB,
            GeneticStructure genStruct,
            int numSiteClasses,
            AdvancedAlignment alignment, 
            CodonTable codonTable,
            boolean roundNEBValues
    ){
        
        // necessary for computing EB
        this.codonEB = codonEB;
        this.genStruct = genStruct;
        this.numSiteClasses = numSiteClasses;
        
        this.roundNEBValues = roundNEBValues;
        
        // used for adding extra columns to EB output to help user work out which sites they're looking at
        this.alignment = alignment;
        this.codonTable = codonTable;
    }
    
    
    @Override
    protected String getHeader(){
        List<String> row = new ArrayList<String>();
        row.add(Constants.HEADER); // blank to match the regular expression at the start of each row
        row.add(Constants.GENE);
        row.add(Constants.CODON_SITE);
        for (int iSiteClass = 0; iSiteClass < this.numSiteClasses; iSiteClass++) {
            row.add(Integer.toString(iSiteClass));
        }
        row.add(Constants.NT_SITES);
        row.add(Constants.CODON);
        row.add(Constants.AA);
        return String.join(Constants.DEL, row);
    }
    
    @Override
    protected List<String> getGeneRows(){
        List<String> allGeneRows = new ArrayList<String>();
        
        for (int iGene = 1; iGene < this.genStruct.getNumberOfGenes()+1; iGene++) { // ignore noncoding gene
            allGeneRows.addAll( getGeneCodonEBRows(iGene) );
            allGeneRows.add("");// empty line between genes for readability
        }
        return allGeneRows;
    }
    
    protected List<String> getGeneCodonEBRows(int gene){ // each element in returned list will be a row to print/save
        
        int[] geneNucSites = getGeneContiguousSites(gene, this.genStruct);
        int[][] codons = getCodons(geneNucSites);
        
        List<String> allRows = new ArrayList<String>();
        
        for (int iCodon = 0; iCodon < codons.length; iCodon++) {
            
            List<String> row = new ArrayList<String>();
            row.add(Constants.CLASSES);
            row.add(Integer.toString(gene));
            row.add(Integer.toString(iCodon+1)); // correct for zero based
            
            double Z = this.codonEB.getNormalisationFactor(codons[iCodon]);
            
            for (int iSiteClass = 0; iSiteClass < this.numSiteClasses; iSiteClass++) {
                String nebString;
                
                if (!sitesSequential(codons[iCodon])) {
                    nebString = Constants.NO_DATA;
                }else{
                    double numerator = this.codonEB.getNumerator(codons[iCodon], iSiteClass);
                    double neb = numerator/Z;

                    if (this.roundNEBValues) {
                        nebString = Double.toString(MatrixPrinter.roundDouble(neb, Constants.DEC_PLACES));
                    }else{
                        nebString = Double.toString(neb);
                    }
                }// else
                
                row.add(nebString);
            }//iSiteClass
            
            int[] codonSitesPlus1 = new int[codons[iCodon].length];
            for (int i = 0; i < codonSitesPlus1.length; i++) {
                codonSitesPlus1[i] = codons[iCodon][i]+1;
            }
            
            row.add( ArrayPrinter.toString(codonSitesPlus1, Constants.ARGS_DELIMITER));
            
            // make string representation of codon in representativeSequence
            char[] bases = new char[codons[iCodon].length];
            for (int iCodonPosition = 0; iCodonPosition < bases.length; iCodonPosition++) {
                bases[iCodonPosition] = this.alignment.getData(this.representativeSequence, codons[iCodon][iCodonPosition]);
            }
            row.add( new String(bases) );
                        
            // make string representation of corresponding amino acid
            char aa = codonTable.getAminoAcidChar(bases);
            row.add( Character.toString(aa) );
            
            String joined = String.join(Constants.DEL, row);
            allRows.add( joined );
        }// iCodon
        
        return allRows;
    }

}
