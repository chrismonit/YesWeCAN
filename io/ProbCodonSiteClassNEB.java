/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.io;

import java.util.ArrayList;
import java.util.List;
import pal.datatype.CodonTable;
import pal.tree.Tree;
import yeswecan.Constants;
import yeswecan.model.empiricalbayes.CodonEmpiricalBayesCalculator;
import yeswecan.model.empiricalbayes.CodonNEBCalculator;
import yeswecan.model.empiricalbayes.EmpiricalBayesCalculator;
import yeswecan.model.likelihood.ProbMatrixGenerator;
import yeswecan.model.matrices.CANMatrixFreqProducts;
import yeswecan.model.submodels.CANModelFrequenciesMix;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;
import yeswecan.utils.ArrayPrinter;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author cmonit1
 */
public class ProbCodonSiteClassNEB extends CodonSiteClass {
    
    protected AdvancedAlignment alignment;
    protected Tree tree;
    protected GeneticStructure genStruct;
    protected CANModelFrequenciesMix canModel;
    
    protected CodonFrequencies[] codonFrequenciesArray;
    protected CodonTable codonTable;
    
    protected int numSiteClasses;
    
    protected CodonNEBCalculator codonNEB;
    
    //protected ProbMatrixGenerator[][][][][] pMatGenes;
    protected boolean roundNEBValues;
    protected int representativeSequence = Constants.DISPLAY_SEQUENCE_INDEX;
    
    protected CodonEmpiricalBayesCalculator codonEB;
    
    public ProbCodonSiteClassNEB(
            AdvancedAlignment alignment, Tree tree, 
            GeneticStructure genStruct, CANModelFrequenciesMix canModel,
            CodonFrequencies[] codonFrequenciesArray, CodonTable codonTable,
            int numSiteClasses, boolean roundNEBValues
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
                new CodonNEBCalculator(
                        this.alignment, this.tree, this.genStruct, 
                        this.canModel, this.numSiteClasses);
        
        CANMatrixFreqProducts[][][][][] Q_matrices = 
                EmpiricalBayesCalculator.getQMatrices(
                        this.genStruct, this.canModel, this.codonFrequenciesArray, 
                        this.codonTable, this.numSiteClasses);
        
        //this.pMatGenes = EmpiricalBayesCalculator.createProbMatrixGenerators(Q_matrices);
        this.roundNEBValues = roundNEBValues;
        
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
            allGeneRows.addAll( getGeneCodonNEBRows(iGene) );
            allGeneRows.add("");// empty line between genes for readability
        }
        return allGeneRows;
    }
    
    protected List<String> getGeneCodonNEBRows(int gene){ // each element in returned list will be a row to print/save
        
        int[] geneNucSites = getGeneContiguousSites(gene, this.genStruct);
        int[][] codons = getCodons(geneNucSites);
        
        List<String> allRows = new ArrayList<String>();
        
        for (int iCodon = 0; iCodon < codons.length; iCodon++) {
            
            List<String> row = new ArrayList<String>();
            row.add(Constants.CLASSES);
            row.add(Integer.toString(gene));
            row.add(Integer.toString(iCodon+1)); // correct for zero based
            
            double Z = this.codonNEB.getNormalisationFactor(codons[iCodon], this.pMatGens);
            
            for (int iSiteClass = 0; iSiteClass < this.numSiteClasses; iSiteClass++) {
                String nebString;
                
                if (!sitesSequential(codons[iCodon])) {
                    nebString = Constants.NO_DATA;
                }else{
                    double numerator = this.codonNEB.getNumerator(
                            codons[iCodon], iSiteClass, this.pMatGens);
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
