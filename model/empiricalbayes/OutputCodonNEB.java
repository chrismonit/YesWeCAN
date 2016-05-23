/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.empiricalbayes;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import pal.datatype.CodonTable;
import pal.tree.Tree;
import yeswecan.Constants;
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
public class OutputCodonNEB {
    
    protected AdvancedAlignment alignment;
    protected Tree tree;
    protected GeneticStructure genStruct;
    protected CANModelFrequenciesMix canModel;
    
    protected CodonFrequencies[] codonFrequenciesArray;
    protected CodonTable codonTable;
    
    protected int numSiteClasses;
    
    protected CodonNaiveEmpiricalBayesCalculator codonNEB;
    
    protected ProbMatrixGenerator[][][][][] pMatGenes;
    protected boolean roundNEBValues;
    protected int representativeSequence = Constants.DISPLAY_SEQUENCE_INDEX;
    
    public OutputCodonNEB(
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
                new CodonNaiveEmpiricalBayesCalculator(
                        this.alignment, this.tree, this.genStruct, 
                        this.canModel, this.numSiteClasses);
        
        CANMatrixFreqProducts[][][][][] Q_matrices = 
                EmpiricalBayesCalculator.getQMatrices(
                        this.genStruct, this.canModel, this.codonFrequenciesArray, 
                        this.codonTable, this.numSiteClasses);
        
        this.pMatGenes = EmpiricalBayesCalculator.createProbMatrixGenerators(Q_matrices);
        this.roundNEBValues = roundNEBValues;
        
    }
    
    public void print(){
        System.out.println( getHeader() );
        for (String row : getAllGenesCodonNEBRows()){
            System.out.println(row);
        }
    }
    
    public void write(String filePath){
    
         try{
            FileWriter fileWriter = new FileWriter(filePath);
            BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
            
            bufferedWriter.write( getHeader() );
            bufferedWriter.newLine();
            for (String row : getAllGenesCodonNEBRows()){
                bufferedWriter.write(row);
                bufferedWriter.newLine();
            }

            bufferedWriter.close();
        
        }catch (IOException e){
            e.printStackTrace();
        }
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
        row.add(Constants.CODON);
        row.add(Constants.AA);
        return String.join(Constants.DEL, row);
    }
    
    public List<String> getAllGenesCodonNEBRows(){
        List<String> allGeneRows = new ArrayList<String>();
        
        for (int iGene = 1; iGene < this.genStruct.getNumberOfGenes()+1; iGene++) { // ignore noncoding gene
            allGeneRows.addAll( getGeneCodonNEBRows(iGene) );
            allGeneRows.add("");// empty line between genes for readability
        }
        return allGeneRows;
    }
    
    public List<String> getGeneCodonNEBRows(int gene){ // each element in returned list will be a row to print/save
        
        int[] geneNucSites = getGeneContiguousSites(gene);
        int[][] codons = getCodons(geneNucSites);
        
        List<String> allRows = new ArrayList<String>();
        
        for (int iCodon = 0; iCodon < codons.length; iCodon++) {
            
            List<String> row = new ArrayList<String>();
            row.add(Constants.CLASSES);
            row.add(Integer.toString(gene));
            row.add(Integer.toString(iCodon+1)); // correct for zero based
            
            double Z = this.codonNEB.getNormalisationFactor(codons[iCodon], this.pMatGenes);
            
            for (int iSiteClass = 0; iSiteClass < this.numSiteClasses; iSiteClass++) {
                String nebString;
                
                if (!sitesSequential(codons[iCodon])) {
                    nebString = Constants.NO_DATA;
                }else{
                    double numerator = this.codonNEB.getNumerator(
                    codons[iCodon], iSiteClass, this.pMatGenes);
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
    
    
    protected static int[] intArray(List<Integer> intList){
        int[] array = new int[intList.size()];
        for (int i = 0; i < intList.size(); i++) {
            array[i] = intList.get(i);
        }
        return array;
    }
    
    public static boolean sitesSequential(int[] codon){
        return (codon[1] == codon[0]+1 && codon[2] == codon[0]+2);
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
    public static int[][] getCodons(int[] sites){
        
        int nCodons = sites.length/3; // integer division. Only considers number of complete triplets
        
        int[][] codons = new int[nCodons][3];
        
        for (int iSite = 0; iSite < sites.length-2; iSite+=3) {
            int codonIndex = iSite/3;
            codons[codonIndex] = new int[]{ sites[iSite], sites[iSite+1], sites[iSite+2] };
        }
        return codons;
    }
    
    
}
