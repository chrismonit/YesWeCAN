/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.run;

import java.util.ArrayList;
import java.util.Collections;
import pal.datatype.CodonTable;
import pal.datatype.CodonTableFactory;
import pal.tree.Tree;
import yeswecan.Constants;
import yeswecan.io.CommandArgs;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;

/**
 *
 * @author cmonit1
 */
public class RunCANFreqProductsMix extends RunModel {
    
    private CommandArgs comArgs;
    private AdvancedAlignment alignment;
    private Tree tree;
    
    private GeneticStructure genStruct;
    
    private CodonTable codonTable;
    private CodonFrequencies[] codonFrequenciesArray;
    
    private int numSiteClasses;
    
    public RunCANFreqProductsMix(AdvancedAlignment alignment, Tree tree, CommandArgs input, int model){
        this.comArgs = input;
        this.alignment = alignment;
        this.tree = tree;
        
        this.genStruct = new GeneticStructure(
                this.comArgs.aFrame(),
                this.comArgs.bFrame(),
                this.comArgs.cFrame(),
                this.comArgs.lengths()
        );
    
        this.codonFrequenciesArray = new CodonFrequencies[this.genStruct.getNumberOfGenes()+1];

        this.codonFrequenciesArray[0] = new CodonFrequencies(); // 1/64 for frames where there is no gene

        // use same codonFrequencies reference for all coding genes. Could make gene-specific if desired
        CodonFrequencies codonFrequencies = new CodonFrequencies(input.getCodonFrequencyPath());        
        System.out.println(Constants.CODON_FREQ_PATH+Constants.DEL+input.getCodonFrequencyPath());
        for (int iGene = 0; iGene < this.genStruct.getNumberOfGenes(); iGene++) {
            this.codonFrequenciesArray[iGene+1] = codonFrequencies;
        }
        
        this.codonTable = CodonTableFactory.createUniversalTranslator();
        
        if (alignment != null) { // in Simulate class, I make instances of this class but pass in null alignment and tree
            if (this.alignment.getLength() != this.genStruct.getTotalLength()){
                throw new RuntimeException("Sum of partition lengths (-l) and alignment length differ! Please check that partition lengths are correct.");
            }
        }
                
        this.numSiteClasses = numberSiteClasses(model);
        
    }// constructor
    
    
    @Override
    public String[] getHeader(){
        ArrayList<String> columns = new ArrayList<String>();
        Collections.addAll(columns, "model", "lnL", "kappa", "sc");
        for (int iGene = 0; iGene < this.comArgs.getGeneNumber(); iGene++) {
            for (int jClass = 0; jClass < this.numSiteClasses; jClass++) {
                columns.add(Integer.toString(iGene+1) + Constants.WITIHIN_FIELD_SEPARATOR + Constants.OMEGA_STRING + Integer.toString(jClass));
            }
            
            for (int jClass = 0; jClass < this.numSiteClasses; jClass++) {
                columns.add(Integer.toString(iGene+1) + Constants.WITIHIN_FIELD_SEPARATOR + Constants.PROB_STRING + Integer.toString(jClass));
            }
        }
        return columns.toArray(new String[columns.size()]);
    }
    
    //TODO
    @Override
    public double[] getInitialValues(){ // NB first element does not contain lnL
        return new double[]{0.1};
    }
    
    //TODO
    @Override
    public double[] calculate(){
        return new double[]{.2};
    }
    
    // TODO
    @Override
    public double[] fit(){
        return new double[]{0.3};
    }
    
    // TODO refactor. Identical method in SimFreqsMix
    public static int numberSiteClasses(int mixtureModel){ // based on near identical method in RunCANMixture
        int numSiteClasses = -1;

        if (mixtureModel == Constants.CODON_FREQ_MIX2_IDENTIFIER)
            numSiteClasses = Constants.NUM_M2_SITE_CLASSES;
        else
            numSiteClasses = Constants.NUM_M1_SITE_CLASSES;
        return numSiteClasses;
    }

    
}
