/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.empiricalbayes;


import pal.datatype.CodonTable;
import pal.tree.Tree;
import yeswecan.model.functions.CANFunctionFreqProductsMix;
import yeswecan.model.likelihood.ProbMatrixFactory;
import yeswecan.model.likelihood.ProbMatrixGenerator;
import yeswecan.model.matrices.CANMatrixFreqProducts;
import yeswecan.model.submodels.CANModelFrequenciesMix;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author cmonit1
 */
public class NaiveEmpiricalBayesCalculator extends EmpiricalBayesCalculator {
    
    protected double[][][] probValues;
    
    protected AdvancedAlignment alignment;
    protected Tree tree;
    protected GeneticStructure genStruct;
    protected CANModelFrequenciesMix canModel;
    
    protected CodonFrequencies[] codonFrequenciesArray;
    protected CodonTable codonTable;
    
    protected static int NUM_SITE_TYPES = 3;
    protected int numSiteClasses;
    
    public NaiveEmpiricalBayesCalculator(
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
        
        // initialise Q matrices
        
        CANMatrixFreqProducts[][][][][] Q_matrices = getQMatrices(
                        genStruct, canModel, codonFrequenciesArray, 
                        codonTable, numSiteClasses);
        
        
        ProbMatrixGenerator[][][][][] pMatGens = createProbMatrixGenerators(Q_matrices);
        
        
        //call method which is responsible for computing probValues
        // then this can be overriden is derived classes
    }
    
    // untested
    protected static CANMatrixFreqProducts[][][][][] getQMatrices( 
            GeneticStructure genStruct, CANModelFrequenciesMix canModel, CodonFrequencies[] codonFrequencies,
            CodonTable codonTable, int numSiteClasses){
        
        CANMatrixFreqProducts[][][][][] Q_matrices = CANFunctionFreqProductsMix.createUnscaledMatrices(
                        genStruct, canModel, codonFrequencies, 
                        codonTable, numSiteClasses);
        
        double nu = CANFunctionFreqProductsMix.computeNu(genStruct, Q_matrices, canModel, numSiteClasses);
        
        CANFunctionFreqProductsMix.scaleMatrices(Q_matrices, genStruct.getNumberOfPartitions(), numSiteClasses, nu);
        
        CANFunctionFreqProductsMix.scaleMatrices(Q_matrices, genStruct.getNumberOfPartitions(), numSiteClasses, canModel.getScaling().get());
        
        return Q_matrices;
    }
    
    // untested
    protected static ProbMatrixGenerator[][][][][] createProbMatrixGenerators(CANMatrixFreqProducts[][][][][] Q_matrices){
        int nPartitions = Q_matrices.length;
        int nFrameAClasses = Q_matrices[0].length; // not sure this is right!
        int nFrameBClasses = Q_matrices[0][0].length;
        int nFrameCClasses = Q_matrices[0][0][0].length;
        int nSiteTypes = Q_matrices[0][0][0][0].length;
        
        ProbMatrixGenerator[][][][][] pMatGens = new ProbMatrixGenerator[nPartitions][nFrameAClasses][nFrameBClasses][nFrameCClasses][nSiteTypes];
    
        for (int iPartition = 0; iPartition < nPartitions; iPartition++) {
                
                for (int iSiteClassA = 0; iSiteClassA < nFrameAClasses; iSiteClassA++) {
                    for (int iSiteClassB = 0; iSiteClassB < nFrameBClasses; iSiteClassB++) {
                        for (int iSiteClassC = 0; iSiteClassC < nFrameCClasses; iSiteClassC++) {
                            
                            for (int iSiteType = 0; iSiteType < nSiteTypes; iSiteType++) {
                            
                                pMatGens[iPartition][iSiteClassA][iSiteClassB][iSiteClassC][iSiteType] = 
                                        ProbMatrixFactory.getPGenerator( Q_matrices[iPartition][iSiteClassA][iSiteClassB][iSiteClassC][iSiteType] );
                            }// iSiteType
                    }// iSiteClassC
                }// iSiteClass B
                
            }// iSiteClassA
        }// iPartition
    
        return pMatGens;
    }
        
    // L
    
    
    // Z
    
    // numerator
    
    
    
    @Override
    public double[][][] getEBValues(){
        return this.probValues;
    }
    
    
    
}
