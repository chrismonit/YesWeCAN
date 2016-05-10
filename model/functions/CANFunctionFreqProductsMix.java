/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.functions;

import org.apache.commons.math3.analysis.MultivariateFunction;
import pal.datatype.CodonTable;
import pal.tree.Tree;
import swmutsel.model.parameters.Omega;
import yeswecan.model.matrices.CANMatrixFreqProducts;
import yeswecan.model.submodels.CANModelFrequenciesMix;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;

/**
 *
 * @author cmonit1
 */
public class CANFunctionFreqProductsMix { //implements MultivariateFunction { //TODO
    private AdvancedAlignment alignment;
    private Tree tree;
    private GeneticStructure genStruct;
    private CANModelFrequenciesMix canModel;
    
    private CodonFrequencies[] codonFrequenciesArray;
    private CodonTable codonTable;
    
    private static int NUM_SITE_TYPES = 3;
    private int numSiteClasses;
    
    
    public CANFunctionFreqProductsMix(
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
    }
    
    
    public static CANMatrixFreqProducts[][][][][] createUnscaledMatrices(
                GeneticStructure genStruct, CANModelFrequenciesMix canModel, 
                CodonFrequencies[] codonFrequenciesArray,
                CodonTable codonTable, int numSiteClasses
        ){
        
        CANMatrixFreqProducts[][][][][] Q_matrices = 
                new CANMatrixFreqProducts[genStruct.getNumberOfPartitions()][numSiteClasses][numSiteClasses][numSiteClasses][NUM_SITE_TYPES];
        
        for (int iPartition = 0; iPartition < genStruct.getNumberOfPartitions(); iPartition++) {
            
            int[] genes = genStruct.getGenesByPartition(iPartition);
            Omega[] omegas = new Omega[3];
            CodonFrequencies[] geneCodonFrequenciesArray = new CodonFrequencies[3];
            
            for (int iFrame = 0; iFrame < 3; iFrame++) {
                geneCodonFrequenciesArray[iFrame] = codonFrequenciesArray[genes[iFrame]];
            }

            for (int iSiteClassA = 0; iSiteClassA < numSiteClasses; iSiteClassA++){
                for (int iSiteClassB = 0; iSiteClassB < numSiteClasses; iSiteClassB++){
                    for (int iSiteClassC = 0; iSiteClassC < numSiteClasses; iSiteClassC++){
                    
                        omegas[0] = canModel.getGeneAndSiteClassOmega(genes[0], iSiteClassA);
                        omegas[1] = canModel.getGeneAndSiteClassOmega(genes[1], iSiteClassB);
                        omegas[2] = canModel.getGeneAndSiteClassOmega(genes[2], iSiteClassC);

                        for (int iSiteType = 0; iSiteType < 3; iSiteType++) {
                            
                            Q_matrices[iPartition][iSiteClassA][iSiteClassB][iSiteClassC][iSiteType] = new CANMatrixFreqProducts(
                                canModel.getKappa(), iSiteType, omegas,
                                geneCodonFrequenciesArray, codonTable
                            );

                        }// for iSiteType

                    }// siteclassC
                    
                }// siteClassB
            
            }// siteClassA

        }// for iPartition
        return Q_matrices;
    }
    
    
    
}
