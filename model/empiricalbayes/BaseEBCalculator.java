/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.empiricalbayes;

import pal.datatype.CodonTable;
import yeswecan.model.functions.CANFunctionFreqProductsMix;
import yeswecan.model.likelihood.ProbMatrixFactory;
import yeswecan.model.likelihood.ProbMatrixGenerator;
import yeswecan.model.matrices.CANMatrixFreqProducts;
import yeswecan.model.submodels.CANModelFrequenciesMix;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;

/**
 *
 * @author cmonit1
 */
public abstract class BaseEBCalculator {
    
    
    
    public static int[] otherIntegers(int[] allInts, int intToExclude ){
        int[] remainingInts = new int[allInts.length-1];
        int adjustment = 0;
        for (int i = 0; i < allInts.length; i++) {
            if (allInts[i] == intToExclude){
                adjustment = 1;
                continue;
            }else{
                remainingInts[i-adjustment] = allInts[i];
            }
        }
        return remainingInts;
    }
    
    
    public static CANMatrixFreqProducts[][][][][] getQMatrices( 
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
    
    
    public static ProbMatrixGenerator[][][][][] createProbMatrixGenerators(CANMatrixFreqProducts[][][][][] Q_matrices){
        int nPartitions = Q_matrices.length;
        int nFrameAClasses = Q_matrices[0].length; 
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
    
}
