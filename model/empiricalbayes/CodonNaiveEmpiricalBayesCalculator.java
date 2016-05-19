/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.empiricalbayes;

import pal.datatype.CodonTable;
import pal.tree.Tree;
import swmutsel.model.parameters.Probabilities;
import yeswecan.model.likelihood.LikelihoodCalculator;
import yeswecan.model.likelihood.ProbMatrixGenerator;
import yeswecan.model.submodels.CANModelFrequenciesMix;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;

/**
 *
 * @author cmonit1
 */
public class CodonNaiveEmpiricalBayesCalculator extends EmpiricalBayesCalculator {
    
    protected double[][][] probValues;
    
    protected AdvancedAlignment alignment;
    protected Tree tree;
    protected GeneticStructure genStruct;
    protected CANModelFrequenciesMix canModel;
    
    protected CodonFrequencies[] codonFrequenciesArray;
    protected CodonTable codonTable;
    
    protected int numSiteClasses;
    
    public CodonNaiveEmpiricalBayesCalculator(
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

        //this.probValues = computeProbValues();
        
    }
    
    protected ProbMatrixGenerator[] getCodonSiteProbMatrices(
            ProbMatrixGenerator[][][][][] pMatGens,
            int partition,
            int iSiteClassV, int iSiteClassW, int iSiteClassX,
            int iSiteClassY, int iSiteClassZ, 
            int codonSite0Type, int codonSite1Type, int codonSite2Type
    ){
        int aFrameClass, bFrameClass, cFrameClass;
                            
        // for site C0
        if (codonSite0Type == 0) { // C0 is alpha site
            aFrameClass = iSiteClassV;
            bFrameClass = iSiteClassX;
            cFrameClass = iSiteClassY;
        }else if (codonSite0Type == 1){ // C0 is beta site
            aFrameClass = iSiteClassY;
            bFrameClass = iSiteClassV;
            cFrameClass = iSiteClassX;
        }else{ // (vFrame == 2). C0 is gamma site
            aFrameClass = iSiteClassX;
            bFrameClass = iSiteClassY;
            cFrameClass = iSiteClassV;
        }

        ProbMatrixGenerator P_c0 = pMatGens[partition][aFrameClass][bFrameClass][cFrameClass][codonSite0Type];

        // site C1
        if (codonSite1Type == 1) { // C1 is beta site
            aFrameClass = iSiteClassV;
            bFrameClass = iSiteClassW;
            cFrameClass = iSiteClassY;
        }else if (codonSite1Type == 2){ // C1 is gamma site
            aFrameClass = iSiteClassY;
            bFrameClass = iSiteClassV;
            cFrameClass = iSiteClassW;
        }else{ // codonSite1Type == 0. C1 is alpha site
            aFrameClass = iSiteClassW;
            bFrameClass = iSiteClassY;
            cFrameClass = iSiteClassV;
        }

        ProbMatrixGenerator P_c1 = pMatGens[partition][aFrameClass][bFrameClass][cFrameClass][codonSite1Type];

        // site C2
        if (codonSite2Type == 2) { // C2 is gamma site
            aFrameClass = iSiteClassV;
            bFrameClass = iSiteClassW;
            cFrameClass = iSiteClassZ;
        }else if (codonSite2Type == 0){ // C2 is alpha
            aFrameClass = iSiteClassZ;
            bFrameClass = iSiteClassV;
            cFrameClass = iSiteClassW;
        }else{ // codonSite2Type == 1. //C2 is beta
            aFrameClass = iSiteClassW;
            bFrameClass = iSiteClassZ;
            cFrameClass = iSiteClassV;
        }

        ProbMatrixGenerator P_c2 = pMatGens[partition][aFrameClass][bFrameClass][cFrameClass][codonSite2Type];
        
        return new ProbMatrixGenerator[]{ P_c0, P_c1, P_c2 };
    }
    
    
    protected static int[][] geneFrames = new int[][]{
        { 0, 1, 2 },
        { 1, 2, 0 },
        { 2, 0, 1 }
    };
    
    public double getNormalisationFactor(int codonSite0, ProbMatrixGenerator[][][][][] pMatGens){
        double sum = 0.0;
        for (int iSiteClassV = 0; iSiteClassV < this.numSiteClasses; iSiteClassV++) {
            sum += getNumerator(codonSite0, iSiteClassV, pMatGens);
        }
        return sum;
    }
    
    public double getNumerator(int codonSite0, int codonVSiteClass, ProbMatrixGenerator[][][][][] pMatGens){
        int codonSite0Type = codonSite0 % 3;
        int codonSite1Type = (codonSite0+1) % 3;
        int codonSite2Type = (codonSite0+2) % 3;
        
        int partition = this.genStruct.getPartitionIndex(codonSite0);
        
        int vFrame = geneFrames[codonSite0Type][0];
        int wxFrame = geneFrames[codonSite0Type][1];
        int yzFrame = geneFrames[codonSite0Type][2];
        
        // NB this is codonSite0 only, assume same genes present at other three sites
        int[] genes = this.genStruct.getGenes(codonSite0);
        
        Probabilities vProbs = this.canModel.getProbabilities().get( genes[vFrame] );
        double vProbSiteClass = vProbs.get()[codonVSiteClass];
        // NB these are the same:
        Probabilities wProbs = this.canModel.getProbabilities().get( genes[wxFrame] );
        Probabilities xProbs = this.canModel.getProbabilities().get( genes[wxFrame] );
        // NB these are the same:
        Probabilities yProbs = this.canModel.getProbabilities().get( genes[yzFrame] );
        Probabilities zProbs = this.canModel.getProbabilities().get( genes[yzFrame] );
                
        double sum = 0.0;
        
        for (int iSiteClassW = 0; iSiteClassW < this.numSiteClasses; iSiteClassW++) {
            for (int iSiteClassX = 0; iSiteClassX < this.numSiteClasses; iSiteClassX++) {
                for (int iSiteClassY = 0; iSiteClassY < this.numSiteClasses; iSiteClassY++) {
                    for (int iSiteClassZ = 0; iSiteClassZ < this.numSiteClasses; iSiteClassZ++) {

                        double probProduct = 
                                wProbs.get()[iSiteClassW] * xProbs.get()[iSiteClassX] * yProbs.get()[iSiteClassY] * zProbs.get()[iSiteClassZ];

                        double likelihoodProduct = 1.0;

                        ProbMatrixGenerator[] sitePMatGens = getCodonSiteProbMatrices(pMatGens, partition,
                                codonVSiteClass, iSiteClassW, iSiteClassX, iSiteClassY, iSiteClassZ,
                                codonSite0Type, codonSite1Type, codonSite2Type);

                        ProbMatrixGenerator P_c0 = sitePMatGens[0];
                        ProbMatrixGenerator P_c1 = sitePMatGens[1];
                        ProbMatrixGenerator P_c2 = sitePMatGens[2];

                        likelihoodProduct *= LikelihoodCalculator.calculateSiteLikelihood(this.alignment, this.tree, codonSite0, P_c0, 1.0);                           
                        likelihoodProduct *= LikelihoodCalculator.calculateSiteLikelihood(this.alignment, this.tree, codonSite0+1, P_c1, 1.0);
                        likelihoodProduct *= LikelihoodCalculator.calculateSiteLikelihood(this.alignment, this.tree, codonSite0+2, P_c2, 1.0);

                        double contrib = probProduct * likelihoodProduct; 
                        sum += contrib;

                    }// iSiteClassZ
                }// iSiteClassY                    
            }// iSiteClassX                
        }// iSiteClassW
                
        return vProbSiteClass * sum;
    }// getNormalisationFactor
    
    
    @Override
    public double[][][] getEBValues(){
        return this.probValues;
    }
}
