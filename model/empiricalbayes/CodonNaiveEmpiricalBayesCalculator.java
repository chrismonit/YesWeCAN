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
import yeswecan.utils.ArrayPrinter;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author cmonit1
 * 
 * Terminology:
 * 4 codons overlap with codon of interest (codon V) across the 3 frames
 * Note that codons W&X and Y&Z do not overlap, but are contiguous codons in the same frame.
 * 
 * |_012_| sites in codon
   |_VVV_| V frame
   |__WWW| WX frame
   |XX___| WX frame
   |YYY__| YZ frame
   |___ZZ| YZ frame
 */
public class CodonNaiveEmpiricalBayesCalculator  {
    
    protected double[][][] probValues;
    
    protected AdvancedAlignment alignment;
    protected Tree tree;
    protected GeneticStructure genStruct;
    protected CANModelFrequenciesMix canModel;
    
    
    protected int numSiteClasses;
    protected int[][][] geneFramesByCodon;
    
    public CodonNaiveEmpiricalBayesCalculator(
            AdvancedAlignment alignment, Tree tree, 
            GeneticStructure genStruct, CANModelFrequenciesMix canModel,
            //CodonFrequencies[] codonFrequenciesArray, CodonTable codonTable,
            int numSiteClasses
    ){
        this.alignment = alignment;
        this.tree = tree;
        this.genStruct = genStruct;
        
        this.canModel = canModel;
        // NB 0th omega is fixed to 1.0 for neutral evolution
           
        
        this.numSiteClasses = numSiteClasses;
        this.geneFramesByCodon = getGeneFramesByCodon();
        
    }
    
    protected static ProbMatrixGenerator[] getCodonSiteProbMatrices(
            ProbMatrixGenerator[][][][][] pMatGens,
            int iSiteClassV, int iSiteClassW, int iSiteClassX,
            int iSiteClassY, int iSiteClassZ, 
            int codonSite0Type, int codonSite1Type, int codonSite2Type,
            int codonSite0Partition, int codonSite1Partition, int codonSite2Partition
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

        ProbMatrixGenerator P_c0 = pMatGens[codonSite0Partition][aFrameClass][bFrameClass][cFrameClass][codonSite0Type];

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

        ProbMatrixGenerator P_c1 = pMatGens[codonSite1Partition][aFrameClass][bFrameClass][cFrameClass][codonSite1Type];

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

        ProbMatrixGenerator P_c2 = pMatGens[codonSite2Partition][aFrameClass][bFrameClass][cFrameClass][codonSite2Type];
        
        return new ProbMatrixGenerator[]{ P_c0, P_c1, P_c2 };
    }
    
    
    protected static int[][][] getGeneFramesByCodon(){
        int[][][] array = new int[3][3][3];
        
        array[0][0] = new int[]{ 0,1,2 }; // codon site 0 and it is alpha site, v=frame0, xw=frame1 and yz=frame2
        array[1][0] = new int[]{ 2,0,1 };
        array[2][0] = new int[]{ 1,2,0 };

        array[0][1] = new int[]{ 1,2,0 };
        array[1][1] = new int[]{ 0,1,2 };
        array[2][1] = new int[]{ 2,0,1 };

        array[0][2] = new int[]{ 2,0,1 };
        array[1][2] = new int[]{ 1,2,0 };// codonsite 1 and it is gamma site, v=frame1, xw=frame2 and yz=frame0
        array[2][2] = new int[]{ 0,1,2 };
        return array;
    }
    
    
    // e.g. getFrame(1, 2, 2) returns the frame for codonPos=1, siteType=gamma and desired frame is yz
    public int getFrame(int codonPosition, int siteType, int codonFrame){
        return this.geneFramesByCodon[codonPosition][siteType][codonFrame];
    }
    
    public double getNormalisationFactor(int[] codonSites, ProbMatrixGenerator[][][][][] pMatGens){
        double sum = 0.0;
        for (int iSiteClassV = 0; iSiteClassV < this.numSiteClasses; iSiteClassV++) {
            sum += getNumerator(codonSites, iSiteClassV, pMatGens);
        }
        return sum;
    }
    

    
    public double getNumerator(int[] codonSites, int codonVSiteClass, ProbMatrixGenerator[][][][][] pMatGens){
        
        int codonSite0Type = codonSites[0] % 3;
        int codonSite1Type = codonSites[1] % 3;
        int codonSite2Type = codonSites[2] % 3;
        
        int codonSite0Partition = this.genStruct.getPartitionIndex(codonSites[0]);
        int codonSite1Partition = this.genStruct.getPartitionIndex(codonSites[1]);
        int codonSite2Partition = this.genStruct.getPartitionIndex(codonSites[2]);
        
        int vFrame = getFrame(0, codonSite0Type, 0 );
        int wxFrame = getFrame(1, codonSite1Type, 1 );
        int yzFrame = getFrame(2, codonSite2Type, 2 );
                
        int[] codonSite0Genes = this.genStruct.getGenes(codonSites[0]);
        int[] codonSite1Genes = this.genStruct.getGenes(codonSites[1]);
        int[] codonSite2Genes = this.genStruct.getGenes(codonSites[2]);
        
        //ArrayPrinter.print(codonSite2Genes, ",");
        
        /*
            We assume the 4 codons overlapping with codonV are complete codons
            i.e. no genes start/end WITHIN codons w, x, y or z.
            Different genes may be present in the same frame,
            e.g. if codons w and x are from different genes,
            with a partition boundary between codonSite0 and codonSite1.
            (ignoring codons at splice boundaries can be justified since 
            splice sites are conserved at nt level)
        */
        int vGene = codonSite0Genes[vFrame]; // same gene at all three codon sites, by definition
        int wGene = codonSite1Genes[wxFrame]; // the same as codonSite2Genes[wxFrame], following assumption described above
        int xGene = codonSite0Genes[wxFrame];
        int yGene = codonSite0Genes[yzFrame]; // the same as codonSite1Genes[yzFrame]
        int zGene = codonSite2Genes[yzFrame];
                
        Probabilities vProbs = this.canModel.getProbabilities().get( vGene );
        double vProbSiteClass = vProbs.get()[codonVSiteClass];
        Probabilities wProbs = this.canModel.getProbabilities().get( wGene );
        Probabilities xProbs = this.canModel.getProbabilities().get( xGene );
        Probabilities yProbs = this.canModel.getProbabilities().get( yGene );
        Probabilities zProbs = this.canModel.getProbabilities().get( zGene );
//        System.out.println("v "+vProbs);
//        System.out.println("w "+wProbs);
//        System.out.println("x "+xProbs);
//        System.out.println("y "+yProbs);
//        System.out.println("z "+zProbs);
        
        double sum = 0.0;
        
        for (int iSiteClassW = 0; iSiteClassW < this.numSiteClasses; iSiteClassW++) {
            for (int iSiteClassX = 0; iSiteClassX < this.numSiteClasses; iSiteClassX++) {
                for (int iSiteClassY = 0; iSiteClassY < this.numSiteClasses; iSiteClassY++) {
                    for (int iSiteClassZ = 0; iSiteClassZ < this.numSiteClasses; iSiteClassZ++) {

                        double probProduct = 
                                wProbs.get()[iSiteClassW] * xProbs.get()[iSiteClassX] * yProbs.get()[iSiteClassY] * zProbs.get()[iSiteClassZ];

                        double likelihoodProduct = 1.0;

                        ProbMatrixGenerator[] sitePMatGens = getCodonSiteProbMatrices(pMatGens,
                                codonVSiteClass, iSiteClassW, iSiteClassX, iSiteClassY, iSiteClassZ,
                                codonSite0Type, codonSite1Type, codonSite2Type,
                                codonSite0Partition, codonSite1Partition, codonSite2Partition);

                        ProbMatrixGenerator P_c0 = sitePMatGens[0];
                        ProbMatrixGenerator P_c1 = sitePMatGens[1];
                        ProbMatrixGenerator P_c2 = sitePMatGens[2];
                        

                        likelihoodProduct *= LikelihoodCalculator.calculateSiteLikelihood(this.alignment, this.tree, codonSites[0], P_c0, 1.0);                           
                        likelihoodProduct *= LikelihoodCalculator.calculateSiteLikelihood(this.alignment, this.tree, codonSites[1], P_c1, 1.0);
                        likelihoodProduct *= LikelihoodCalculator.calculateSiteLikelihood(this.alignment, this.tree, codonSites[2], P_c2, 1.0);

                        double contrib = probProduct * likelihoodProduct; 
                        sum += contrib;

                    }// iSiteClassZ
                }// iSiteClassY                    
            }// iSiteClassX                
        }// iSiteClassW
                
        return vProbSiteClass * sum;
    }// getNormalisationFactor
    
    

}