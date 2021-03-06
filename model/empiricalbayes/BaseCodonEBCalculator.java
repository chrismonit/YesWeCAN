/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.empiricalbayes;

import yeswecan.model.likelihood.ProbMatrixGenerator;

/**
 *
 * @author cmonit1
 * 
 *  * Terminology:
 * 4 codons overlap with codon of interest (codon V) across the 3 frames
 * Note that codons W&X and Y&Z do not overlap, but are contiguous codons in the same frame.
 * 
 * |_012_| sites in codon
   |_VVV_| V frame
   |__WWW| WX frame
   |XX___| WX frame
   |YYY__| YZ frame
   |___ZZ| YZ frame
 * 
 */
public abstract class BaseCodonEBCalculator extends BaseEBCalculator {
    
    public abstract double getNormalisationFactor(int[] codonSites);
    
    public abstract double getNumerator(int[] codonSites, int codonVSiteClass);

    
    
    //returns an array of ProbMatrixGenerators, one for each nuc site in the codon
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
    
    
    
}
