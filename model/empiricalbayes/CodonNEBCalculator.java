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
import yeswecan.model.matrices.CANMatrixFreqProducts;
import yeswecan.model.submodels.CANModelFrequenciesMix;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;

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
public class CodonNEBCalculator extends CodonEmpiricalBayesCalculator {
    
    protected double[][][] probValues;
    
    protected AdvancedAlignment alignment;
    protected Tree tree;
    protected GeneticStructure genStruct;
    protected CANModelFrequenciesMix canModel;
    
    
    protected int numSiteClasses;
    protected int[][][] geneFramesByCodon;
    
    protected ProbMatrixGenerator[][][][][] pMatGens;
    
    public CodonNEBCalculator(
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
           
        
        this.numSiteClasses = numSiteClasses;
        this.geneFramesByCodon = getGeneFramesByCodon();
        
        CANMatrixFreqProducts[][][][][] Q_matrices = getQMatrices(
                this.genStruct, this.canModel, codonFrequenciesArray, 
                codonTable, this.numSiteClasses );
        
        this.pMatGens = createProbMatrixGenerators(Q_matrices);
    }
    
    
    
    // TODO could move to parent class
    // e.g. getFrame(1, 2, 2) returns the frame for codonPos=1, siteType=gamma and desired frame is yz
    public int getFrame(int codonPosition, int siteType, int codonFrame){
        return this.geneFramesByCodon[codonPosition][siteType][codonFrame];
    }
    
    @Override
    public double getNormalisationFactor(int[] codonSites){
        double sum = 0.0;
        for (int iSiteClassV = 0; iSiteClassV < this.numSiteClasses; iSiteClassV++) {
            sum += getNumerator(codonSites, iSiteClassV);
        }
        return sum;
    }
    

    @Override
    public double getNumerator(int[] codonSites, int codonVSiteClass){
        
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

                        ProbMatrixGenerator[] sitePMatGens = getCodonSiteProbMatrices(this.pMatGens,
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
