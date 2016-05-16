/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.empiricalbayes;


import pal.datatype.CodonTable;
import pal.tree.Tree;
import yeswecan.Constants;
import yeswecan.model.functions.CANFunctionFreqProductsMix;
import yeswecan.model.likelihood.LikelihoodCalculator;
import yeswecan.model.likelihood.ProbMatrixFactory;
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
        

        
        this.probValues = computeProbValues();
        
        
        
        //call method which is responsible for computing probValues
        // then this can be overriden is derived classes
    }
    
    
    // untested
    protected double getNormalisationFactor( 
            ProbMatrixGenerator[][][][][] pMatGens,
            int site, int[] genes, int partition, int siteType ){
        
        double Z = 0.0;
        for (int iSiteClassA = 0; iSiteClassA < this.numSiteClasses; iSiteClassA++) {
            for (int iSiteClassB = 0; iSiteClassB < this.numSiteClasses; iSiteClassB++) {
                for (int iSiteClassC = 0; iSiteClassC < this.numSiteClasses; iSiteClassC++) {

                    double pA = this.canModel.getProbability(genes[0], iSiteClassA); 
                    double pB = this.canModel.getProbability(genes[1], iSiteClassB);
                    double pC = this.canModel.getProbability(genes[2], iSiteClassC);

                    ProbMatrixGenerator P = pMatGens[partition][iSiteClassA][iSiteClassB][iSiteClassC][siteType];

                    double contrib = pA * pB* pC * LikelihoodCalculator.calculateSiteLikelihood(this.alignment, this.tree, site, P, 1.0);
                    Z += contrib;

                }// C
            }// B
        }// A
        return Z;
    }
    

    
    protected double[][][] computeProbValues(){
        CANMatrixFreqProducts[][][][][] Q_matrices = getQMatrices(
                this.genStruct, this.canModel, this.codonFrequenciesArray, 
                this.codonTable, this.numSiteClasses);
        
        ProbMatrixGenerator[][][][][] pMatGens = createProbMatrixGenerators(Q_matrices);
        
        double[][][] probValues = new double[this.genStruct.getTotalLength()][this.genStruct.getNumberOfGenes()+1][this.numSiteClasses]; // +1 to include noncoding gene
        
        for (int iSite = 0; iSite < this.genStruct.getTotalLength(); iSite++) {
            
            int partition = this.genStruct.getPartitionIndex(iSite);
            int siteType = iSite % 3;
            int[] genes = genStruct.getGenes(iSite); // the genes present in the three frames in this partition

            // compute denominator (normalisation factor)
            
            // untested
            double Z = getNormalisationFactor(pMatGens, iSite, genes, partition, siteType);
            

            // compute numerator
            
            for (int iGeneX = 0; iGeneX < this.genStruct.getNumberOfGenes(); iGeneX++) {
                System.out.println("iGene "+iGeneX);
                if (this.genStruct.containsGene(partition, iGeneX)){
                    
                    int[] otherGenes = otherIntegers(genes, iGeneX);
                    int otherGeneY = otherGenes[0];
                    int otherGeneZ = otherGenes[1];
                    
                    int iGeneXFrame = this.genStruct.getFrame(partition, iGeneX);
                    int otherGeneYFrame = this.genStruct.getFrame(partition, otherGeneY);
                    int otherGeneZFrame = this.genStruct.getFrame(partition, otherGeneZ);
                    
                    int[] geneFramesXYZ = new int[]{ iGeneXFrame, otherGeneYFrame, otherGeneZFrame };
                    
                    for (int iSiteClassX = 0; iSiteClassX < this.numSiteClasses; iSiteClassX++) {
                        
                        double p_gene_classX = this.canModel.getProbability(iGeneX, iSiteClassX); // p_{siteclassX}^{gene}
                        
                        double sum = 0.0;
                        
                        for (int iSiteClassY = 0; iSiteClassY < this.numSiteClasses; iSiteClassY++) { //otherGenes[0]
                            for (int iSiteClassZ = 0; iSiteClassZ < this.numSiteClasses; iSiteClassZ++) { //otherGenes[1]
                                                                
                                int[] siteClassFrameOrdered = new int[3];
                                
                                siteClassFrameOrdered[ geneFramesXYZ[0] ] = iSiteClassX;  
                                siteClassFrameOrdered[ geneFramesXYZ[1] ] = iSiteClassY;
                                siteClassFrameOrdered[ geneFramesXYZ[2] ] = iSiteClassZ;
                                                                
                                int aFrameClass = siteClassFrameOrdered[0];
                                int bFrameClass = siteClassFrameOrdered[1];
                                int cFrameClass = siteClassFrameOrdered[2];
                                
                                /* need to match these with the siteclass indexers above
                                e.g. if iGeneX is in frame A, aFrameClass=iSiteClassX
                                */
                                
                                ProbMatrixGenerator P = pMatGens[partition][aFrameClass][bFrameClass][cFrameClass][siteType];
                                
                                double p_gene_classY = this.canModel.getProbability(otherGeneYFrame, iSiteClassY);
                                double p_gene_classZ = this.canModel.getProbability(otherGeneZFrame, iSiteClassZ);
                                
                                double contrib = p_gene_classY * p_gene_classZ * LikelihoodCalculator.calculateSiteLikelihood(this.alignment, this.tree, iSite, P, 1.0);
                                sum += contrib;
                            }// Z
                        }// Y
                        
                        probValues[iSite][iGeneX][iSiteClassX] = (p_gene_classX * sum) / Z;
                        
                    }// for iSiteClass
                    
                    
                }// if gene present at this site
                else{
                    
                    for (int iSiteClass = 0; iSiteClass < this.numSiteClasses; iSiteClass++) {
                        probValues[iSite][iGeneX][iSiteClass] = Constants.NO_GENE_VALUE; // gene is not present so there's no NEB value to give
                    }
                }
                
            }// iGene

        }// iSite
        
        return probValues;
    }
    
    
    protected static int[] otherIntegers(int[] allInts, int intToExclude ){
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
    
    
    protected static ProbMatrixGenerator[][][][][] createProbMatrixGenerators(CANMatrixFreqProducts[][][][][] Q_matrices){
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
        
    // L
    
    
    // Z
    
    // numerator
    
    
    
    @Override
    public double[][][] getEBValues(){
        return this.probValues;
    }
    
    public static void main(String[] args){
        System.out.println("hello");
        
        int[] genes = new int[]{ 2,3,4 };
        int gene = 3;
        int[] indexers = new int[]{ 0, 1, 2 };
        

        
    }
    

}
