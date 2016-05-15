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
    
    
//    protected double getNumerator(
//            ProbMatrixGenerator[][][][][] pMatGens,
//            int site, int[] genes, int partition, int siteType
//    ){
//        
//    }
    
    protected static int[] orderIndexers(int[] genes, int gene, int[] indexers){
        //int[] otherGenes = otherIntegers(genes, gene); // index0=Y, index1=Z
        int[] orderedIndexers = new int[genes.length];
        
        if (genes[0] == gene) {
            orderedIndexers[0] = indexers[0];
            orderedIndexers[1] = indexers[1];
            orderedIndexers[2] = indexers[2];
        }else if (genes[1] == gene){
            orderedIndexers[0] = indexers[1];
            orderedIndexers[1] = indexers[0];
            orderedIndexers[2] = indexers[2];
            
        }else if (genes[2] == gene){
            orderedIndexers[0] = indexers[1];
            orderedIndexers[1] = indexers[2];
            orderedIndexers[2] = indexers[0];
        }
        return orderedIndexers;
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
            
            for (int iGene = 0; iGene < this.genStruct.getNumberOfGenes(); iGene++) {
                
                if (this.genStruct.containsGene(partition, iGene)){
                    
                    int[] otherGenes = otherIntegers(genes, iGene);
                    
                    
                    //ArrayPrinter.print(otherGenes, ",");
                    
                    for (int iSiteClassX = 0; iSiteClassX < this.numSiteClasses; iSiteClassX++) {
                        
                        double p_gene_classX = this.canModel.getProbability(iGene, iSiteClassX); // p_{siteclassX}^{gene}
                        
                        double sum = 0.0;
                        
                        for (int iSiteClassY = 0; iSiteClassY < this.numSiteClasses; iSiteClassY++) { //otherGenes[0]
                            for (int iSiteClassZ = 0; iSiteClassZ < this.numSiteClasses; iSiteClassZ++) { //otherGenes[1]
                                
                                int[] indexers = new int[]{ iSiteClassX, iSiteClassY, iSiteClassZ };
                                int[] orderedIndexersByFrame = orderIndexers(genes, iGene, indexers);
                                
                                // iSiteClass index should be for whichever FRAME iGene is in
                                
                                int aFrameClass = orderedIndexersByFrame[0];
                                int bFrameClass = orderedIndexersByFrame[1];
                                int cFrameClass = orderedIndexersByFrame[2];
                                
                                /* need to match these with the siteclass indexers above
                                e.g. if iGene is in frame A, aFrameClass=iSiteClass
                                                                
                                */
                                
                                ProbMatrixGenerator P = pMatGens[partition][aFrameClass][bFrameClass][cFrameClass][siteType];
                                
                                double p_gene_classY = this.canModel.getProbability(otherGenes[0], iSiteClassY);
                                double p_gene_classZ = this.canModel.getProbability(otherGenes[1], iSiteClassZ);
                                
                                double contrib = p_gene_classY * p_gene_classZ * LikelihoodCalculator.calculateSiteLikelihood(this.alignment, this.tree, iSite, P, 1.0);
                                sum += contrib;
                            }// Z
                        }// Y
                        
                        probValues[iSite][iGene][iSiteClassX] = 3.1416;//(p_gene_class * sum) / Z;
                        
                    }// for iSiteClass
                    
                    
                }// if gene present at this site
                else{
                    
                    for (int iSiteClass = 0; iSiteClass < this.numSiteClasses; iSiteClass++) {
                        probValues[iSite][iGene][iSiteClass] = 2.7183;  // Constants.NO_GENE_VALUE; // gene is not present so there's no NEB value to give
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
        
        int[] result = orderIndexers(genes, gene, indexers);
        ArrayPrinter.print(result, ",");
        
    }
    

}
