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
    public double getNormalisationFactor( 
            ProbMatrixGenerator[][][][][] pMatGens,
            int site, int[] genes, int partition, int siteType ){
        System.out.println("Start Z calc");
        double Z = 0.0;
        for (int iSiteClassA = 0; iSiteClassA < this.numSiteClasses; iSiteClassA++) {
            for (int iSiteClassB = 0; iSiteClassB < this.numSiteClasses; iSiteClassB++) {
                for (int iSiteClassC = 0; iSiteClassC < this.numSiteClasses; iSiteClassC++) {

                    double pA = this.canModel.getProbability(genes[0], iSiteClassA); 
                    double pB = this.canModel.getProbability(genes[1], iSiteClassB);
                    double pC = this.canModel.getProbability(genes[2], iSiteClassC);
                    
                    ProbMatrixGenerator P = pMatGens[partition][iSiteClassA][iSiteClassB][iSiteClassC][siteType];
                    double L = LikelihoodCalculator.calculateSiteLikelihood(this.alignment, this.tree, site, P, 1.0);
                    double contrib = pA * pB* pC * L;
                    Z += contrib;
                    System.out.println("partition "+partition+" aFrameClass "+iSiteClassA+" bFrameClass "+iSiteClassB+" cFrameClass "+iSiteClassC+" siteType "+siteType+" L "+L+" pA "+pA+" pB "+pB+" pC "+pC);
                }// C
            }// B
            System.out.println("");
        }// A
        System.out.println("end Z calc");
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
            ArrayPrinter.print(genes, ",");
            System.out.println("\niSite "+iSite+"================================================\n");
            double Z = getNormalisationFactor(pMatGens, iSite, genes, partition, siteType); // should this be inside the gene loop?            
            // compute numerator
            
            for (int iGeneX = 0; iGeneX < this.genStruct.getNumberOfGenes()+1; iGeneX++) { // + because 0 is noncoding gene and getNumberOfGenes() only counts the coding genes
                //System.out.println("iGene "+iGeneX);
                if (this.genStruct.containsGene(partition, iGeneX)){
                    
                    int[] otherGenes = otherIntegers(genes, iGeneX);
                    int otherGeneY = otherGenes[0];
                    int otherGeneZ = otherGenes[1];
                    System.out.println("gene "+iGeneX+" others "+ArrayPrinter.toString(otherGenes, ","));
                    
                    int iGeneXFrame = this.genStruct.getFrame(partition, iGeneX);
                    int otherGeneYFrame = this.genStruct.getFrame(partition, otherGeneY);
                    int otherGeneZFrame = this.genStruct.getFrame(partition, otherGeneZ);
                    
                    System.out.println("frame(x)="+iGeneXFrame+" frame(Y)="+otherGeneYFrame+" frame(z)="+otherGeneZFrame);
                    
                    int[] geneFramesXYZ = new int[]{ iGeneXFrame, otherGeneYFrame, otherGeneZFrame };
                    
                    double sumX = 0.0; // for testing
                    System.out.println("Start iSiteClassX loop");
                    for (int iSiteClassX = 0; iSiteClassX < this.numSiteClasses; iSiteClassX++) {
                        System.out.println("iSiteClassX "+iSiteClassX);
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
                                //System.out.println("iSiteClassX="+iSiteClassX+" iSiteClassY="+iSiteClassY+" iSiteClassZ="+iSiteClassZ);
                                //System.out.println("");
                                //double t = 0.0001;
                                //MatrixPrinter.PrintMatrix(P.getP(t).getData(), "t="+t);
                                
                                double p_gene_classY = this.canModel.getProbability(otherGeneY, iSiteClassY);
                                double p_gene_classZ = this.canModel.getProbability(otherGeneZ, iSiteClassZ);
                                
                                double L = LikelihoodCalculator.calculateSiteLikelihood(this.alignment, this.tree, iSite, P, 1.0);
                                double contrib = p_gene_classY * p_gene_classZ * L;
                                System.out.println("iGene="+iGeneX+" p="+partition+" a="+aFrameClass+" b="+bFrameClass+" c="+cFrameClass+" t="+siteType+" L "+L+" px "+p_gene_classX+" py "+p_gene_classY+" pz "+p_gene_classZ);

                                sum += contrib;
                            }// Z
                        }// Y
                        double numerator = (p_gene_classX * sum);
                        
                        System.out.println("numerator sum "+sum);
                        System.out.println("Z             "+Z);
                        
                        sumX += numerator;
                        
                        probValues[iSite][iGeneX][iSiteClassX] = numerator / Z;
                        if (numerator > Z+Constants.EPSILON) {
                            
                            System.out.println("\nnumerator "+numerator);
                            System.out.println("Z "+Z);
                            System.out.println("ratio "+(numerator/Z));
                            System.out.println("");
                            System.exit(1);
                        }//if 

                    }// for iSiteClassX
                    //System.out.println("Start iSiteClassX loop");

                    System.out.println("sumX "+sumX);
                    System.out.println("   Z "+Z);
                    System.out.println("");
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
