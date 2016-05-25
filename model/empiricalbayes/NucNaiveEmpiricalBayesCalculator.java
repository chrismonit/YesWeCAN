/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.empiricalbayes;


import pal.datatype.CodonTable;
import pal.tree.Tree;
import yeswecan.Constants;
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
 */
public class NucNaiveEmpiricalBayesCalculator extends EmpiricalBayesCalculator {
    
    protected double[][][] probValues;
    
    protected AdvancedAlignment alignment;
    protected Tree tree;
    protected GeneticStructure genStruct;
    protected CANModelFrequenciesMix canModel;
    
    protected CodonFrequencies[] codonFrequenciesArray;
    protected CodonTable codonTable;
    
    protected int numSiteClasses;
    
    public NucNaiveEmpiricalBayesCalculator(
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

        this.probValues = computeProbValues();
        
    }
    
    
    
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
                    double L = LikelihoodCalculator.calculateSiteLikelihood(this.alignment, this.tree, site, P, 1.0);
                    double contrib = pA * pB* pC * L;
                    Z += contrib;
                }// C
            }// B
        }// A
        
        return Z;
    }

    protected double[][][] computeProbValues(){
        
        // make Q matrices, based on MLEs
        CANMatrixFreqProducts[][][][][] Q_matrices = getQMatrices(
                this.genStruct, this.canModel, this.codonFrequenciesArray, 
                this.codonTable, this.numSiteClasses);
        
        // make P matrix generators, from these MLE Qs
        ProbMatrixGenerator[][][][][] pMatGens = createProbMatrixGenerators(Q_matrices);
        
        double[][][] probValues = new double[this.genStruct.getTotalLength()][this.genStruct.getNumberOfGenes()+1][this.numSiteClasses]; // +1 to include noncoding gene
        
        for (int iSite = 0; iSite < this.genStruct.getTotalLength(); iSite++) {
            
            int partition = this.genStruct.getPartitionIndex(iSite);
            int siteType = iSite % 3;
            int[] genes = genStruct.getGenes(iSite); // the genes present in the three frames in this partition
            
            double Z = getNormalisationFactor(pMatGens, iSite, genes, partition, siteType); // prob of observing data at this site, marginalising over all gene site classes       
            
            // compute numerator for each gene which is present at this site
            for (int iGeneX = 0; iGeneX < this.genStruct.getNumberOfGenes()+1; iGeneX++) { // + because 0 is noncoding gene and getNumberOfGenes() only counts the coding genes
                
                if (this.genStruct.containsGene(partition, iGeneX)){
                    
                    int[] otherGenes = otherIntegers(genes, iGeneX); // the genes present in the other two frames
                    int otherGeneY = otherGenes[0];
                    int otherGeneZ = otherGenes[1];
                    
                    // determine frames of the three genes present
                    int iGeneXFrame = this.genStruct.getFrame(partition, iGeneX);
                    int otherGeneYFrame = this.genStruct.getFrame(partition, otherGeneY);
                    int otherGeneZFrame = this.genStruct.getFrame(partition, otherGeneZ);
                    int[] geneFramesXYZ = new int[]{ iGeneXFrame, otherGeneYFrame, otherGeneZFrame };
                    
                    double sumX = 0.0; // for error checking. Want to make sure marginalising over gene X equals Z
                    
                    // iterating over each of the site classes for gene X
                    for (int iSiteClassX = 0; iSiteClassX < this.numSiteClasses; iSiteClassX++) {
                        double p_gene_classX = this.canModel.getProbability(iGeneX, iSiteClassX); // p_{siteclassX}^{geneX}
                        
                        double sum = 0.0;
                        
                        for (int iSiteClassY = 0; iSiteClassY < this.numSiteClasses; iSiteClassY++) { //otherGenes[0]
                            for (int iSiteClassZ = 0; iSiteClassZ < this.numSiteClasses; iSiteClassZ++) { //otherGenes[1]
                                                                
                                int[] siteClassFrameOrdered = new int[3];
                                
                                /* need to match accessors for the site class dimensions in the pMatGens array
                                with the siteclass indexers in the for loops
                                e.g. if iGeneX is in frame A, aFrameClass=iSiteClassX
                                    if geneZ is in frame B, bFrameClass=iSiteClassZ
                                and if geneY is in frameC, cFrameClass=iSiteClassY
                                */
                                
                                siteClassFrameOrdered[ geneFramesXYZ[0] ] = iSiteClassX; // Notice that frames for genes X, Y and Z are ordered the same way in geneFramesXYZ, above
                                siteClassFrameOrdered[ geneFramesXYZ[1] ] = iSiteClassY;
                                siteClassFrameOrdered[ geneFramesXYZ[2] ] = iSiteClassZ;
                                                                
                                int aFrameClass = siteClassFrameOrdered[0];
                                int bFrameClass = siteClassFrameOrdered[1];
                                int cFrameClass = siteClassFrameOrdered[2];
                                                                                                
                                ProbMatrixGenerator P = pMatGens[partition][aFrameClass][bFrameClass][cFrameClass][siteType];

                                
                                double p_gene_classY = this.canModel.getProbability(otherGeneY, iSiteClassY);
                                double p_gene_classZ = this.canModel.getProbability(otherGeneZ, iSiteClassZ);
                                
                                double L = LikelihoodCalculator.calculateSiteLikelihood(this.alignment, this.tree, iSite, P, 1.0);
                                double contrib = p_gene_classY * p_gene_classZ * L;

                                sum += contrib;
                            }// Z
                        }// Y
                        double nebNumerator = (p_gene_classX * sum);
                                                
                        sumX += nebNumerator;
                        probValues[iSite][iGeneX][iSiteClassX] = nebNumerator / Z;
                        
                        // sanity check
                        if (nebNumerator > Z+Constants.EPSILON) {
                            System.out.println("Error computing NEB. (nebNumerator > Z+Constants.EPSILON). Exiting.");
                            System.out.println("\n NEB numerator "+nebNumerator);
                            System.out.println("Z                "+Z);
                            System.out.println("ratio "+(nebNumerator/Z));
                            System.exit(1);
                        }//if 

                    }// for iSiteClassX

                    if (sumX < Z-Constants.EPSILON || sumX > Z+Constants.EPSILON) {
                        System.out.println("Error computing NEB. (sumX < Z-Constants.EPSILON || sumX > Z+Constants.EPSILON)  Exiting.");
                        System.out.println("sumX "+sumX);
                        System.out.println("   Z "+Z);
                        System.exit(1);
                    }
                    
                }// if gene present at this site
                else{
                    for (int iSiteClass = 0; iSiteClass < this.numSiteClasses; iSiteClass++) {
                        probValues[iSite][iGeneX][iSiteClass] = Constants.NO_GENE_VALUE; // gene is not present so there's no NEB value to give
                    }
                }
                
            }// iGene
        }// iSite
        
        return probValues;
    }//computeProbValues
    

    @Override
    public double[][][] getEBValues(){
        return this.probValues;
    }

}
