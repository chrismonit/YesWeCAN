/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.functions;

import java.util.List;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.exception.MaxCountExceededException;
import pal.datatype.CodonTable;
import pal.tree.Tree;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.Parameter;
import swmutsel.model.parameters.Probabilities;
import static yeswecan.model.functions.CANFunctionFreqProducts.computeNu;
import static yeswecan.model.functions.CANFunctionFreqProducts.createUnscaledMatrices;
import static yeswecan.model.functions.CANFunctionFreqProducts.scaleMatrices;
import yeswecan.model.likelihood.LikelihoodCalculator;
import yeswecan.model.likelihood.ProbMatrixFactory;
import yeswecan.model.likelihood.ProbMatrixGenerator;
import yeswecan.model.matrices.CANMatrixFreqProducts;
import yeswecan.model.submodels.CANModelFrequenciesMix;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;
import yeswecan.phylo.States;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author cmonit1
 */
public class CANFunctionFreqProductsMix implements MultivariateFunction {
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
    
    
    
    public static double computeNu(
            GeneticStructure genStruct, CANMatrixFreqProducts[][][][][] Q_matrices_unscaled,
            CANModelFrequenciesMix canMix, int numSiteClasses
    ){
        
        double nuDenominator = 0.0;
        
        for (int iPartition = 0; iPartition < genStruct.getNumberOfPartitions(); iPartition++) {
            
            int[] genes = genStruct.getGenesByPartition(iPartition);
            for (int iSiteType = 0; iSiteType < 3; iSiteType++) {
                
                double rateSum = 0.0;
                for (int iSiteClassA = 0; iSiteClassA < numSiteClasses; iSiteClassA++){
                    for (int iSiteClassB = 0; iSiteClassB < numSiteClasses; iSiteClassB++){
                        for (int iSiteClassC = 0; iSiteClassC < numSiteClasses; iSiteClassC++){
                            
                            rateSum += canMix.getProbability(genes[0], iSiteClassA) *
                                       canMix.getProbability(genes[1], iSiteClassB) *
                                       canMix.getProbability(genes[2], iSiteClassC) *
                                       Q_matrices_unscaled[iPartition][iSiteClassA][iSiteClassB][iSiteClassC][iSiteType].getTotalRate();
                            //System.out.println("iPart "+iPartition+"\ttype "+iSiteType+"\tA "+iSiteClassA+"\tB "+iSiteClassB+"\tC "+iSiteClassC+"\t\t"+Q_matrices_unscaled[iPartition][iSiteClassA][iSiteClassB][iSiteClassC][iSiteType].getTotalRate());

                        }// C
                    }//  B
                }// A
                
                double numSitesPartType = (double)genStruct.getSiteTypeCount(iPartition, iSiteType); // number of sites of type iSiteType in partition iPartition
                
                nuDenominator += numSitesPartType * rateSum;
            } // iSiteType
        }// iPartition

        return (double)genStruct.getTotalLength() / nuDenominator;
    }
    
    
    public static void scaleMatrices(CANMatrixFreqProducts[][][][][] Q_matrices, int numPartitions,
            int numSiteClasses, double scalar){
        // scales matrices 'in place', i.e. without creating new array
        
        for (int iPartition = 0; iPartition < numPartitions; iPartition++) {
            for (int iSiteType = 0; iSiteType < 3; iSiteType++) {
                
                for (int iSiteClassA = 0; iSiteClassA < numSiteClasses; iSiteClassA++){
                    for (int iSiteClassB = 0; iSiteClassB < numSiteClasses; iSiteClassB++){
                        for (int iSiteClassC = 0; iSiteClassC < numSiteClasses; iSiteClassC++){
                
                            for (int i = 0; i < States.NT_STATES; i++) {
                                for (int j = 0; j < States.NT_STATES; j++) {
                                    Q_matrices[iPartition][iSiteClassA][iSiteClassB][iSiteClassC][iSiteType].multiplyEntry(i, j, scalar);

                                }// column j
                            }// row i

                        }// C
                    }// B
                }// A
                
            }// iSiteType
        }// iPartition
    }
    
    @Override 
    public double value(double[] point) {
        
        Mapper.setOptimisable(this.canModel.getParameters(), point);

        CANMatrixFreqProducts[][][][][] Q_matrices = createUnscaledMatrices(
                this.genStruct, this.canModel, this.codonFrequenciesArray, 
                this.codonTable, this.numSiteClasses);
        
        double nu = computeNu(this.genStruct, Q_matrices, this.canModel, this.numSiteClasses);
        scaleMatrices(Q_matrices, this.genStruct.getNumberOfPartitions(), this.numSiteClasses, nu);
        
        // TODO we probably don't want to bother with this scaling in future, in which case we can just remove this line
        scaleMatrices(Q_matrices, this.genStruct.getNumberOfPartitions(), this.numSiteClasses, this.canModel.getScaling().get());
        
        double totalLogL = 0.0;
        
        for (int iSite = 0; iSite < this.genStruct.getTotalLength(); iSite++) {
            int partition = this.genStruct.getPartitionIndex(iSite);
            int siteType = iSite % 3;
            int[] genes = genStruct.getGenes(iSite); // the genes present in the three frames in this partition

            double siteL = 0.0;
            // iterate over site classes, for each gene
            for (int iSiteClassA = 0; iSiteClassA < this.numSiteClasses; iSiteClassA++) {

                for (int iSiteClassB = 0; iSiteClassB < this.numSiteClasses; iSiteClassB++) {

                    for (int iSiteClassC = 0; iSiteClassC < this.numSiteClasses; iSiteClassC++) {

                        CANMatrixFreqProducts Q = Q_matrices[partition][iSiteClassA][iSiteClassB][iSiteClassC][siteType];
                        
                        ProbMatrixGenerator P;
                        try{
                            P = ProbMatrixFactory.getPGenerator(Q);
                        }
                        catch(MaxCountExceededException e){
                            System.out.println("Eigendecomposition failure: " + iSite);

                            for (Parameter p : this.canModel.getParameters()){
                                System.out.println(p.toString());
                            }

                            return Double.NEGATIVE_INFINITY;
                        }
                                                
                        // pF = p^{F}_{i} (prob for frame F's omega being that of site class i)
                        double pA = this.canModel.getProbability(genes[0], iSiteClassA); 
                        double pB = this.canModel.getProbability(genes[1], iSiteClassB);
                        double pC = this.canModel.getProbability(genes[2], iSiteClassC);
                        
                        double contrib = pA * pB* pC * LikelihoodCalculator.calculateSiteLikelihood(this.alignment, this.tree, iSite, P, 1.0);

                        siteL += contrib;
                    }// C
                }// B
            }// A
            
            totalLogL += Math.log(siteL);
        }// iSite
        
        return totalLogL;
    }// value
    
}
