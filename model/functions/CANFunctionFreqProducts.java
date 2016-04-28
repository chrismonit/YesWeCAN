/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.functions;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.exception.MaxCountExceededException;
import pal.datatype.CodonTable;
import pal.tree.Tree;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.Parameter;
import yeswecan.model.likelihood.LikelihoodCalculator;
import yeswecan.model.likelihood.ProbMatrixFactory;
import yeswecan.model.likelihood.ProbMatrixGenerator;
import yeswecan.model.matrices.CANMatrixFreqProducts;
import yeswecan.model.submodels.CANModelFrequencies;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;
import yeswecan.phylo.States;

/**
 *
 * @author cmonit1
 */
public class CANFunctionFreqProducts implements MultivariateFunction {
    
    private AdvancedAlignment alignment;
    private Tree tree;
    private GeneticStructure genStruct;
    private CANModelFrequencies canModel;
    
    private CodonFrequencies[] codonFrequenciesArray;
    private CodonTable codonTable;
    
    public CANFunctionFreqProducts(
            AdvancedAlignment alignment, Tree tree, 
            GeneticStructure genStruct, CANModelFrequencies canModel,
            CodonFrequencies[] codonFrequenciesArray, CodonTable codonTable
    ){
        this.alignment = alignment;
        this.tree = tree;
        this.genStruct = genStruct;
        
        this.canModel = canModel;
        // NB 0th omega is fixed to 1.0 for neutral evolution
           
        this.codonFrequenciesArray = codonFrequenciesArray;
        this.codonTable = codonTable;
    }
    
    
    public static CANMatrixFreqProducts[][] createUnscaledMatrices(
                GeneticStructure genStruct, CANModelFrequencies canModel, 
                CodonFrequencies[] codonFrequenciesArray,
                CodonTable codonTable
        ){
        
        CANMatrixFreqProducts[][] Q_matrices = new CANMatrixFreqProducts[genStruct.getNumberOfPartitions()][3];
        
        for (int iPartition = 0; iPartition < genStruct.getNumberOfPartitions(); iPartition++) {
            for (int iSiteType = 0; iSiteType < 3; iSiteType++) {
                // create unscaled Q matrix ( Q_{0} )
                int[] genes = genStruct.getGenesByPartition(iPartition);
                Omega[] omegas = new Omega[3];
                CodonFrequencies[] geneCodonFrequenciesArray = new CodonFrequencies[3];
                
                for (int iFrame = 0; iFrame < 3; iFrame++) {
                    omegas[iFrame] = canModel.getOmegas().get(genes[iFrame]);
                    geneCodonFrequenciesArray[iFrame] = codonFrequenciesArray[genes[iFrame]];
                }
                
                Q_matrices[iPartition][iSiteType] = new CANMatrixFreqProducts(
                    canModel.getKappa(), iSiteType, omegas,
                    geneCodonFrequenciesArray, codonTable
                );
                
            }// for iSiteType
        }// for iPartition
        return Q_matrices;
    }
    
    // TODO this is copied directly from CANFunctionSum. Should refactor so they share methods
    public static double computeNu(GeneticStructure genStruct, CANMatrixFreqProducts[][] Q_matrices_unscaled){
        /* \nu = 
            frac{n_{\el}}{ sum_{z} \sum_{k} n_{z,k} r_{z,k} }
            where $z$ represents partition, $k$ site type and $n_{\el}$ is the total number of sites in alignment
        */
        double nuDenominator = 0.0;
        
        for (int iPartition = 0; iPartition < genStruct.getNumberOfPartitions(); iPartition++) {
            for (int iSiteType = 0; iSiteType < 3; iSiteType++) {
                double rate = Q_matrices_unscaled[iPartition][iSiteType].getTotalRate();
                double siteTypeCount = (double)genStruct.getSiteTypeCount(iPartition, iSiteType);
                
                nuDenominator += siteTypeCount * rate;
            } // iSiteType
        }// iPartition
        
        return (double)genStruct.getTotalLength() / nuDenominator;
    }
    
    // TODO this is copied directly from CANFunctionSum. Should refactor so they share methods
    public static void scaleMatrices(GeneticStructure genStruct, CANMatrixFreqProducts[][] Q_matrices, double scalar){
        // scales matrices 'in place', i.e. without creating new array
        
        for (int iPartition = 0; iPartition < genStruct.getNumberOfPartitions(); iPartition++) {
            for (int iSiteType = 0; iSiteType < 3; iSiteType++) {
                
                for (int i = 0; i < States.NT_STATES; i++) {
                    for (int j = 0; j < States.NT_STATES; j++) {
                        Q_matrices[iPartition][iSiteType].multiplyEntry(i, j, scalar);
                        
                    }// column j
                }// row i
                
            }// iSiteType
        }// iPartition
    }
    
    
    @Override
    public double value(double[] point) {
        //System.out.println("value");
        // NB 0th omega is fixed to 1.0 for neutral evolution
        
        Mapper.setOptimisable(this.canModel.getParameters(), point);

        CANMatrixFreqProducts[][] Q_matrices = createUnscaledMatrices(this.genStruct, 
                this.canModel, this.codonFrequenciesArray, this.codonTable);
        
        double nu = computeNu(this.genStruct, Q_matrices);

        scaleMatrices(this.genStruct, Q_matrices, nu);
        
        //MatrixPrinter.PrintMatrix(Q_matrices[0][0].getData(), "after nu scaling");
        
        scaleMatrices(this.genStruct, Q_matrices, this.canModel.getScaling().get());
        
        //MatrixPrinter.PrintMatrix(Q_matrices[0][0].getData(), "after sc scaling");

        
        double lnL = 0.0;
        
        for (int iSite = 0; iSite < this.alignment.getLength(); iSite++) {
            int partition = this.genStruct.getPartitionIndex(iSite);
            int siteType = iSite % 3;
            
            CANMatrixFreqProducts Q = Q_matrices[partition][siteType];
            
            ProbMatrixGenerator P;
            try{
                P = ProbMatrixFactory.getPGenerator(Q);
            }
            catch(MaxCountExceededException e){
                System.out.println("Eigendecomposition failure: " + iSite);
                
                for (Parameter p : this.canModel.getParameters()){
                    System.out.println(p.toString());
                }
                //MatrixPrinter.PrintMatrix(Q.getData(), "Q");
                //P = null;
                //e.printStackTrace();
                //System.exit(1);
                return Double.NEGATIVE_INFINITY;
            }
            
            double siteL = LikelihoodCalculator.calculateSiteLikelihood(this.alignment, this.tree, iSite, P, 1.0);
            
            lnL += Math.log(siteL);
            
        }// for iSite
        return lnL;
        
    }
    
    
    
    
    
    
}// class
