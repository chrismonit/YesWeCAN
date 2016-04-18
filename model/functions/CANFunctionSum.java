/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.model.functions;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.exception.MaxCountExceededException;
import pal.tree.Tree;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.Parameter;
import yeswecan.model.codonawareness.CodonSum;
import yeswecan.model.likelihood.LikelihoodCalculator;
import yeswecan.model.likelihood.ProbMatrixFactory;
import yeswecan.model.likelihood.ProbMatrixGenerator;
import yeswecan.model.matrices.CANMatrixSum;
import yeswecan.model.submodels.CANModelSum;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.GeneticStructure;
import yeswecan.phylo.States;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class CANFunctionSum implements MultivariateFunction {
    
    private AdvancedAlignment alignment;
    private Tree tree;
    private GeneticStructure genStruct;
    private CANModelSum canModelSum;
    private CodonSum codonSum;
    
    private CANMatrixSum[][] Q_matrices;
    
    public CANFunctionSum(
            AdvancedAlignment alignment, Tree tree, 
            GeneticStructure genStruct, CANModelSum canModelSum, 
            CodonSum codonSum){
        this.alignment = alignment;
        this.tree = tree;
        this.genStruct = genStruct;
        
        this.canModelSum = canModelSum;
        this.codonSum = codonSum;
        // NB 0th omega is fixed to 1.0 for neutral evolution
        
        this.Q_matrices = new CANMatrixSum[this.genStruct.getNumberOfPartitions()][3];
    }
    
    
    private void createUnscaledMatrices(){
        
        for (int iPartition = 0; iPartition < this.genStruct.getNumberOfPartitions(); iPartition++) {
            for (int iSiteType = 0; iSiteType < 3; iSiteType++) {
                // create unscaled Q matrix ( Q_{0} )
                int[] genes = this.genStruct.getGenesByPartition(iPartition);
                Omega aOmega = this.canModelSum.getOmegas().get(genes[0]);
                Omega bOmega = this.canModelSum.getOmegas().get(genes[1]);
                Omega cOmega = this.canModelSum.getOmegas().get(genes[2]);
                
                this.Q_matrices[iPartition][iSiteType] = new CANMatrixSum(
                    this.canModelSum.getKappa(),
                    iSiteType,
                    aOmega, bOmega, cOmega,
                    this.canModelSum.getScaling(),
                    this.codonSum
                );
                
            }// for iSiteType
        }// for iPartition
    }
    
    private static double computeNu(GeneticStructure genStruct, CANMatrixSum[][] Q_matrices_unscaled){
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
    
    private void scaleMatrices(double nu){
        for (int iPartition = 0; iPartition < this.genStruct.getNumberOfPartitions(); iPartition++) {
            for (int iSiteType = 0; iSiteType < 3; iSiteType++) {
                
                for (int i = 0; i < States.NT_STATES; i++) {
                    for (int j = 0; j < States.NT_STATES; j++) {
                        this.Q_matrices[iPartition][iSiteType].multiplyEntry(i, j, nu);
                    }// column j
                }// row i
                
            }// iSiteType
        }// iPartition
    }
    
    @Override
    public double value(double[] point) {
        //System.out.println("value");
        // NB 0th omega is fixed to 1.0 for neutral evolution
        
        Mapper.setOptimisable(this.canModelSum.getParameters(), point);

        createUnscaledMatrices();
        
        double nu = computeNu(this.genStruct, this.Q_matrices);

        scaleMatrices(nu);
        
        double lnL = 0.0;
        
        for (int iSite = 0; iSite < this.alignment.getLength(); iSite++) {
            int partition = this.genStruct.getPartitionIndex(iSite);
            int siteType = iSite % 3;
            
            CANMatrixSum Q = this.Q_matrices[partition][siteType];
            
            ProbMatrixGenerator P;
            try{
                P = ProbMatrixFactory.getPGenerator(Q);
            }
            catch(MaxCountExceededException e){
                System.out.println("Eigendecomposition failure: " + iSite);
                
                for (Parameter p : this.canModelSum.getParameters()){
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
    
    
}
