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
import yeswecan.model.matrices.NewCodonAwareMatrix;
import yeswecan.model.submodels.CANModelSum;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.GeneticStructure;

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
    }
    

    @Override
    public double value(double[] point) {
        //System.out.println("value");
        // NB 0th omega is fixed to 1.0 for neutral evolution
        
        Mapper.setOptimisable(this.canModelSum.getParameters(), point);
        
        //double branchScaling = this.canModel.getScaling().get();
        
        double lnL = 0.0;
        
        for (int iSite = 0; iSite < this.alignment.getLength(); iSite++) {
            // determine which process is active 
            int siteType = iSite % 3;
            
            int[] genes = genStruct.getGenes(iSite); // the genes present in the three frames in this partition
            
            //get the right omegas from the list of parameters

            Omega aOmega = this.canModelSum.getOmegas().get(genes[0]);
            Omega bOmega = this.canModelSum.getOmegas().get(genes[1]);
            Omega cOmega = this.canModelSum.getOmegas().get(genes[2]);
            
            // make rate matrix, make p matrix, compute lnL for site
            
            NewCodonAwareMatrix Q = new NewCodonAwareMatrix(
                    this.canModelSum.getKappa(),
                    siteType,
                    aOmega, bOmega, cOmega,
                    this.canModelSum.getScaling()
            );
            
            //MatrixPrinter.PrintMatrix(Q.getData(), "Q", "");
            //System.out.println("w1: " + this.canModelSum.getOmegas().get(1).toString());
            
            
            
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
            //System.out.println("site_"+iSite + "\t" + sitelnL);
            
            lnL += Math.log(siteL);
            
        }// for iSite
        return lnL;
        
    }
    
    
}
