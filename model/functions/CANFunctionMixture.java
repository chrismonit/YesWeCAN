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
import yeswecan.model.matrices.CodonAwareMatrix;
import yeswecan.model.LikelihoodCalculator;
import yeswecan.model.ProbMatrixFactory;
import yeswecan.model.ProbMatrixGenerator;
import yeswecan.model.ratioscaling.ProportionScaler;
import yeswecan.model.submodels.CANModelMixture;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.GeneticStructure;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class CANFunctionMixture implements MultivariateFunction {

    private AdvancedAlignment alignment;
    private Tree tree;
    private GeneticStructure genStruct;
    private CANModelMixture canModel;
    
    private int numSiteClasses;
  
    public CANFunctionMixture(AdvancedAlignment alignment, Tree tree, GeneticStructure genStruct, CANModelMixture can, int numSiteClasses){
        this.alignment = alignment;
        this.tree = tree;
        this.genStruct = genStruct;
        
        this.canModel = can;
        this.numSiteClasses = numSiteClasses;
        // NB 0th omega is fixed to 1.0 for neutral evolution
        
    }
    
 
    @Override
    public double value(double[] point) {
        //System.out.println("value");
        // NB 0th omega is fixed to 1.0 for neutral evolution

        Mapper.setOptimisable(this.canModel.getParameters(), point);
        

        
        double totalLogL = 0.0; // lnL for all sites
        
        for (int iSite = 0; iSite < this.alignment.getLength(); iSite++) {
            //System.out.println("\nsite "+iSite);
            
            // determine which process is active at this site
            int siteType = iSite % 3;
            int[] genes = genStruct.getGenes(iSite); // the genes present in the three frames in this partition
            
            // compute likelihood for this site
            
            double siteL = 0.0;
            // iterate over site classes, for each gene
            for (int iSiteClassA = 0; iSiteClassA < this.numSiteClasses; iSiteClassA++) {
                
//                Omega aOmega = this.canModel.getOmega(genes[0], iSiteClassA);
//                double pA = this.canModel.getProbability(genes[0], iSiteClassA); 
                for (int iSiteClassB = 0; iSiteClassB < this.numSiteClasses; iSiteClassB++) {
 //                   Omega bOmega = this.canModel.getOmega(genes[1], iSiteClassB);
 //                   double pB = this.canModel.getProbability(genes[1], iSiteClassB);
                    
                    for (int iSiteClassC = 0; iSiteClassC < this.numSiteClasses; iSiteClassC++) {
                        
                        //get the right omegas from the list of parameters
                        Omega aOmega = this.canModel.getOmega(genes[0], iSiteClassA);
                        Omega bOmega = this.canModel.getOmega(genes[1], iSiteClassB);
                        Omega cOmega = this.canModel.getOmega(genes[2], iSiteClassC);
                        
                        CodonAwareMatrix Q = new CodonAwareMatrix(
                            this.canModel.getKappa(),
                            this.canModel.getPi(),
                            new ProportionScaler(),
                            siteType,
                            aOmega, bOmega, cOmega,
                            this.canModel.getScaling()
                        );
                        
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

                        double contrib = pA * pB* pC * LikelihoodCalculator.calculateSiteLikelihood(alignment, tree, iSite, P, 1.0);

                        //System.out.println(contrib);
                        siteL += contrib;
                    } //iSiteClassC
                } // iSiteClassB
            }// iSiteClassA
            //System.out.println("siteL: " + siteL);
            totalLogL += Math.log(siteL);
            
        }// for iSite
        
        return totalLogL;
        
    } // value
    

}// class
