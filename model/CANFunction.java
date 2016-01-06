/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model;

import java.util.ArrayList;
import org.apache.commons.math3.analysis.MultivariateFunction;
import pal.tree.Tree;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.Parameter;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.GeneticStructure;
import yeswecan.utils.ArrayPrinter;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class CANFunction implements MultivariateFunction {

    private AdvancedAlignment alignment;
    private Tree tree;
    private GeneticStructure genStruct;
    private CANModel canModel;
    
  
    public CANFunction(AdvancedAlignment alignment, Tree tree, GeneticStructure genStruct){
        this.alignment = alignment;
        this.tree = tree;
        this.genStruct = genStruct;
        
        ArrayList<Omega> omegas = new ArrayList<Omega>();
        Omega neutral = new Omega(1.0);
        neutral.setOptimisable(false);
        omegas.add(neutral);
        
        for (int i = 0; i < this.genStruct.getNumberOfGenes(); i++) {
            omegas.add(new Omega(-1.0));
        }
        
        this.canModel = new CANModel(
                new TsTvRatioAdvanced(-1.0),
                new BaseFrequencies(), // default
                new BranchScaling(-1.0),
                omegas
        );
        
        
        
        // NB 0th omega is fixed to 1.0 for neutral evolution
        
        /* NB while this instance of CANmodel is defined here with parameter
        (the default) values, they are never used. The only way to get anything
        out from this class is through the value() method, which always populates
        the CANModel instance anew 
        */

    }
    
    
    
    @Override
    public double value(double[] point) {
        
        // NB 0th omega is fixed to 1.0 for neutral evolution
        //ArrayPrinter.print(point, ",");
        
        Mapper.setOptimisable(this.canModel.getParameters(), point);
        double branchScaling = this.canModel.getScaling().get();
        
        //System.out.println("pi: " + this.canModel.getPi().toString());
      
        
        double lnL = 0.0;
        
        for (int iSite = 0; iSite < this.alignment.getLength(); iSite++) {
            // determine which process is active 
            int siteType = iSite % 3;
            
            int[] genes = genStruct.getGenes(iSite); // the genes present in the three frames in this partition
            
            //get the right omegas from the list of parameters

            Omega aOmega = this.canModel.getOmegas().get(genes[0]);
            Omega bOmega = this.canModel.getOmegas().get(genes[1]);
            Omega cOmega = this.canModel.getOmegas().get(genes[2]);
            
            
            // make rate matrix, make p matrix, compute lnL for site
            
            CodonAwareMatrix Q = new CodonAwareMatrix(
                    this.canModel.getKappa(),
                    this.canModel.getPi(),
                    new ProportionScaler(),
                    siteType,
                    aOmega, bOmega, cOmega
            );
            
            //System.out.println("w1: " + this.canModel.getOmegas().get(1).toString());
            
            ProbMatrixGenerator P;
            try{
                P = ProbMatrixFactory.getPGenerator(Q);
            }
            catch(Exception e){
                System.out.println("\t\tFAILURE: " + iSite);
                
                for (Parameter p : this.canModel.getParameters()){
                    System.out.println(p.toString());
                }
                //MatrixPrinter.PrintMatrix(Q.getData(), "Q");
                //P = null;
                //e.printStackTrace();
                //System.exit(1);
                return Double.NEGATIVE_INFINITY;
            }
            
            
            double sitelnL = LogLikelihoodCalculator.calculateSiteLogLikelihood(this.alignment, this.tree, iSite, P, branchScaling);
            //System.out.println("site_"+iSite + "\t" + sitelnL);
            
            lnL += sitelnL;
            
        }// for iSite
        
        return lnL;
        
    }
    
    
    
    
    
}
