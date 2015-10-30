
/*
May want this to be able to satisfy more than one optimisation type, e.g.
linear and gradient (will need a class which computes gradients if the latter)

Given the param values, will have to construct a RateMatrix here
Can pass alignment and tree to constructor

*/

package yeswecan.model;

import java.util.List;
import org.apache.commons.math3.analysis.MultivariateFunction;
import pal.tree.Tree;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.Parameter;
import swmutsel.model.parameters.Mapper;
import yeswecan.Constants;

import yeswecan.model.LogLikelihoodCalculator;
import yeswecan.model.RateMatrix;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.utils.MatrixPrinter;



/**
 *
 * @author cmonit1
 */
public class Function implements MultivariateFunction {
    
    

    
    private AdvancedAlignment alignment;
    private Tree tree;
    private MutationModel mutModel;
    
    public Function(AdvancedAlignment alignment, Tree tree){
        this.alignment = alignment;
        this.tree = tree;
        
        /* NB while this instance of MutationModel is defined here with parameter
        (the default) values, they are never used. The only way to get anything
        out from this class is through the value() method, which always populates
        the MutationModel instance anew 
        */
        this.mutModel = new MutationModel(new TsTvRatioAdvanced(Constants.DEFAULT_KAPPA), 
                new BaseFrequencies(Constants.DEFAULT_PI));
        
        
    }
    
    private int i = 0; // debugging
    
    public double value(double[] point){
        
        
        /* point vector arrives with values in optim space.
            need to map back to real parameter space and then calculate likelihood
        */
 
        Mapper.setOptimisable(this.mutModel.getParameters(), point);


        //System.out.print("i: " + i +"\t");
        i++;
        
        boolean verbose = false;
        
        if (verbose) {
            System.out.print(this.mutModel.getKappa().toString() + "\t");
    //        System.out.println("point[0]: " + point[0]);
    //        System.out.println("");

//            System.out.println(mutModel.getPi().toString());
//
//            System.out.println("length: " + point.length);
//
//            String pointString = "";
//                for (int i = 0; i < point.length; i++) {
//                    pointString += "," + Double.toString(point[i]);
//            }
//            System.out.println("pointString: " + pointString);
//            //MatrixPrinter.PrintMatrix(Q.getData(), "Q: " + i);
//
//            System.out.println("");
//            System.out.println("");

        }
        
        
        // make Q matrix
        RateMatrix Q = new RateMatrix(this.mutModel.getKappa(), this.mutModel.getPi());
        
        
        
        //make P matrix generator
        ProbMatrixGenerator P = ProbMatrixFactory.getPGenerator(Q);
        
        
        // can then compute likelihood

        double lnL = 0.0;
        for (int iSite=0; iSite < this.alignment.getLength(); iSite++){
            lnL += LogLikelihoodCalculator.calculateSiteLogLikelihood(this.alignment, this.tree, iSite, P);
        }
        

        //System.out.println("lnL: " + lnL);
        
        //System.out.println("");

        return lnL;
    }// value
    
}// class
