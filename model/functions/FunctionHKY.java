
/*
May want this to be able to satisfy more than one optimisation type, e.g.
linear and gradient (will need a class which computes gradients if the latter)

Given the param values, will have to construct a RateMatrix here
Can pass alignment and tree to constructor

*/

package yeswecan.model.functions;

import yeswecan.model.submodels.HKYModel;
import java.util.List;
import org.apache.commons.math3.analysis.MultivariateFunction;
import pal.tree.Tree;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Parameter;
import swmutsel.model.parameters.TsTvRatio;
import yeswecan.Constants;

import yeswecan.model.likelihood.LikelihoodCalculator;
import yeswecan.model.likelihood.LikelihoodCalculator;
import yeswecan.model.likelihood.ProbMatrixFactory;
import yeswecan.model.likelihood.ProbMatrixGenerator;
import yeswecan.model.submodels.HKYModel;
import yeswecan.model.matrices.RateMatrix;
import yeswecan.model.matrices.RateMatrix;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author cmonit1
 */
public class FunctionHKY implements MultivariateFunction {
    
    private AdvancedAlignment alignment;
    private Tree tree;
    private HKYModel mutModel;
    
    public FunctionHKY(AdvancedAlignment alignment, Tree tree){
        this.alignment = alignment;
        this.tree = tree;
        
        /* NB while this instance of HKYModel is defined here with parameter
        (the default) values, they are never used. The only way to get anything
        out from this class is through the value() method, which always populates
        the HKYModel instance anew 
        */
        this.mutModel = new HKYModel(new TsTvRatioAdvanced(TsTvRatioAdvanced.getDefault()), 
                new BaseFrequencies(BaseFrequencies.getDefault()));
        
        
    }
    

    public double value(double[] point){
  
        /* point vector arrives with values in optim space.
            need to map back to real parameter space and then calculate likelihood
        */
 
        Mapper.setOptimisable(this.mutModel.getParameters(), point);

        // make Q matrix
        RateMatrix Q = new RateMatrix(this.mutModel.getKappa(), this.mutModel.getPi());
        
        //make P matrix generator
        ProbMatrixGenerator P = ProbMatrixFactory.getPGenerator(Q);        
        
        // can then compute likelihood

        double lnL = 0.0;
        for (int iSite=0; iSite < this.alignment.getLength(); iSite++){
            double siteL = LikelihoodCalculator.calculateSiteLikelihood(this.alignment, this.tree, iSite, P, 1.0);
            //System.out.println("site_"+iSite + "\t" + sitelnL);
            
            lnL += Math.log(siteL);
        }

        return lnL;
    }// value
    
}// class
