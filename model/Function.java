
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
    
    
    
    public double value(double[] point){
        
        
        /* point vector arrives with values in optim space.
            need to map back to real parameter space and then calculate likelihood
        */
        
        Mapper.setOptimisable(mutModel.getParameters(), point);
        
        
        // make Q matrix
        RateMatrix Q = new RateMatrix(mutModel.getKappa(), mutModel.getPi());
        
        //TsTvRatioAdvanced k = parameters.get(0);
        //RateMatrix Q = new RateMatrix(new TsTvRatioAdvanced(k), parameters.get(1));
        
        //make P matrix generator
        
        ProbMatrixGenerator P = ProbMatrixFactory.getPGenerator(Q);
        
        
        // can then compute likelihood
        
        
        double lnL = 0.0;
        for (int iSite=0; iSite < alignment.getLength(); iSite++){
            lnL += LogLikelihoodCalculator.calculateSiteLogLikelihood(alignment, tree, iSite, P);
        }
        
        return lnL;
    }
    
}
