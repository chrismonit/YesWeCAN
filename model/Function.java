
/*
May want this to be able to satisfy more than one optimisation type, e.g.
linear and gradient (will need a class which computes gradients if the latter)

Given the param values, will have to construct a RateMatrix here
Can pass alignment and tree to constructor

*/

package yeswecan.model;

import org.apache.commons.math3.analysis.MultivariateFunction;
import pal.tree.Tree;
import yeswecan.model.RateMatrix;
import yeswecan.model.parameters.BaseFrequencies;
import yeswecan.model.parameters.TsTvRatio;
import yeswecan.model.LogLikelihoodCalculator;

import yeswecan.phylo.AdvancedAlignment;



/**
 *
 * @author cmonit1
 */
public class Function implements MultivariateFunction {
    
    
    // have constructor take mapping instance and use this to get param values from point in value method
    // lets assume for now we're just making an HKY85 and using fixed branch lengths
    // point[0] == kappa
    // point[1-3] == first 3 base sequences 
    
    /*
    * mapping instance could be passed into constructor, itself having been initialised so
    * it knows how many parameters are being optimised.
       therefore it know how long the point array should be and what values are what
    */
    
    private AdvancedAlignment alignment;
    private Tree tree;
    
    public Function(AdvancedAlignment alignment, Tree tree){
        this.alignment = alignment;
        this.tree = tree;
    }
    
    
    public double value(double[] point){
        
        //make RateMatrix
        // temporary, need to make use of some sort of mapping object
        RateMatrix Q = new RateMatrix(new TsTvRatio(point[0]), 
                new BaseFrequencies(new double[]{point[1], point[2], point[3], (1.0-(point[1]+point[2]+point[3]))}));
        
        
        double lnL = 0.0;
        for (int iSite=0; iSite < alignment.getLength(); iSite++){
            lnL += LogLikelihoodCalculator.calculateSiteLogLikelihood(alignment, tree, iSite, Q);
        }
        
        return lnL;
    }
    
}
