/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.optim;
import java.util.List;
import org.apache.commons.math3.analysis.MultivariateFunction;
//import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimplePointChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Parameter;
import yeswecan.Constants;



/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class Optimise {

    private ObjectiveFunction myOF;
    private SimplePointChecker myChecker;
    private MaxEval maxeval;
    private MaxIter maxiter;       

    
     public Optimise(){
        this(Constants.RELATIVE_TRESHOLD, Constants.ABSOLUTE_THRESHOLD, Constants.MAXIMUM_EVALUATIONS, Constants.MAXIMUM_ITERATIONS);
    }
	
    public Optimise(int relativeThreshold, double absoluteThreashold, int maxEval, int maxIter){
        this.myChecker = new SimplePointChecker(relativeThreshold, absoluteThreashold);
        this.maxeval = new MaxEval(maxEval);
        this.maxiter = new MaxIter(maxIter);
    }

    // optimise multivariate function using Nelder-Mead simplex algorithm
    public PointValuePair optNMS(MultivariateFunction function, List<Parameter> initial){
        
        double[] initialValues = Mapper.getOptimisable(initial);
        int numFreeParameters = initialValues.length;
        
        SimplexOptimizer mySO = new SimplexOptimizer(myChecker);
        NelderMeadSimplex myNMS = new NelderMeadSimplex(numFreeParameters);
        
        InitialGuess guess = new InitialGuess(initialValues);
        return mySO.optimize( maxeval, maxiter, myNMS, myOF, guess );

    }//opt
    
    
  
}//class

