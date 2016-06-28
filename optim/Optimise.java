/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.optim;
//import java.util.List;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimplePointChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import swmutsel.model.parameters.Mapper;
//import swmutsel.model.parameters.Parameter;
import yeswecan.Constants;
import yeswecan.model.submodels.SubstitutionModel;
import yeswecan.utils.ArrayPrinter;



/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class Optimise {

    private ObjectiveFunction myOF;
    private SimplePointChecker myChecker;
    private MaxEval maxeval;
    private MaxIter maxiter;       

    
     public Optimise(){ //default constructor
        this(Constants.RELATIVE_TRESHOLD, 
                Constants.ABSOLUTE_THRESHOLD, 
                Constants.MAXIMUM_EVALUATIONS, 
                Constants.MAXIMUM_ITERATIONS);
    }
	
    public Optimise(int relativeThreshold, double absoluteThreashold, int maxEval, int maxIter){ // constructor allowing you to specifiy this stuff
        this.myChecker = new SimplePointChecker(relativeThreshold, absoluteThreashold); // these numbers influence the optimisation process in some way. I don't understand how. I use the numbers Asif uses in swmutsel
        this.maxeval = new MaxEval(maxEval);
        this.maxiter = new MaxIter(maxIter);
    }

    // optimise multivariate function using Nelder-Mead simplex algorithm
    // may want to call this multiple times to repeat otpimisation with different start points, perhaps by peturbing initial params each time. Asif recommends this.
    public SubstitutionModel optNMS(MultivariateFunction function, SubstitutionModel model){
        
        // 'function' is my own class which implements the MultivariateFunction interface
        
        // values in initialModel are in parameter space. Map to optimisation space:
        double[] initialValues = Mapper.getOptimisable(model.getParameters());// might not be happy about model being declared as SubstitutionModel, which is abstract
        int numFreeParameters = initialValues.length; // optimiser needs to know how many variables its dealing with
        
        SimplexOptimizer mySO = new SimplexOptimizer(myChecker);
        NelderMeadSimplex myNMS = new NelderMeadSimplex(numFreeParameters);
        
        myOF = new ObjectiveFunction(function); // ObjectiveFunction is effectively a wrapper. I don't know why it is necessary, but it is.
        
        InitialGuess guess = new InitialGuess(initialValues); // starting values
        
        // let's go:
        PointValuePair result = mySO.optimize( maxeval, maxiter, myNMS, myOF, guess );
        
//        System.out.println("");
//        System.out.println("Optimise. MLE opt space:");
//        ArrayPrinter.print(result.getPoint(), ",");
        
        Mapper.setOptimisable(model.getParameters(), result.getPoint()); // map values from optimisation space back to parameter space
        model.setLnL(result.getValue()); 
        
        return model; // model is effectively a datastructure, containing parameter values and (by the end) the lnL for that combination of parameters
    }//opt
    
    
    
    public SubstitutionModel optBOBYQA(MultivariateFunction function, SubstitutionModel model){
        double[] initialValues = Mapper.getOptimisable(model.getParameters());// might not be happy about model being declared as SubstitutionModel, which is abstract
        int numFreeParameters = initialValues.length; // optimiser needs to know how many variables its dealing with
        
        // values Richard uses:
        int numberOfInterpolationPoints = 2*numFreeParameters;
        double initialTrustRegionRadius = 0.01;
        double stoppingTrustRegionRadius = 1.0E-6;
        
        BOBYQAOptimizer bobyqa = new BOBYQAOptimizer(numberOfInterpolationPoints, initialTrustRegionRadius, stoppingTrustRegionRadius);

        myOF = new ObjectiveFunction(function);

        InitialGuess guess = new InitialGuess(initialValues); // starting values
        
        // I don't know if all of these arguments are necessary
        PointValuePair result = bobyqa.optimize( maxeval, maxiter, myOF, guess, GoalType.MAXIMIZE );

        Mapper.setOptimisable(model.getParameters(), result.getPoint()); // map values from optimisation space back to parameter space
        model.setLnL(result.getValue()); 

        return model; 

    }

  
}//class

