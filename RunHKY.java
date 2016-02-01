/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan;

import java.util.ArrayList;
import pal.tree.Tree;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.Mapper;
import yeswecan.cli.CommandArgs;
import yeswecan.model.SubstitutionModel;
import yeswecan.model.hky.FunctionHKY;
import yeswecan.model.hky.HKYModel;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.optim.Optimise;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.ReorderFrequencies;
import yeswecan.phylo.States;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class RunHKY extends RunModel {
    //using protected fields allows subclasses to access these fields
    protected CommandArgs comArgs;
    protected AdvancedAlignment alignment;
    protected Tree tree;
    
    
    
    public RunHKY(AdvancedAlignment alignment, Tree tree, CommandArgs input){
        this.alignment = alignment;
        this.tree = tree;
        this.comArgs = input;
    }
    
    public static HKYModel makeHKY(CommandArgs comArgs){
        TsTvRatioAdvanced kappa = new TsTvRatioAdvanced(comArgs.kappa());
        
        double[] frequencies = new double[States.NT_STATES]; // will be in correct order, whatever that may be
        
        if (Boolean.parseBoolean(comArgs.tcag())){
            frequencies = ReorderFrequencies.pamlToAlpha(comArgs.pi());
        }
        else{
            frequencies = comArgs.pi();
        }

        BaseFrequencies pi = new BaseFrequencies(frequencies);
        return new HKYModel(kappa, pi);
    }
    
    @Override
    public String[] getHeader(){
        return new String[]{ "kappa", "A", "C", "G", "T" };
    }
    
    @Override
    public double[] calculate(){
        HKYModel hky = makeHKY(this.comArgs);
        double[] optimisableParams = Mapper.getOptimisable(hky.getParameters()); // map parameters to optimisation space, so FunctionHKY.value can use them
        FunctionHKY calculator = new FunctionHKY(this.alignment, this.tree);
        ArrayList<Double> values = RunModel.getParameterValues(hky.getParameters());
        double lnL = calculator.value(optimisableParams);
        values.add(lnL);
        double[] resultArray = new double[values.size()];
        for (int i = 0; i < values.size(); i++) {
            resultArray[i] = values.get(i);
        }
        return resultArray;
    }
    
        
    // start the optimisation
    @Override
    public double[] fit(){
        FunctionHKY optFunction = new FunctionHKY(this.alignment, this.tree);
        Optimise opt = new Optimise();
        HKYModel result = (HKYModel)opt.optNMS(optFunction, makeHKY(this.comArgs));
        double lnL = result.getLnL();
        ArrayList<Double> values = RunModel.getParameterValues(result.getParameters());
        values.add(lnL);
        double[] resultArray = new double[values.size()];
        for (int i = 0; i < values.size(); i++) {
            resultArray[i] = values.get(i);
        }
        return resultArray;
    }
    
    
    
}// class
