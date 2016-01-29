/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan;

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
public class RunHKY {
    private CommandArgs comArgs;
    private AdvancedAlignment alignment;
    private Tree tree;
    
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
        //return Arrays.asList(kappa, pi);        
    }
    
    
    public static void calculateFixed(HKYModel hky, Tree tree, AdvancedAlignment alignment){
        double[] optimisableParams = Mapper.getOptimisable(hky.getParameters()); // map parameters to optimisation space, so FunctionHKY.value can use them
        FunctionHKY calculator = new FunctionHKY(alignment, tree);
        double lnL = calculator.value(optimisableParams);
        System.out.println("lnL: " + lnL + " "); // better to have it print the input parameters too, so you can see input and output together
    }
    

    
    
    // start the optimisation
    public void fit(){
        FunctionHKY optFunction = new FunctionHKY(this.alignment, this.tree);
        Optimise opt = new Optimise();
        SubstitutionModel result = opt.optNMS(optFunction, makeHKY(this.comArgs));
        
        System.out.println("opt lnL: "+result.getLnL());
        System.out.println( result.getParameters().get(0).toString());
        System.out.println(result.getParameters().get(1).toString());
    }
    
    
    
}// class
