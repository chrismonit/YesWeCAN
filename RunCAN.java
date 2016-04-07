/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import pal.tree.Tree;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Omega;
import yeswecan.cli.CommandArgs;
import yeswecan.model.CodonScaler;
import yeswecan.model.ProportionScaler;
import yeswecan.model.RatioScaler;
import yeswecan.model.functions.CANFunction;
import yeswecan.model.can.CANModel;
import yeswecan.model.hky.HKYModel;
import yeswecan.optim.Optimise;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class RunCAN extends RunModel {
    
    private CommandArgs comArgs;
    private AdvancedAlignment alignment;
    private Tree tree;
    
    protected GeneticStructure genStruct;
    
    public RunCAN(AdvancedAlignment alignment, Tree tree, CommandArgs input){
        
        this.comArgs = input;
        this.alignment = alignment;
        this.tree = tree;
        
        this.genStruct = new GeneticStructure(this.comArgs.aFrame(),
                                                            this.comArgs.bFrame(),
                                                            this.comArgs.cFrame(),
                                                            this.comArgs.lengths());
    }
    
    @Override
    public String[] getHeader(){
        ArrayList<String> columns = new ArrayList<String>();
        Collections.addAll(columns, "model", "lnL", "kappa", "A", "C", "G", "T", "sc", "0_w");
        for (int i = 0; i < this.comArgs.getGeneNumber(); i++) {
            columns.add(Integer.toString(i+1) + "_" +Constants.OMEGA_STRING); // +1 for zero based correction
        }
        return columns.toArray(new String[columns.size()]);
    }
    
    @Override
    public double[] getInitialValues(){ // NB first element does not contain lnL
        ArrayList<Double> values = RunModel.getParameterValues(makeCAN(this.comArgs).getParameters());
        double[] resultArray = new double[values.size()+2];
        resultArray[0] = Constants.CAN0_IDENTIFIER;
        resultArray[1] = Double.NaN;
        
        for (int i = 0; i < values.size(); i++) {
            resultArray[i+2] = values.get(i);
        }
        return resultArray;
    }
    
    
    // need to set whether these are fixed or not at this point
    public static CANModel makeCAN(CommandArgs comArgs){        

        HKYModel hky = RunHKY.makeHKY(comArgs);

        BranchScaling scaling = new BranchScaling(comArgs.scaling());
        if (comArgs.fix().contains(Constants.FIX_SCALING)) {
            scaling.setOptimisable(false);
        }
        
        List<Omega> omegas = new ArrayList<Omega>();
        
        Omega neutralOmega = new Omega(1.0); // for frames where there is no gene
        neutralOmega.setOptimisable(false); // never want this to change in optimisation
        omegas.add(neutralOmega);
        
//         positions of omegas correspond to the genes they represent (i.e. gene 1 omega is first)
        for (int i = 0; i < comArgs.omegas().length; i++) {
            double omegaValue = comArgs.omegas()[i];
            Omega w = new Omega( omegaValue );
            if ( comArgs.fix().contains(Integer.toString(i+1)) ) { // +1 because 0th omega is neutral
                w.setOptimisable(false);
            }
            omegas.add(w);
        }
        return new CANModel(hky, scaling, omegas);
        
    }
    
//    private RatioScaler getRatioScaler(CommandArgs comArgs){
//        switch (comArgs.getRatioScalingMethod()){
//            case Constants.PROPORTION_SCALER_IDENTIFIER: return new ProportionScaler();
//            case Constants.CODON_SCALER_IDENTIFIER: return new CodonScaler(new CodonFrequencies(comArgs.getCodonFrequencyPath()));
//            // any others we wish to add later
//            default: throw new RuntimeException("Unrecognised argument to ratio scaling option");      
//        }
//    }
    
    private static RatioScaler getRatioScaler(CommandArgs comArgs){
        RatioScaler scalerToReturn;
        int identifier;
        
        if (comArgs.getRatioScalingMethod() == Constants.PROPORTION_SCALER_IDENTIFIER){
            scalerToReturn = new ProportionScaler();
            identifier = Constants.PROPORTION_SCALER_IDENTIFIER;
        
        }else if (comArgs.getRatioScalingMethod() == Constants.CODON_SCALER_IDENTIFIER){
            scalerToReturn = new CodonScaler(new CodonFrequencies(comArgs.getCodonFrequencyPath()));
            identifier = Constants.CODON_SCALER_IDENTIFIER;
            
            System.out.println(Constants.CODON_FREQ_PATH + Constants.DEL + comArgs.getCodonFrequencyPath());
            
        }else{
            throw new RuntimeException("Unrecognised argument to ratio scaling option");
        }
        
        System.out.println(Constants.P_N + Constants.DEL + identifier + Constants.DEL + scalerToReturn.toString());
        return scalerToReturn;
                    
    }
    
    
    @Override
    public double[] fit(){
  
        CANModel can = makeCAN(this.comArgs);
        CANFunction optFunction = new CANFunction(this.alignment, this.tree, genStruct, can, getRatioScaler(this.comArgs));
        Optimise opt = new Optimise();
        CANModel result = (CANModel)opt.optNMS(optFunction, can);
        
        ArrayList<Double> values = RunModel.getParameterValues(result.getParameters());
        values.add(0, result.getLnL()); // prepend
        double[] resultArray = new double[values.size()+1];
        resultArray[0] = Constants.CAN0_IDENTIFIER;
        for (int i = 0; i < values.size(); i++) {
            resultArray[i+1] = values.get(i);
        }
        return resultArray;
 
    }
    
 
    @Override
    public double[] calculate(){
        CANModel can = makeCAN(this.comArgs);
        double[] optimisableParams = Mapper.getOptimisable(can.getParameters()); // map parameters to optimisation space, so FunctionHKY.value can use them
        CANFunction calculator = new CANFunction(this.alignment, this.tree, this.genStruct, can, getRatioScaler(this.comArgs));
        
        ArrayList<Double> values = RunModel.getParameterValues(can.getParameters());
        double lnL = calculator.value(optimisableParams);
        values.add(0, lnL);
        double[] resultArray = new double[values.size()+1];
        resultArray[0] = Constants.CAN0_IDENTIFIER;
        
        for (int i = 0; i < values.size(); i++) {
            resultArray[i+1] = values.get(i);
        }
        return resultArray;
    }
    
    
}// class
