/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.run;


import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import pal.datatype.CodonTableFactory;
import pal.tree.Tree;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Omega;
import yeswecan.Constants;
import yeswecan.cli.CommandArgs;
import yeswecan.model.codonawareness.CodonSum;
import yeswecan.model.functions.CANFunctionSum;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.model.submodels.CANModelSum;
import yeswecan.optim.Optimise;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class RunCANSum extends RunModel {
    
    private CommandArgs comArgs;
    private AdvancedAlignment alignment;
    private Tree tree;
    
    private GeneticStructure genStruct;
    
    private CodonSum codonSum;
    
    public RunCANSum(AdvancedAlignment alignment, Tree tree, CommandArgs input){
        
        this.comArgs = input;
        this.alignment = alignment;
        this.tree = tree;
        
        this.genStruct = new GeneticStructure(
                this.comArgs.aFrame(),
                this.comArgs.bFrame(),
                this.comArgs.cFrame(),
                this.comArgs.lengths()
        );
        
        this.codonSum = new CodonSum(
                new CodonFrequencies(input.getCodonFrequencyPath()), 
                CodonTableFactory.createUniversalTranslator()
        );
        System.out.println(Constants.CODON_FREQ_PATH+Constants.DEL+input.getCodonFrequencyPath());
    }
    
    
    @Override
    public String[] getHeader(){
        ArrayList<String> columns = new ArrayList<String>();
        Collections.addAll(columns, "model", "lnL", "kappa", "sc", "0_w");
        for (int i = 0; i < this.comArgs.getGeneNumber(); i++) {
            columns.add(Integer.toString(i+1) + "_" +Constants.OMEGA_STRING); // +1 for zero based correction
        }
        return columns.toArray(new String[columns.size()]);
    }
    
    
    @Override
    public double[] getInitialValues(){ // NB first element does not contain lnL
        ArrayList<Double> values = RunModel.getParameterValues(makeCANSum(this.comArgs).getParameters());
        double[] resultArray = new double[values.size()+2];
        resultArray[0] = Constants.CAN_SUM_IDENTIFIER;
        resultArray[1] = Double.NaN;
        
        for (int i = 0; i < values.size(); i++) {
            resultArray[i+2] = values.get(i);
        }
        return resultArray;
    }
    
    @Override
    public double[] calculate(){
        CANModelSum canSum = makeCANSum(this.comArgs);
        double[] optimisableParams = Mapper.getOptimisable(canSum.getParameters()); // map parameters to optimisation space, so FunctionHKY.value can use them
        CANFunctionSum calculator = new CANFunctionSum(
                this.alignment, this.tree, this.genStruct, 
                canSum, this.codonSum
        );
        
        ArrayList<Double> values = RunModel.getParameterValues(canSum.getParameters());
        double lnL = calculator.value(optimisableParams);
        values.add(0, lnL);
        
        double[] resultArray = new double[values.size()+1];
        resultArray[0] = Constants.CAN_SUM_IDENTIFIER;
        
        for (int i = 0; i < values.size(); i++) {
            resultArray[i+1] = values.get(i);
        }
        return resultArray;
    }
    
    @Override 
    public double[] fit(){
  
        CANModelSum canSum = makeCANSum(this.comArgs);
        CANFunctionSum optFunction = new CANFunctionSum(
                this.alignment, this.tree, genStruct, canSum, 
                this.codonSum
        );
        Optimise opt = new Optimise();
        CANModelSum result = (CANModelSum)opt.optNMS(optFunction, canSum);
        
        ArrayList<Double> values = RunModel.getParameterValues(result.getParameters());
        values.add(0, result.getLnL()); // prepend
        double[] resultArray = new double[values.size()+1];
        resultArray[0] = Constants.CAN_SUM_IDENTIFIER;
        for (int i = 0; i < values.size(); i++) {
            resultArray[i+1] = values.get(i);
        }
        return resultArray;
 
    }
    

    
    protected static CANModelSum makeCANSum(CommandArgs comArgs){
    
        TsTvRatioAdvanced kappa = new TsTvRatioAdvanced(comArgs.kappa());
        if (comArgs.fix().contains(Constants.FIX_KAPPA)) {
            kappa.setOptimisable(false);
        }
        
        BranchScaling scaling = new BranchScaling(comArgs.scaling());
        if (comArgs.fix().contains(Constants.FIX_SCALING)){
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
        
        return new CANModelSum(kappa, scaling, omegas );
    }
    
}
