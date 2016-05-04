/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.run;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import pal.datatype.CodonTable;
import pal.datatype.CodonTableFactory;
import pal.tree.Tree;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Omega;
import yeswecan.Constants;
import yeswecan.cli.CommandArgs;
import yeswecan.model.functions.CANFunctionFreqProducts;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.model.submodels.CANModelFrequencies;
import yeswecan.optim.Optimise;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.GeneticStructure;

/**
 *
 * @author cmonit1
 */
public class RunCANFreqProducts extends RunModel {
    
    private CommandArgs comArgs;
    private AdvancedAlignment alignment;
    private Tree tree;
    
    private GeneticStructure genStruct;
    
    private CodonTable codonTable;
    private CodonFrequencies[] codonFrequenciesArray;
    
    public RunCANFreqProducts(AdvancedAlignment alignment, Tree tree, CommandArgs input){
        this.comArgs = input;
        this.alignment = alignment;
        this.tree = tree;
        
        this.genStruct = new GeneticStructure(
                this.comArgs.aFrame(),
                this.comArgs.bFrame(),
                this.comArgs.cFrame(),
                this.comArgs.lengths()
        );
    
        this.codonFrequenciesArray = new CodonFrequencies[this.genStruct.getNumberOfGenes()+1];

        this.codonFrequenciesArray[0] = new CodonFrequencies(); // 1/64 for frames where there is no gene

        // use same codonFrequencies reference for all coding genes. Could make gene-specific if desired
        CodonFrequencies codonFrequencies = new CodonFrequencies(input.getCodonFrequencyPath());        
        System.out.println(Constants.CODON_FREQ_PATH+Constants.DEL+input.getCodonFrequencyPath());
        for (int iGene = 0; iGene < this.genStruct.getNumberOfGenes(); iGene++) {
            this.codonFrequenciesArray[iGene+1] = codonFrequencies;
        }
        
        this.codonTable = CodonTableFactory.createUniversalTranslator();
        
        if (alignment != null) { // in Simulate class, I make instances of this class but pass in null alignment and tree
            if (this.alignment.getLength() != this.genStruct.getTotalLength()){
                throw new RuntimeException("Sum of partition lengths (-l) and alignment length differ! Please check that partition lengths are correct.");
            }
        }
        
        
    }// constructor
    
    
    // copied from RunCANSum. should refactor
    @Override
    public String[] getHeader(){
        ArrayList<String> columns = new ArrayList<String>();
        Collections.addAll(columns, "model", "lnL", "kappa", "sc", "0_w");
        for (int i = 0; i < this.comArgs.getGeneNumber(); i++) {
            columns.add(Integer.toString(i+1) + "_" +Constants.OMEGA_STRING); // +1 for zero based correction
        }
        return columns.toArray(new String[columns.size()]);
    }
    
    // TODO copied from RunCANSum. refactor
    @Override
    public double[] getInitialValues(){ // NB first element does not contain lnL
        ArrayList<Double> values = RunModel.getParameterValues(makeCAN(this.comArgs).getParameters());
        double[] resultArray = new double[values.size()+2];
        resultArray[0] = Constants.CAN_FREQ_IDENTIFIER;
        resultArray[1] = Double.NaN;
        
        for (int i = 0; i < values.size(); i++) {
            resultArray[i+2] = values.get(i);
        }
        return resultArray;
    }
    
    
    // TODO copied from RUNCANSum. refactor
    public static CANModelFrequencies makeCAN(CommandArgs comArgs){
    
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
        
        return new CANModelFrequencies(kappa, scaling, omegas );
    }
    
    
    @Override
    public double[] calculate(){
        CANModelFrequencies canModel = makeCAN(this.comArgs);
        double[] optimisableParams = Mapper.getOptimisable(canModel.getParameters()); // map parameters to optimisation space, so FunctionHKY.value canModel use them
        CANFunctionFreqProducts calculator = new CANFunctionFreqProducts(
                this.alignment, this.tree, this.genStruct, 
                canModel, this.codonFrequenciesArray, this.codonTable
        );
        
        ArrayList<Double> values = RunModel.getParameterValues(canModel.getParameters());
        double lnL = calculator.value(optimisableParams);
        values.add(0, lnL);
        
        double[] resultArray = new double[values.size()+1];
        resultArray[0] = Constants.CAN_FREQ_IDENTIFIER;
        
        for (int i = 0; i < values.size(); i++) {
            resultArray[i+1] = values.get(i);
        }
        return resultArray;
    }
    
    
    
    @Override 
    public double[] fit(){
  
        CANModelFrequencies canModel = makeCAN(this.comArgs);
        CANFunctionFreqProducts optFunction = new CANFunctionFreqProducts(
                this.alignment, this.tree, genStruct, canModel, 
                this.codonFrequenciesArray, this.codonTable
        );
        Optimise opt = new Optimise();
        CANModelFrequencies result = (CANModelFrequencies)opt.optNMS(optFunction, canModel);
        
        ArrayList<Double> values = RunModel.getParameterValues(result.getParameters());
        values.add(0, result.getLnL()); // prepend
        double[] resultArray = new double[values.size()+1];
        resultArray[0] = Constants.CAN_FREQ_IDENTIFIER;
        for (int i = 0; i < values.size(); i++) {
            resultArray[i+1] = values.get(i);
        }
        return resultArray;
 
    }
}
