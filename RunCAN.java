/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan;

import java.util.ArrayList;
import java.util.List;
import pal.tree.Tree;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.Parameter;
import yeswecan.cli.CommandArgs;
import yeswecan.model.CANFunction;
import yeswecan.model.CANModel;
import yeswecan.model.SubstitutionModel;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.optim.Optimise;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.GeneticStructure;
import yeswecan.phylo.ReorderFrequencies;
import yeswecan.phylo.States;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class RunCAN {
    
    private CommandArgs comArgs;
    private AdvancedAlignment alignment;
    private Tree tree;
    
    private GeneticStructure genStruct;
    
    public RunCAN(AdvancedAlignment alignment, Tree tree, CommandArgs input){
        this.alignment = alignment;
        this.tree = tree;
        this. comArgs = input;
        
        this.genStruct = new GeneticStructure(this.comArgs.aFrame(),
                                                            this.comArgs.bFrame(),
                                                            this.comArgs.cFrame(),
                                                            this.comArgs.lengths());
    }
    
        // need to set whether these are fixed or not at this point
    public CANModel makeCAN(){        

        TsTvRatioAdvanced kappa = new TsTvRatioAdvanced(this.comArgs.kappa());
        if (this.comArgs.fix().contains(Constants.FIX_KAPPA)) {
            kappa.setOptimisable(false);
        }
      
        double[] frequencies = new double[States.NT_STATES]; // will be in correct order, whatever that may be
        
        if (Boolean.parseBoolean(this.comArgs.tcag())){
            frequencies = ReorderFrequencies.pamlToAlpha(this.comArgs.pi());
        }
        else{
            frequencies = this.comArgs.pi();
        }

        BaseFrequencies pi = new BaseFrequencies(frequencies);
        
        if (this.comArgs.fix().contains(Constants.FIX_FREQUENCIES)) {
            pi.setOptimisable(false);
        }
        
        BranchScaling scaling = new BranchScaling(this.comArgs.scaling());
        if (this.comArgs.fix().contains(Constants.FIX_SCALING)) {
            scaling.setOptimisable(false);
        }
        
                
        List<Omega> omegas = new ArrayList<Omega>();
        
        Omega neutralOmega = new Omega(1.0); // for frames where there is no gene
        neutralOmega.setOptimisable(false); // never want this to change in optimisation
        omegas.add(neutralOmega);
        
//         positions of omegas correspond to the genes they represent (i.e. gene 1 omega is first)
        for (int i = 0; i < this.comArgs.omegas().length; i++) {
            double omegaValue = this.comArgs.omegas()[i];
            Omega w = new Omega( omegaValue );
            if ( this.comArgs.fix().contains(Integer.toString(i+1)) ) { // +1 because 0th omega is neutral
                w.setOptimisable(false);
            }
            omegas.add(w);
        }

        return new CANModel(kappa, pi, scaling, omegas);
        
        
    }
    
    
    
        // start the optimisation
    public void fitCAN(){
        // need code which produces report about this analysis
        // i.e. the initial params, which are being fixed, the genetic structure, name of input aln and tree
        
        System.out.println("Codon Aware Nucleotide model. Optimising...\n");

        CANModel can = makeCAN(); // only one instance of CANModel is ever created. First populated with initial parameter values and, by end of optimisation, populated with MLEs
        CANFunction optFunction = new CANFunction(this.alignment, this.tree, genStruct, can);
        Optimise opt = new Optimise();
        SubstitutionModel result = opt.optNMS(optFunction, can);
        
        String delim = "\t";
        
        StringBuilder header = new StringBuilder("Header"+delim);
        StringBuilder mles = new StringBuilder("MLEs"+delim);
        
        
        for (Parameter p : result.getParameters()) {
            
            String fixedOrFree = "";
            if (p.isOptimisable()) {
                fixedOrFree = "[Free]";
            }else{
                fixedOrFree = "[Fixed]";
            }
            
            header.append(p.getArgument()+fixedOrFree+delim);
            mles.append(p.toString().replaceAll("[^0-9E,.-]", "")+delim);
            //mles.append(p.toString()+delim);
        }
        
        System.out.println("Tree"+delim+this.comArgs.tree());
        System.out.println("Aln"+delim+this.comArgs.alignment());
        System.out.println(header.toString());
        System.out.println(mles.toString());
        System.out.println("lnL"+delim+result.getLnL());
    }
    
    
    
    public void calculateFixed(CANModel can, Tree tree, AdvancedAlignment alignment){
        
        System.out.println("Calculating with fixed values. Input from CLI:");
        for (Parameter p : can.getParameters()){
            System.out.println(p.toString());
        }
        System.out.println("");

        double[] optimisableParams = Mapper.getOptimisable(can.getParameters()); // map parameters to optimisation space, so FunctionHKY.value can use them
        
        CANFunction calculator = new CANFunction(alignment, tree, genStruct, can);
        double lnL = calculator.value(optimisableParams);
        System.out.println("lnL: " + lnL + " "); // better to have it print the input parameters too, so you can see input and output together
    }
    
    
}// class
