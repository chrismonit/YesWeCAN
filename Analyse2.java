/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import pal.alignment.AlignmentReaders;
import pal.alignment.SimpleAlignment;
import pal.datatype.Nucleotides;
import pal.tree.ReadTree;
import pal.tree.Tree;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.Parameter;
import swmutsel.model.parameters.TsTvRatio;
import yeswecan.cli.CommandArgs;
import yeswecan.model.CANFunction;
import yeswecan.model.CANModel;
import yeswecan.model.FunctionHKY;
import yeswecan.model.MutationModel;
import yeswecan.model.SubstitutionModel;
import yeswecan.optim.Optimise;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.GeneticStructure;
import yeswecan.phylo.ReorderFrequencies;
import yeswecan.phylo.States;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class Analyse2 {
    
    public static void main(String[] args) {
        
        new Analyse2(args);

    }// main
    
    private AdvancedAlignment alignment;
    private Tree tree;
    private CommandArgs comArgs;
    
    private GeneticStructure genStruct;
    
    public Analyse2(String[] args){
        this.comArgs = new CommandArgs();
        JCommander jcom = new JCommander(this.comArgs);
        
        try{
            jcom.parse(args);
        }
        catch(ParameterException ex){
            System.out.println(ex.getMessage());
        }
        
        loadData(this.comArgs.alignment(), this.comArgs.tree());
 
        this.genStruct = new GeneticStructure(this.comArgs.aFrame(),
                                                            this.comArgs.bFrame(),
                                                            this.comArgs.cFrame(),
                                                            this.comArgs.lengths());
        
        
        // calculate lnL with fixed params
        if (Boolean.parseBoolean(this.comArgs.fix())){
            calculateFixed(
                    canParameters(),
                    this.tree,
                    this.alignment
            ); 
        }
        else{
            fitCAN();
        }
        
    }
    
    //TODO make more sophistcated exceptions to help user find problems. Separate tree and alignment reading in 
    // TODO make this able to read either fasta or phylip alignments
    public void loadData(String alignmentPath, String treePath){
        try{
            this.alignment = new AdvancedAlignment(
                                new SimpleAlignment(
                                        AlignmentReaders.readFastaSequences(new FileReader(alignmentPath), new Nucleotides())));
            this.tree = new ReadTree(treePath);
        }
        catch(Exception e){
            System.out.println(Constants.ERROR_PREFIX + "Unable to load alignment or tree file(s)");
            e.printStackTrace();
            System.exit(1);
        }

    }//loadData
    
    
    // need to set whether these are fixed or not at this point
    public List<Parameter> canParameters(){
        ArrayList<Parameter> parameters = new ArrayList<Parameter>();
        
        TsTvRatio kappa = new TsTvRatio(this.comArgs.kappa());
        parameters.add(kappa);
      
        double[] frequencies = new double[States.NT_STATES]; // will be in correct order, whatever that may be
        
        if (Boolean.parseBoolean(this.comArgs.tcag())){
            frequencies = ReorderFrequencies.pamlToAlpha(this.comArgs.pi());
        }
        else{
            frequencies = this.comArgs.pi();
        }

        BaseFrequencies pi = new BaseFrequencies(frequencies);
        parameters.add(pi);
        
        BranchScaling scaling = new BranchScaling(this.comArgs.scaling());
        parameters.add(scaling);
        
        Omega neutral = new Omega(1.0); // for frames where there is no gene
        neutral.setOptimisable(false); // don't want this to change in optimisation
        parameters.add(neutral);
        
        // if a given omega is meant to be fixed, can do so here
        for (double w : this.comArgs.omegas()) {
            parameters.add(new Omega(w));
        }
        
        return parameters;
    }
    
    
    public void calculateFixed(List<Parameter> parameters, Tree tree, AdvancedAlignment alignment){
        
        //NB this might get confusing with some parameters being fixed and others not...
        double[] optimisableParams = Mapper.getOptimisable(parameters); // map parameters to optimisation space, so FunctionHKY.value can use them
        
        CANFunction calculator = new CANFunction(alignment, tree, genStruct, new CANModel(parameters));
        double lnL = calculator.value(optimisableParams);
        System.out.println("lnL: " + lnL + " "); // better to have it print the input parameters too, so you can see input and output together
    }
    

    
    
    // start the optimisation
    public void fitCAN(){
        System.out.println("hello");
        List<Parameter> parameters = canParameters();
        for (Parameter p : parameters) {
            System.out.println(p.toString());
        }
        CANFunction optFunction = new CANFunction(this.alignment, this.tree, genStruct, new CANModel(parameters));
        Optimise opt = new Optimise();
        SubstitutionModel result = opt.optNMS(optFunction, new CANModel(canParameters()));
        
        System.out.println("opt lnL: "+result.getLnL());
        for (Parameter p : result.getParameters()) {
            System.out.println(p.toString());
        }
    }
    
    
}
