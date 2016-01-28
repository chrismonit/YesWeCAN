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
import swmutsel.model.parameters.Probabilities;
import yeswecan.cli.CommandArgs;
import yeswecan.model.CANFunction;
import yeswecan.model.CANModel;
import yeswecan.model.CANModelOriginal;
import yeswecan.model.SubstitutionModel;
import yeswecan.model.parameters.OmegaNegative;
import yeswecan.model.parameters.OmegaPositive;
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
                
        loadData(this.comArgs.alignment(), this.comArgs.tree(), Boolean.parseBoolean(this.comArgs.phy()));
 
        this.genStruct = new GeneticStructure(this.comArgs.aFrame(),
                                                            this.comArgs.bFrame(),
                                                            this.comArgs.cFrame(),
                                                            this.comArgs.lengths());
        
        
        // calculate lnL with fixed params
        if (this.comArgs.fix().contains("all")){
            calculateFixed(
                    makeCAN(),
                    this.tree,
                    this.alignment
            ); 
        }
        else{
            fitCAN();
        }
        
    }
    
    //TODO make more sophistcated exceptions to help user find problems. Separate tree and alignment reading in 
    // could do with a clever catch for when the wrong format is presented (ie incongruous with the -phy value
    // the exception thrown is an index out of range in the lnL calculator
    public void loadData(String alignmentPath, String treePath, Boolean readPhylip){
        try{
            SimpleAlignment simple;
            if (readPhylip)
                simple = new SimpleAlignment(AlignmentReaders.readPhylipClustalAlignment(new FileReader(alignmentPath), new Nucleotides()));
            else
                simple = new SimpleAlignment(AlignmentReaders.readFastaSequences(new FileReader(alignmentPath), new Nucleotides()));
            
            this.alignment = new AdvancedAlignment(simple);
                                
            this.tree = new ReadTree(treePath);
        }
        catch(Exception e){
            System.out.println(Constants.ERROR_PREFIX + "Unable to load alignment or tree file(s)");
            e.printStackTrace();
            System.exit(1);
        }

    }//loadData
    
    
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
        
        // 'omega distribution' refers to the set of parameter objects which define the omegas for a given gene
        // each gene has its own set of proabitlities and w_0, w_1 and w_2, which are bundled into a list for that gene
        
        // first doing this for M2a, then will have to do separately for M1a
        ArrayList<ArrayList<Parameter>> allOmegaDistributions = new ArrayList<ArrayList<Parameter>>();
        
        ArrayList<Parameter> neutral = new ArrayList<Parameter>(); // for no-gene frames
        Probabilities neutralProbs = new Probabilities(new double[]{0.0, 1.0, 0.0}); // all density on w_1 class, because w=1
        OmegaNegative neutralNeg = new OmegaNegative(0.0); // should never be used
        Omega neutralOmegaNeutral = new Omega(1.0);
        OmegaPositive neutralPos = new OmegaPositive(0.0); // should never be used
        
        neutralProbs.setOptimisable(false);
        neutralNeg.setOptimisable(false);
        neutralOmegaNeutral.setOptimisable(false);
        neutralPos.setOptimisable(false);
        
        neutral.add(neutralPos);
        neutral.add(neutralNeg);
        neutral.add(neutralOmegaNeutral);
        neutral.add(neutralPos);
        
        allOmegaDistributions.add(neutral); // 0th element is for neutral evolution, in gene-less frames
        // iterating over genes
        for (int i = 0; i < this.comArgs.getGeneNumber(); i++) {
            ArrayList<Parameter> omegaDistribution = new ArrayList<Parameter>();
            Probabilities geneProbs = new Probabilities(new double[]{ this.comArgs.prob0()[i], this.comArgs.prob1()[i], this.comArgs.prob2()[i]}); 
            OmegaNegative geneNeg = new OmegaNegative(this.comArgs.omega0()[i]);
            Omega geneNeutral = new Omega(1.0);
            geneNeutral.setOptimisable(false); // fix w_1 = 1 for every gene
            OmegaPositive genePos = new OmegaPositive(this.comArgs.omega2()[i]);
            
            omegaDistribution.add(geneProbs);
            omegaDistribution.add(geneNeg);
            omegaDistribution.add(geneNeutral);
            omegaDistribution.add(genePos);
            
            allOmegaDistributions.add(omegaDistribution);
        }
        
        
        
//        
//        List<Omega> omegas = new ArrayList<Omega>();
//        
//        Omega neutralOmega = new Omega(1.0); // for frames where there is no gene
//        neutralOmega.setOptimisable(false); // never want this to change in optimisation
//        omegas.add(neutralOmega);
        
////         positions of omegas correspond to the genes they represent (i.e. gene 1 omega is first)
//        for (int i = 0; i < this.comArgs.omegas().length; i++) {
//            double omegaValue = this.comArgs.omegas()[i];
//            Omega w = new Omega( omegaValue );
//            if ( this.comArgs.fix().contains(Integer.toString(i+1)) ) { // +1 because 0th omega is neutral
//                w.setOptimisable(false);
//            }
//            omegas.add(w);
//        }

        
        
        return new CANModel(kappa, pi, scaling, allOmegaDistributions);
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
    

    
    
    // start the optimisation
    public void fitCAN(){
        // need code which produces report about this analysis
        // i.e. the initial params, which are being fixed, the genetic structure, name of input aln and tree
        
        System.out.println("Codon Aware Nucleotide model. Optimising...\n");

        CANModel can = makeCAN(); // only one instance of CANModelOriginal is ever created. First populated with initial parameter values and, by end of optimisation, populated with MLEs
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
    
}
