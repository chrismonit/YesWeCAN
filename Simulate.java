/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan;

import yeswecan.run.RunModel;
import yeswecan.run.RunCANMixture;
import yeswecan.run.RunCAN;
import yeswecan.run.RunHKY;
import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import java.util.Random;
import pal.alignment.Alignment;
import pal.datatype.CodonTableFactory;
import pal.tree.ReadTree;
import pal.tree.Tree;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.Parameter;
import yeswecan.io.CommandArgs;
import yeswecan.model.codonawareness.CodonSum;
import yeswecan.model.submodels.CANModel;
import yeswecan.model.submodels.CANModelMixture;
import yeswecan.model.submodels.HKYModel;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.model.submodels.CANModelFrequencies;
import yeswecan.phylo.CodonFrequencies;
import yeswecan.phylo.FastaWriter;
import yeswecan.phylo.GeneticStructure;
import yeswecan.run.RunCANSum;
import yeswecan.sim.SimCAN;
import yeswecan.sim.SimCANMixture;
import yeswecan.sim.SimCANSum;
import yeswecan.sim.SimFreqProductsMix;
import yeswecan.sim.SimFreqs;
import yeswecan.sim.SimFreqsMix;
import yeswecan.sim.SimHKY;
import yeswecan.utils.ArrayPrinter;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class Simulate {
    public static void main(String[] args){
        new Simulate(args);
    }
    
    private Tree tree;
    private CommandArgs comArgs;
    
    public Simulate(String[] args){        
        this.comArgs = new CommandArgs();
        JCommander jcom = new JCommander(this.comArgs);
        
        try{
            jcom.parse(args);
        }
        catch(ParameterException ex){
            System.out.println(ex.getMessage());
        }
        
        this.tree = loadTree(this.comArgs.tree());
        
        
        System.out.println("Simulation");
        
        System.out.println(Constants.TREE_PATH + Constants.DEL + this.comArgs.tree());
        System.out.println(Constants.ALIGN_PATH + Constants.DEL + this.comArgs.alignment());
        
        Alignment result = null;
        
        Random rand = new Random();
        //RunModel run;
        if (this.comArgs.getModel() == Constants.HKY_IDENTIFIER){
            HKYModel hky = new HKYModel(new TsTvRatioAdvanced(this.comArgs.kappa()), new BaseFrequencies(this.comArgs.pi()));
            
            // NB if -l has more than one length given, hky only uses the first in the list
            SimHKY simHky = new SimHKY(this.tree, rand, hky, 
                    this.comArgs.lengths()[0], this.comArgs.verbose());
            
            printReport(new RunHKY(null, null, this.comArgs));
            result = simHky.simulate();
            
        }
        else if (this.comArgs.getModel() == Constants.CAN0_IDENTIFIER){
            CANModel can = RunCAN.makeCAN(this.comArgs);
            SimCAN simCan = new SimCAN(
                    this.tree, rand, can, makeGenStruct(this.comArgs), this.comArgs.verbose()
            );
            
            printReport(new RunCAN(null, null, this.comArgs));
            result = simCan.simulate();
            
        }
        else if (this.comArgs.getModel() == Constants.M1_IDENTIFIER || this.comArgs.getModel() == Constants.M2_IDENTIFIER ){
            CANModelMixture canMix = RunCANMixture.makeMixture(this.comArgs, this.comArgs.getModel(), Constants.M2_IDENTIFIER);
            SimCANMixture simMix = new SimCANMixture(
                    this.tree, rand, canMix, makeGenStruct(this.comArgs), this.comArgs.verbose()
            );
            
            printReport(new RunCANMixture(null, null, this.comArgs, this.comArgs.getModel()));
            result = simMix.simulate(); 
            
        }
        else if (this.comArgs.getModel() == Constants.CODON_FREQ_IDENTIFIER){
            System.out.println(Constants.CODON_FREQ_PATH + Constants.DEL + comArgs.getCodonFrequencyPath());
            SimFreqs simFreqs = new SimFreqs(this.tree, rand, makeGenStruct(this.comArgs), this.comArgs);
            
            System.out.println(Constants.HEADER + Constants.DEL + String.join(Constants.DEL, simFreqs.getHeader()));
            System.out.println(Constants.SIMULATION + Constants.DEL + ArrayPrinter.toString(simFreqs.getSimParameters(), Constants.DEL));
            
            result = simFreqs.simulate();
            
            System.out.println(Constants.NU + Constants.DEL + simFreqs.getMeanNu()); // NB this must be called after simulate method
        }
        else if (this.comArgs.getModel() == Constants.CODON_FREQ_MIX1_IDENTIFIER || 
                this.comArgs.getModel() == Constants.CODON_FREQ_MIX2_IDENTIFIER
                ){
            
            System.out.println(Constants.CODON_FREQ_PATH + Constants.DEL + comArgs.getCodonFrequencyPath());
            
            SimFreqsMix simFreqsMix = new SimFreqsMix(this.tree, rand, makeGenStruct(this.comArgs), this.comArgs);
            
            simFreqsMix.printTrueSites();
            
            System.out.println(Constants.HEADER + Constants.DEL + String.join(Constants.DEL, simFreqsMix.getHeader()));
            System.out.println(Constants.SIMULATION + Constants.DEL + ArrayPrinter.toString(simFreqsMix.getSimParameters(), Constants.DEL));
            
            result = simFreqsMix.simulate();
            
            System.out.println(Constants.NU + Constants.DEL + simFreqsMix.getMeanNu()); // NB this must be called after simulate method
            
        }
        else if (this.comArgs.getModel() == Constants.CODON_FREQ_MIX2_CANSIM_IDENTIFIER){
            // NB this simulation always expects there to be the number of site classes associated with the alternative model
            // to simulate with the null, just set the probability for the pos sel site class to 0
            System.out.println(Constants.CODON_FREQ_PATH + Constants.DEL + comArgs.getCodonFrequencyPath());
            SimFreqProductsMix simFreqProdMix = new SimFreqProductsMix(
                    this.tree, rand, makeGenStruct(this.comArgs), this.comArgs, this.comArgs.verbose());
            
            for (Parameter p : simFreqProdMix.getModel().getParameters()){ // TODO format parameters properly
                System.out.println(p);
            }
            
            result = simFreqProdMix.simulate();
            
        }
        else if (this.comArgs.getModel() == Constants.CAN_SUM_IDENTIFIER){
            CANModelFrequencies canSum = RunCANSum.makeCANSum(this.comArgs);
            
            CodonSum codonSum = new CodonSum(
                new CodonFrequencies(this.comArgs.getCodonFrequencyPath()), 
                CodonTableFactory.createUniversalTranslator()
            );
            
            SimCANSum simSum = new SimCANSum(
                    this.tree, rand, canSum, makeGenStruct(this.comArgs), 
                    this.comArgs.verbose(), codonSum
            );
            
            printReport(new RunCANSum(null, null, this.comArgs));
            result = simSum.simulate();
        }
        else{
            throw new RuntimeException(Constants.ERROR_PREFIX + "Invalid model argument (-m)");
        }
        
        new FastaWriter().writeFasta(result, this.comArgs.alignment());
        System.out.println("Simulation complete");
    }// constructor
        
    
    private static GeneticStructure makeGenStruct(CommandArgs comArgs){ // convenience method
        GeneticStructure genStruct = new GeneticStructure(
                comArgs.aFrame(),
                comArgs.bFrame(),
                comArgs.cFrame(),
                comArgs.lengths()
        );
        System.out.println(genStruct.toString());
        return genStruct;
    }
        
    private void printReport(RunModel run){
        String[] paramsHeader = run.getHeader();
        double[] paramsValues = run.getInitialValues();
        
        System.out.println(Constants.HEADER + Constants.DEL + String.join(Constants.DEL, paramsHeader));
        System.out.println(Constants.SIMULATION + Constants.DEL + ArrayPrinter.toString(paramsValues, Constants.DEL));
        System.out.println("");
    }
    
    
    
    
    
    
    
    
    private Tree loadTree(String treePath){  // a box for the boring exception code for reading in tree
        try{
            return new ReadTree(treePath);
        }
        catch(Exception e){
            System.out.println("Failed to read in tree for Simulator");
            e.printStackTrace();
            System.exit(1);
        }
        return null;
    }
    
    
}
