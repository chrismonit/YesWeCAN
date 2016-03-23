/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import java.util.Random;
import pal.alignment.Alignment;
import pal.tree.ReadTree;
import pal.tree.Tree;
import swmutsel.model.parameters.BaseFrequencies;
import yeswecan.cli.CommandArgs;
import yeswecan.model.can.CANModel;
import yeswecan.model.canmix.CANModelMixture;
import yeswecan.model.hky.HKYModel;
import yeswecan.model.parameters.TsTvRatioAdvanced;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.FastaWriter;
import yeswecan.phylo.GeneticStructure;
import yeswecan.sim.SimCAN;
import yeswecan.sim.SimCANMixture;
import yeswecan.sim.SimFreqs;
import yeswecan.sim.SimHKY;
import yeswecan.sim.SimModel;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class Simulate {
    public static void main(String[] args){
        new Simulate(args);
    }
    
    private AdvancedAlignment alignment;
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
        
        // need to print a report of the params used
        // and say exactly what model used
        
        
        Alignment result = null;
        
        Random rand = new Random();
        //RunModel run;
        if (this.comArgs.getModel() == Constants.HKY_IDENTIFIER){
            // print input params using RunHKY.getInititalValues
            HKYModel hky = new HKYModel(
                    new TsTvRatioAdvanced(this.comArgs.kappa()), 
                    new BaseFrequencies(this.comArgs.pi())
            );
            // NB if -l has more than one length given, hky only uses the first in the list
            SimHKY simHky = new SimHKY(this.tree, rand, hky, 
                    this.comArgs.lengths()[0], this.comArgs.verbose());
            result = simHky.simulate();

        }
        else if (this.comArgs.getModel() == Constants.CAN0_IDENTIFIER){
            CANModel can = RunCAN.makeCAN(this.comArgs);
            SimCAN simCan = new SimCAN(
                    this.tree, rand, can, makeGenStruct(this.comArgs), this.comArgs.verbose()
            );
            result = simCan.simulate();
        }
        else if (this.comArgs.getModel() == Constants.M1_IDENTIFIER || this.comArgs.getModel() == Constants.M2_IDENTIFIER ){
            CANModelMixture canMix = RunCANMixture.makeMixture(this.comArgs, this.comArgs.getModel());
            SimCANMixture simMix = new SimCANMixture(
                    this.tree, rand, canMix, makeGenStruct(this.comArgs), this.comArgs.verbose()
            );
            result = simMix.simulate(); 
        }
        else if (this.comArgs.getModel() == Constants.CODON_FREQ_IDENTIFIER){
            SimFreqs simFreqs = new SimFreqs(this.tree, rand, makeGenStruct(this.comArgs), this.comArgs);
            result = simFreqs.simulate();
        }
        else{
            throw new RuntimeException(Constants.ERROR_PREFIX + "Invalid model argument (-m)");
        }
        
        new FastaWriter().writeFasta(result, this.comArgs.alignment());        
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
