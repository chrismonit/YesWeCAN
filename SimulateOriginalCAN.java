/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import pal.alignment.Alignment;
import pal.tree.ReadTree;
import pal.tree.Tree;
import yeswecan.cli.CommandArgs;
import yeswecan.phylo.FastaWriter;
import yeswecan.phylo.GeneticStructure;
import yeswecan.sim.Simulator;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class SimulateOriginalCAN {
    
    public static void main(String[] args){
        new SimulateOriginalCAN(args);
    }
    
    private CommandArgs comArgs;
    
    
    public SimulateOriginalCAN(String[] args){
        System.out.println("Simulating");
        this.comArgs = new CommandArgs();
        JCommander jcom = new JCommander(this.comArgs);
        
        try{
            jcom.parse(args);
        }
        catch(ParameterException ex){
            System.out.println(ex.getMessage());
        }
        
        Tree tree = loadTree(this.comArgs.tree());
        
        
        GeneticStructure genStruct = new GeneticStructure(this.comArgs.aFrame(),
                                                            this.comArgs.bFrame(),
                                                            this.comArgs.cFrame(),
                                                            this.comArgs.lengths());
        
        //System.out.println(genStruct.toString());
                
        double kappa = this.comArgs.kappa();
        double[] baseFrequencies = this.comArgs.pi();
        double[] omegaValues = this.comArgs.omegas();
        double scaling = this.comArgs.scaling();
        
        Simulator sim = new Simulator(tree, genStruct, kappa, baseFrequencies, omegaValues, scaling);
        Alignment result = sim.simulate();
        
        
        new FastaWriter().writeFasta(result, this.comArgs.alignment());
        System.out.println("Done");
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
