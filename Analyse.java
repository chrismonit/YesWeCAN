/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import java.io.FileReader;
import pal.alignment.AlignmentReaders;
import pal.alignment.SimpleAlignment;
import pal.datatype.Nucleotides;
import pal.tree.ReadTree;
import pal.tree.Tree;
import yeswecan.cli.CommandArgs;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.GeneticStructure;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class Analyse {
    public static void main(String[] args){
        new Analyse(args);
    }
    
    private AdvancedAlignment alignment;
    private Tree tree;
    private CommandArgs comArgs;
    
    private GeneticStructure genStruct;
    
    public Analyse(String[] args){
        
        // stuff that needs to hapen regardless of model being used:
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
        
        // TODO determine model to use and invoke relevant Run class
        //
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
    
    
}
