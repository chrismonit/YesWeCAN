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
import yeswecan.utils.ArrayPrinter;

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
        
    public Analyse(String[] args){
        
        this.comArgs = new CommandArgs();
        JCommander jcom = new JCommander(this.comArgs);
        
        try{
            jcom.parse(args);
        }
        catch(ParameterException ex){
            System.out.println(ex.getMessage());
        }
                
        loadData(this.comArgs.alignment(), this.comArgs.tree(), this.comArgs.phy());
 
        RunModel run;
        
        switch (this.comArgs.getModel()){
            case Constants.HKY_IDENTIFIER: run = new RunHKY(alignment, tree, this.comArgs);
                break;
            case Constants.CAN0_IDENTIFIER: run = new RunCAN(alignment, tree, this.comArgs);
                break;
            case Constants.M1_IDENTIFIER: run = new RunCANMixture(alignment, tree, this.comArgs, Constants.M1_IDENTIFIER);
                break;
            case Constants.M2_IDENTIFIER: run = new RunCANMixture(alignment, tree, this.comArgs, Constants.M2_IDENTIFIER);
                break;
            default: throw new RuntimeException(Constants.ERROR_PREFIX + "Invalid model argument (-m)");

        }
        
        System.out.println(Constants.HEADER + Constants.DEL + String.join(Constants.DEL, run.getHeader()));
        
        if (this.comArgs.fix().contains(Constants.FIX_ALL)){
            double[] result = run.calculate();
            
            System.out.println( Constants.CALC + Constants.DEL +
                ArrayPrinter.toString(result, Constants.DEL) );
        }
        else{
            double[] initial = run.getInitialValues();
            System.out.println(Constants.INITIAL + Constants.DEL + 
                    ArrayPrinter.toString(initial, Constants.DEL));
            
            long start = System.currentTimeMillis();
            
            double[] result = run.fit();
            
            long runTime = System.currentTimeMillis() - start;
            
            System.out.println( Constants.MLE + Constants.DEL +
                    ArrayPrinter.toString(result, Constants.DEL) ); 
            
            System.out.println("");
            
            double seconds = (double)runTime/1000.0;
            double minutes = seconds/60.0;
            double hours = minutes/60.0;
            double days = hours/24.0;
            String[] time = new String[]{ Double.toString(seconds), Double.toString(minutes), Double.toString(hours), Double.toString(days) };
            
            System.out.println(Constants.TIME + Constants.DEL + String.join(Constants.DEL, time));
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
    
    
}
