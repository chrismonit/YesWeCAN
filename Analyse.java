
package yeswecan;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import java.io.FileReader;
import java.util.List;
import pal.alignment.Alignment;
import pal.alignment.AlignmentReaders;
import pal.alignment.SimpleAlignment;
import pal.datatype.Nucleotides;
import pal.tree.ReadTree;
import pal.tree.Tree;
import swmutsel.model.parameters.Parameter;
import yeswecan.cli.CommandArgs;
import yeswecan.phylo.AdvancedAlignment;

/**
 *
 * @author Christopher Monit (c.monit.12@ucl.ac.uk)
 */
public class Analyse {

    /**
     * Read in arguments
     * 
     * data: create tree and alignment instances
     * model: construct list of Parameter instances from input options 
     * (these make up the model, including whether variables are fixed)
     * 
     * fit model to data
     * 
     * produce output
     * 
     */
    public static void main(String[] args) {
        
        new Analyse(args);

    }// main
    
    
    private AdvancedAlignment alignment;
    private Tree tree;
    
    public Analyse(String[] args){
        CommandArgs comArgs = new CommandArgs();
        JCommander jcom = new JCommander(comArgs);
        
        try{
            jcom.parse(args);
        }
        catch(ParameterException ex){
            System.out.println(ex.getMessage());
        }
        
        loadData(comArgs.alignment(), comArgs.tree());
        
    }
    
    //TODO make more sophistcated exceptions to help user find problems. Separate tree and alignment reading in 
    public void loadData(String alignmentPath, String treePath){
        try{
            alignment = new AdvancedAlignment(
                                new SimpleAlignment(
                                        AlignmentReaders.readFastaSequences(new FileReader(alignmentPath), new Nucleotides())));
            tree = new ReadTree(treePath);
        }
        catch(Exception e){
            System.out.println(Constants.ERROR_PREFIX + "Unable to load alignment or tree file(s)");
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    
    
    public List<Parameter> makeModel(CommandArgs comArgs){
        
    
    }
    
    
    
    public void fit(){
        
    }
    
    
}// class
