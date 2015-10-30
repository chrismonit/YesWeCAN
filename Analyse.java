
package yeswecan;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import pal.alignment.Alignment;
import pal.alignment.AlignmentReaders;
import pal.alignment.SimpleAlignment;
import pal.datatype.Nucleotides;
import pal.tree.ReadTree;
import pal.tree.Tree;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Parameter;
import swmutsel.model.parameters.TsTvRatio;
import yeswecan.cli.CommandArgs;
import yeswecan.model.Function;
import yeswecan.phylo.AdvancedAlignment;
import yeswecan.phylo.ReorderFrequencies;
import yeswecan.phylo.States;
import yeswecan.utils.ArrayPrinter;

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
        //TODO
        // params = makeModel(...)
        // fit(data, model)
        
        // calculate lnL with fixed params
        if (Boolean.parseBoolean(comArgs.fix())){
            calculateFixed(
                    makeModel(comArgs),
                    this.tree,
                    this.alignment
            ); 
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
    
    
    // only HKY at the moment
    public static List<Parameter> makeModel(CommandArgs comArgs){
        TsTvRatio kappa = new TsTvRatio(comArgs.kappa());
        
        double[] frequencies = new double[States.NT_STATES]; // will be in correct order, whatever that may be
        
        if (Boolean.parseBoolean(comArgs.tcag())){
            frequencies = ReorderFrequencies.pamlToAlpha(comArgs.pi());
        }
        else{
            frequencies = comArgs.pi();
        }

        BaseFrequencies pi = new BaseFrequencies(frequencies);
        
        return Arrays.asList(kappa, pi);        
    }
    
    public static void calculateFixed(List<Parameter> model, Tree tree, AdvancedAlignment alignment){
        double[] optimisableParams = Mapper.getOptimisable(model); // map parameters to optimisation space, so Function.value can use them
        Function calculator = new Function(alignment, tree);
    
        double lnL = calculator.value(optimisableParams);
        System.out.println("lnL: " + lnL); // better to have it print the input parameters too, so you can see input and output together
    }
    
    
    //TODO
    // start the optimisation
    public void fit(){
        // make new instance of function, which takes in the data
        //use Optimise to fit the model to the data
        // get the MLEs
        // format in some way
    }
    
    
}// class
