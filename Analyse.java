
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
import yeswecan.model.hky.FunctionHKY;
import yeswecan.model.hky.MutationModel;
import yeswecan.model.SubstitutionModel;
import yeswecan.optim.Optimise;
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
 
 data: create tree and alignment instances
 model: construct list of Parameter instances from input options 
 (these make up the model, including whether variables are fixed)
 
 fitHKY model to data
 
 produce output
     * 
     */
    public static void main(String[] args) {
        
        new Analyse(args);

    }// main
    
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
        
        loadData(this.comArgs.alignment(), this.comArgs.tree());
        //TODO
        // params = makeHKY(...)
        // fitHKY(data, model)
        
        // calculate lnL with fixed params
        if (this.comArgs.fix().contains("all")){
            calculateFixed(
                    makeHKY(),
                    this.tree,
                    this.alignment
            ); 
        }
        else{
            fitHKY();
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
    public List<Parameter> makeHKY(){
        TsTvRatio kappa = new TsTvRatio(this.comArgs.kappa());
        
        double[] frequencies = new double[States.NT_STATES]; // will be in correct order, whatever that may be
        
        if (Boolean.parseBoolean(this.comArgs.tcag())){
            frequencies = ReorderFrequencies.pamlToAlpha(this.comArgs.pi());
        }
        else{
            frequencies = this.comArgs.pi();
        }

        BaseFrequencies pi = new BaseFrequencies(frequencies);
        return Arrays.asList(kappa, pi);        
    }
    
    
    public static void calculateFixed(List<Parameter> model, Tree tree, AdvancedAlignment alignment){
        double[] optimisableParams = Mapper.getOptimisable(model); // map parameters to optimisation space, so FunctionHKY.value can use them
        FunctionHKY calculator = new FunctionHKY(alignment, tree);
        double lnL = calculator.value(optimisableParams);
        System.out.println("lnL: " + lnL + " "); // better to have it print the input parameters too, so you can see input and output together
    }
    

    
    
    // start the optimisation
    public void fitHKY(){
        FunctionHKY optFunction = new FunctionHKY(this.alignment, this.tree);
        Optimise opt = new Optimise();
        SubstitutionModel result = opt.optNMS(optFunction, new MutationModel(makeHKY()));
        
        System.out.println("opt lnL: "+result.getLnL());
        System.out.println( result.getParameters().get(0).toString());
        System.out.println(result.getParameters().get(1).toString());
    }
    
    
}// class
