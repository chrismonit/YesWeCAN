/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.cli;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import java.text.MessageFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import yeswecan.Constants;
import yeswecan.phylo.States;
import yeswecan.utils.ArrayPrinter;
import yeswecan.utils.MatrixPrinter;
/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class CommandArgs {
    
//    @Parameter
//    private List<String> parameters = new ArrayList<>();
    
    //data
    
   @Parameter(names = {"-sequences", "-s"}, required = true, description = "Phylip format nucleotide sequence alignment")
    private String alignmentPath = "NO ALIGNMENT PATH";
    
    @Parameter(names = {"-tree", "-t"}, required = true, description = "Newick format tree file")
    private String treePath = "NO TREE PATH";
    

    // initial values for params
    
    @Parameter(names = {"-kappa", "-k"}, required = false, description = "Transition/transversion rate ratio. Default=1")
    private Double kappa = Constants.DEFAULT_KAPPA;
    
    @Parameter(names = {"-frequencies", "-pi"}, required = false, description = "Stationary frequencies for nucleotides, delimited by comma. Default=\"0.25,0.25,0.25,0.25\"")
    private String pi = "0.25,0.25,0.25,0.25";
    
    @Parameter(names = {"-scaling", "-sc"}, required = true, description = "Scaling factor for branch lengths.")
    private double scaling = Constants.DEFAULT_SCALING;
    
    //@Parameter(names = {"-omegas", "-w"}, required = false, description = "Initial omega values for each gene. Delimited by comma")
    //private String omegasArgument = "";
    
    @Parameter(names = {"-omega0", "-w0"}, required = false, description = "Initial values for each gene's w_0. Delimited by comma")
    private String omegaArg0 = "";
    
    @Parameter(names = {"-omega2", "-w2"}, required = false, description = "Initial omega values for each gene's_2 (ignored if running M1). Delimited by comma")
    private String omegaArg2 = "";
    
    //NB w_1 is fixed to 1.0, so doesn't need an argument
    
    @Parameter(names = {"-frameA", "-a"}, required = true, description = "Gene layout for frame A, delimited by comma")
    private String aFrame = "";
    
    @Parameter(names = {"-frameB", "-b"}, required = true, description = "Gene layout for frame B, delimited by comma")
    private String bFrame = "";
    
    @Parameter(names = {"-frameC", "-c"}, required = true, description = "Gene layout for frame C, delimited by comma")
    private String cFrame = "";
    
    
    @Parameter(names = {"-lengths", "-l"}, required = true, description = "Lengths for each partition, delimted by comma")
    private String lengths = ""; // can't have default since don't know how many genes there are
    
   
    
    @Parameter(names = {"-tcag"}, required = false, description = "Input frequencies are ordered: T,C,A,G (PAML style)")
    private String tcag = "false";
    
    @Parameter(names = {"-phy"}, required = false, description = "Sequence data are in Phylip format (assumes Fasta format by default")
    private String phy = "false";
    
    
    
//    @Parameter(names = {"-fix"}, required = false, description = "Parameters to be fixed at initial values")
//    private String fix = "false"; // need to change so no argument required, just looks whether the flag is present or not
//    //private List<String> fix = new ArrayList<>();
    
    @Parameter(names = {"-fix"}, variableArity = true, required = false, description = "Parameters to be fixed at initial values")
    private List<String> fix = new ArrayList<>();
    
    public List<String> fix(){
        return fix;
    }
    
    public String alignment(){
        return alignmentPath;
    }
    
    public String tree(){
        return treePath;
    }
    
    public double scaling(){
        return scaling;
    }
    
    public double kappa(){
        return kappa.doubleValue(); // does doubleValue() here actually do anything??
    }
    
    public double[] pi(){ 
        double[] toReturn = new double[States.NT_STATES];
        String[] piValueStrings;
        
        try{
            piValueStrings = pi.split(Constants.ARGS_DELIMITER);
            if (piValueStrings.length != States.NT_STATES){
                String errorMsg = MessageFormat.format( "-frequences: Must have {0} floating point values, delimited by \"{1}\"", States.NT_STATES, Constants.ARGS_DELIMITER );
                throw new ParameterException(errorMsg);
            }
            
            for (int i = 0; i < States.NT_STATES; i++) {
                toReturn[i] = Double.parseDouble(piValueStrings[i]);
            }
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        return toReturn;
    }// pi
    

    
    public int[] aFrame(){
        return stringToIntArray(this.aFrame, Constants.ARGS_DELIMITER);
    }
    
    public int[] bFrame(){
        return stringToIntArray(this.bFrame, Constants.ARGS_DELIMITER);
    }
    
    public int[] cFrame(){
        return stringToIntArray(this.cFrame, Constants.ARGS_DELIMITER);
    }
    
    public int[] lengths(){
        return stringToIntArray(this.lengths, Constants.ARGS_DELIMITER);
    }
    
    private static int[] stringToIntArray(String argument, String delimiter){
        String[] argumentArray = argument.split(delimiter);
        int[] toReturn = new int[argumentArray.length];
        for (int i = 0; i < argumentArray.length; i++) {
            toReturn[i] = Integer.parseInt(argumentArray[i]);
        }
        return toReturn;
    }
    
    // TODO this is very inefficient
    public int getGeneNumber(){
        return findGeneNumber();
    }
    
    private int findGeneNumber(){
        int[][] allGenes = { 
            stringToIntArray(this.aFrame, Constants.ARGS_DELIMITER),
            stringToIntArray(this.bFrame, Constants.ARGS_DELIMITER),
            stringToIntArray(this.cFrame, Constants.ARGS_DELIMITER)
        };
          
        int max = 0;
        for (int i = 0; i < allGenes.length; i++) {
            for (int j = 0; j < allGenes[0].length; j++) {
                max = Math.max(max, allGenes[i][j]);
            }
        }
        return max;
    }
    
    
//    public double[] omegas(){
//        double[] omegaValues = new double[getGeneNumber()];
//        if ("".equals(this.omegasArgument)){ // no starting omegasArgument have been supplied by user
//            for (int i = 0; i < omegaValues.length; i++) {
//                omegaValues[i] = Constants.DEFAULT_OMEGA;
//            }
//        }
//        else{
//            String[] omegasStrings = this.omegasArgument.split(Constants.ARGS_DELIMITER);
//            for (int i = 0; i < omegaValues.length; i++) {
//                omegaValues[i] = Double.parseDouble(omegasStrings[i]);
//            }
//        }
//        return omegaValues;
//    }
    
    public double[] omega0(){
        return omegas(this.omegaArg0);
    }
    
    public double[] omega2(){
        return omegas(this.omegaArg2);
    }
    
    
    private double[] omegas(String omegasArgument){
        double[] omegaValues = new double[getGeneNumber()];
        if ("".equals(omegasArgument)){ // no starting omegasArgument have been supplied by user
            for (int i = 0; i < omegaValues.length; i++) {
                omegaValues[i] = Constants.DEFAULT_OMEGA;
            }
        }
        else{
            String[] omegasStrings = omegasArgument.split(Constants.ARGS_DELIMITER);
            for (int i = 0; i < omegaValues.length; i++) {
                omegaValues[i] = Double.parseDouble(omegasStrings[i]);
            }
        }
        return omegaValues;
    }
    
    
  public String tcag(){
        return tcag;
    }
  
  public String phy(){
        return phy;
    }
    
//    public List<String> fix(){
//        return fix;
//    }
    
    
}//class
