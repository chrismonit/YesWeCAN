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
import java.util.List;
import yeswecan.Constants;
import yeswecan.phylo.States;
/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class CommandArgs {
    
    
   @Parameter(names = {"-sequences", "-s"}, required = true, description = "Phylip format nucleotide sequence alignment")
    private String alignmentPath = "NO ALIGNMENT PATH";
    
    @Parameter(names = {"-tree", "-t"}, required = true, description = "Newick format tree file")
    private String treePath = "NO TREE PATH";
    

    // initial values for params
    
    @Parameter(names = {"-kappa", "-k"}, required = false, description = "Transition/transversion rate ratio. Default=1")
    private Double kappa = Constants.DEFAULT_KAPPA;
    
    @Parameter(names = {"-frequencies", "-pi"}, required = false, description = "Stationary frequencies for nucleotides, delimited by comma. Default=\"0.25,0.25,0.25,0.25\"")
    private String piString = Constants.DEFAULT_PI;
    
    @Parameter(names = {"-scaling", "-sc"}, required = false, description = "Scaling factor for branch lengths.")
    private double scaling = Constants.DEFAULT_SCALING;
    
    @Parameter(names = {"-model", "-m"}, required = true, description = "Model of selection to use. -1 = HKY, 0 = CAN original, 1 = M1a (negative/neutral), 2 = M2a (negative, neutral/positive)")
    private int model = Constants.CAN0_IDENTIFIER;
    
    
    @Parameter(names = {"-omega", "-w"}, required = false, description = "Initial omega values for each gene. Delimited by comma")
    private String omegasArgument = "";
    
    @Parameter(names = {"-omega0", "-w0"}, required = false, description = "Initial values for each gene's w_0. Delimited by comma")
    private String omegaArg0 = "";
    
    @Parameter(names = {"-omega2", "-w2"}, required = false, description = "Initial omega values for each gene's_2 (ignored if running M1). Delimited by comma")
    private String omegaArg2 = "";
    
    //NB w_1 is fixed to 1.0, so doesn't need an argument
    
    
    @Parameter(names = {"-prob0", "-p0"}, required = false, description = "Initial values for each gene's p_0. Delimited by comma")
    private String probArg0 = "";
    
    @Parameter(names = {"-prob1", "-p1"}, required = false, description = "Initial values for each gene's p_1. Delimited by comma")
    private String probArg1 = "";
    
    @Parameter(names = {"-prob2", "-p2"}, required = false, description = "Initial values for each gene's p_2 (Ignored if runing M2). Delimited by comma")
    private String probArg2 = "";
    
    
    
    @Parameter(names = {"-frameA", "-a"}, required = false, description = "Gene layout for frame A, delimited by comma")
    private String aFrame = "";
    
    @Parameter(names = {"-frameB", "-b"}, required = false, description = "Gene layout for frame B, delimited by comma")
    private String bFrame = "";
    
    @Parameter(names = {"-frameC", "-c"}, required = false, description = "Gene layout for frame C, delimited by comma")
    private String cFrame = "";
    
    @Parameter(names = {"-lengths", "-l"}, required = true, description = "Lengths for each partition, delimted by comma")
    private String lengths = ""; 
    
   
    
    @Parameter(names = {"-tcag"}, required = false, description = "Input frequencies are ordered: T,C,A,G (PAML style)")
    private String tcag = "false";
    
    @Parameter(names = {"-phy"}, required = false, description = "Sequence data are in Phylip format (assumes Fasta format by default")
    private boolean phy = false;
    
    @Parameter(names = {"-verbose", "-v"}, required = false, description = "Print additional information during program run")
    private boolean verbose = false;
    

    @Parameter(names = {"-fix"}, variableArity = true, required = false, description = "Parameters to be fixed at initial values")
    private List<String> fix = new ArrayList<>();
    
    public int getModel(){
        return model;
    }
    
    public List<String> fix(){
        return fix;
    }
    
    public String tcag(){
        return tcag;
    }
  
    public boolean phy(){
        return phy;
    }
    
    public boolean verbose(){
        return phy;
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
    
    private double[] piValues;
    
    
    public double[] pi(){ 
        if (this.piValues == null){
            piValues = new double[States.NT_STATES];
            String[] piValueStrings;

            try{
                piValueStrings = this.piString.split(Constants.ARGS_DELIMITER);
                if (piValueStrings.length != States.NT_STATES){
                    String errorMsg = MessageFormat.format( "-frequences: Must have {0} floating point values, delimited by \"{1}\"", States.NT_STATES, Constants.ARGS_DELIMITER );
                    throw new ParameterException(errorMsg);
                }

                for (int i = 0; i < States.NT_STATES; i++) {
                    piValues[i] = Double.parseDouble(piValueStrings[i]);
                }
            }
            catch(Exception e){
                e.printStackTrace();
                System.exit(1);
            }
        }// if
        return this.piValues; 

    }// piString
    
     
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
    
    private int geneNumber = -1;
    // initially set to -1 above. Don't want to call findGeneNumber each time, 
    // so we set geneNumber value the first (and only) time findGeneNumber is called
    public int getGeneNumber(){
        if (this.geneNumber == -1){
            this.geneNumber = findGeneNumber();
        }
        return this.geneNumber;
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
    

    
//    public double[] geneSpecificParameter(){
//        double[] values = new double[getGeneNumber()];
//        if ("".equals(this.omegasArgument)){ // no starting omegasArgument have been supplied by user
//            for (int i = 0; i < values.length; i++) {
//                values[i] = Constants.DEFAULT_OMEGA;
//            }
//        }
//        else{
//            String[] omegasStrings = this.omegasArgument.split(Constants.ARGS_DELIMITER);
//            for (int i = 0; i < values.length; i++) {
//                values[i] = Double.parseDouble(omegasStrings[i]);
//            }
//        }
//        return values;
//    }
    
    // for CANM
    public double[] omegas(){
        return geneSpecificParameter(this.omegasArgument, Constants.DEFAULT_OMEGA);
    }
    
    // for CAN mixture
    public double[] omega0(){
        return geneSpecificParameter(this.omegaArg0, Constants.DEFAULT_OMEGA_0);
    }
    
    public double[] omega2(){
        return geneSpecificParameter(this.omegaArg2, Constants.DEFAULT_OMEGA_2);
    }
    
    private double[] prob(String argument){
        if (this.model == Constants.M1_IDENTIFIER){
            return geneSpecificParameter(argument, 1.0/(double)Constants.NUM_M1_SITE_CLASSES);
        }
        else if (this.model == Constants.M2_IDENTIFIER){
            return geneSpecificParameter(argument, 1.0/(double)Constants.NUM_M2_SITE_CLASSES);
        }
        else{
            throw new RuntimeException("Model probability values requested from CommandArgs, when model specified does not include probabilities");
        }
    }
    
    public double[] prob0(){
        return prob(this.probArg0);
    }
    
    public double[] prob1(){
        return prob(this.probArg1);
    }
    
    public double[] prob2(){
        return prob(this.probArg2);
    }
    
    // for omegas or prob values. argument is string supplied at CLI after the relevant flag
    private double[] geneSpecificParameter(String argument, double defaultValue){
        
        double[] values = new double[getGeneNumber()];
        try{
            if ("".equals(argument)){ // no starting argument have been supplied by user
                for (int i = 0; i < values.length; i++) {
                    values[i] = defaultValue;
                }
            }
            else{
                String[] omegasStrings = argument.split(Constants.ARGS_DELIMITER);
                for (int i = 0; i < values.length; i++) {
                    values[i] = Double.parseDouble(omegasStrings[i]);
                }
            }
        }catch(ArrayIndexOutOfBoundsException e){
            System.out.println("Could be a mismatch between number of -w or -p arguments specified and the number of genes?");
            e.printStackTrace();
            System.exit(1);
        }
        return values;
    }
    
    

    
    
}//class