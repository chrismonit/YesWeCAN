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
    

    
    
    
    @Parameter(names = {"-fix"}, required = false, description = "Parameters to be fixed at initial values")
    private List<String> fix = new ArrayList<>();
    
    
    public String alignment(){
        return alignmentPath;
    }
    
    public String tree(){
        return treePath;
    }
    
    public double kappa(){
        return kappa.doubleValue();
    }
    
    public double[] pi(){ //just makes a double[] from ArrayList<Double>
        double[] toReturn = new double[States.NT_STATES];
        String[] piValueStrings;
        
        try{
            piValueStrings = pi.split(Constants.PI_DELIMITER);
            if (piValueStrings.length != States.NT_STATES){
                String errorMsg = MessageFormat.format( "-frequences: Must have {0} floating point values, delimited by \"{1}\"", States.NT_STATES, Constants.PI_DELIMITER );
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
    
    public List<String> fix(){
        return fix;
    }
    
    
}//class
