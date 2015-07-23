/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.cli;

import com.beust.jcommander.Parameter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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
    
    @Parameter(names = {"-kappa", "-k"}, description = "Transition/transversion rate ratio. Default=1")
    private Double kappa = 1.0;
    
    @Parameter(names = {"-frequencies", "-pi"}, description = "Stationary frequencies for nucleotides. Default=[0.25,0.25,0.25,0.25]")
    private List<Double> pi = new ArrayList<>(Arrays.asList(0.25, 0.25, 0.25, 0.25));
    
    
    
    @Parameter(names = {"-fix"}, description = "Parameters to be fixed at initial values")
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
        for (int i = 0; i < States.NT_STATES; i++) {
            toReturn[i] = pi.get(i).doubleValue();
        }
        return toReturn;
    }
    
    public List<String> fix(){
        return fix;
    }
    
    
}//class
