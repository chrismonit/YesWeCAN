/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class Constants {
    
    public static double EPSILON = 1e-10;
    
    public static int RELATIVE_TRESHOLD = -1;
    public static double ABSOLUTE_THRESHOLD = 5.0e-6;
    public static int MAXIMUM_EVALUATIONS = 4000;
    public static int MAXIMUM_ITERATIONS = 10000;
    
    public static double DEFAULT_KAPPA = 1.0;
    public static double[] DEFAULT_PI = new double[]{ 0.25, 0.25, 0.25, 0.25 };
    public static double DEFAULT_SCALING = 1.0;
    
    // may want to replace this with default values in the jcommander class
    
    public static String ERROR_PREFIX = "ERROR: ";
    
    public static String ARGS_DELIMITER = ",";
    
    public static char[] FRAMES = {'A', 'B', 'C'};
    
    public static String LAYOUT = "LAYOUT"; // prefix for genetic structure diagram
}
