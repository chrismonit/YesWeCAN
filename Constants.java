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
    public static String DEFAULT_PI = "0.25,0.25,0.25,0.25";
    public static double DEFAULT_SCALING = 1.0;
    public static double DEFAULT_OMEGA = 1.0;
    
    public static int HKY_IDENTIFIER = -1;
    public static int CAN0 = 0;
    public static int M1_IDENTIFIER = 1;
    public static int M2_IDENTIFIER = 2;

    public static int NUM_M1_SITE_CLASSES = 2;
    public static int NUM_M2_SITE_CLASSES = 3;
    
    public static String OMEGA_STRING = "w";
    public static String PROB_STRING = "p";
    public static String SITE_CLASS_0 = "0";
    public static String SITE_CLASS_2 = "2";
    
    public static String ERROR_PREFIX = "ERROR: ";
    
    public static String ARGS_DELIMITER = ",";
    
    public static char[] FRAMES = {'A', 'B', 'C'};
    
    public static String OUTPUT_DELIMITER = "\t";
    public static String LAYOUT = "LAYOUT"; // prefix for genetic structure diagram
    public static String HEADER = "HEADER"; // prefix for header in output
    public static String INITIAL = "INITIAL";
    public static String CALC = "CALC";
    public static String MLE = "MLE";
    public static String WITIHIN_FIELD_SEPARATOR = "_";
    
    // arguments for -fix option
    // e.g. "-fix kappa pi" will fix kappa and pi parameters, but omegas and scaling will be optimised
    public static String FIX_ALL = "all"; // compute lnL only
    public static String FIX_KAPPA = "kappa"; 
    public static String FIX_FREQUENCIES = "pi";
    public static String FIX_SCALING = "sc";
    
}
