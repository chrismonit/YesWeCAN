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
    
    public static final double EPSILON = 1e-15;
    
    public static final int RELATIVE_TRESHOLD = -1;
    public static final double ABSOLUTE_THRESHOLD = 5.0e-6;
    public static final int MAXIMUM_EVALUATIONS = 10000;
    public static final int MAXIMUM_ITERATIONS = 10000;
    
    public static final double DEFAULT_KAPPA = 1.0;
    public static final String DEFAULT_PI = "0.25,0.25,0.25,0.25";
    public static final double DEFAULT_SCALING = 1.0;
    public static final double DEFAULT_OMEGA = 1.0;
    public static final double DEFAULT_OMEGA_0 = 0.99;
    public static final double DEFAULT_OMEGA_2 = 1.01;
    
    public static final int HKY_IDENTIFIER = -1;
    public static final int CAN0_IDENTIFIER = 0;
    public static final int M1_IDENTIFIER = 1;
    public static final int M2_IDENTIFIER = 2;
    public static final int CODON_FREQ_IDENTIFIER = 3;
    public static final int CODON_FREQ_MIX1_IDENTIFIER = 6; // TODO renumber these models!
    public static final int CODON_FREQ_MIX2_IDENTIFIER = 7;
    public static final int CAN_SUM_IDENTIFIER = 4;
    public static final int CAN_FREQ_IDENTIFIER = 5;

    public static final int PROPORTION_SCALER_IDENTIFIER = 0;
    public static final int CODON_SCALER_IDENTIFIER = 1;
    
    public static final int NUM_M1_SITE_CLASSES = 2;
    public static final int NUM_M2_SITE_CLASSES = 3;
    
    public static final String OMEGA_STRING = "w";
    public static final String PROB_STRING = "p";
    public static final String SITE_CLASS_0 = "0";
    public static final String SITE_CLASS_2 = "2";
    
    public static final String ERROR_PREFIX = "ERROR: ";
    
    public static final String ARGS_DELIMITER = ",";
    
    public static final char[] FRAMES = {'A', 'B', 'C'};
    public static final String[] SITE_TYPES = {"ALPHA", "BETA", "GAMMA"};
    
    // for output text
    public static final String DEL = "\t";
    public static final String LAYOUT = "LAYOUT"; // prefix for genetic structure diagram
    public static final String HEADER = "HEADER"; // prefix for header in output
    public static final String INITIAL = "INITIAL";
    public static final String CALC = "CALC";
    public static final String MLE = "MLE";
    public static final String TIME = "TIME";
    public static final String WITIHIN_FIELD_SEPARATOR = "_";
    public static final String TREE_PATH = "TREE_PATH";
    public static final String ALIGN_PATH = "ALIGN_PATH";
    public static final String CODON_FREQ_PATH = "CODON_FREQ_PATH";
    public static final String SIMULATION = "SIM";
    public static final String NU = "NU";
    public static final String P_N = "P_N";
    
    // for site class output tables
    public static final String MAJOR_DELIM = "|";
    public static final String CLASSES = "CLASSES"; // prefix for table rows for true site classes in mixmodel simulation
    public static final String SITES = "SITES";
    public static final int DEC_PLACES = 4;
    public static final String FALSE_STRING = "-";
    public static final String TRUE_STRING = "+";
    public static final double PROB_THRESH = 0.95;
    public static final String THRESH_MARKER = "_X";
    
    
    // arguments for -fix option
    // e.g. "-fix kappa pi" will fix kappa and pi parameters, but omegas and scaling will be optimised
    public static final String FIX_ALL = "all"; // compute lnL only
    public static final String FIX_KAPPA = "kappa"; 
    public static final String FIX_FREQUENCIES = "pi";
    public static final String FIX_SCALING = "sc";
}
