/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.io;

import java.util.ArrayList;
import java.util.List;
import yeswecan.Constants;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author cmonit1
 */
public class ProbSiteClassTable extends SiteClassTable {
    
    private double[][][] siteClassProbs;
    private int decimalPlaces = Constants.DEC_PLACES;
    private boolean addThresholdMarker = false;
    private double PROB_THRESH = Constants.PROB_THRESH;
    private String THRESH_MARKER = Constants.THRESH_MARKER;
    private boolean roundValues = false;
    
    public ProbSiteClassTable(
            double[][][] siteClassProbs, int numberOfSiteClasses, 
            boolean addThresholdMarker, boolean roundValues){
        super(numberOfSiteClasses, siteClassProbs[0].length);
        this.siteClassProbs = siteClassProbs;
        this.addThresholdMarker = addThresholdMarker;
        this.roundValues = roundValues;
        makeTable();
    }
    
    public ProbSiteClassTable(
            double[][][] siteClassProbs, int numberOfSiteClasses,
            boolean addThresholdMarker, boolean roundValues,
            String majorDelimiter, String minorDelimiter
    ){
        this(siteClassProbs, numberOfSiteClasses, addThresholdMarker, roundValues);
        this.MAJOR_DELIMITER = majorDelimiter;
        this.MINOR_DELIMITER = minorDelimiter;
    }

    @Override
    protected void makeTable(){
        super.prependHeader();
        
        for (int iSite = 0; iSite < this.siteClassProbs.length; iSite++) {
            
            List<String> row = new ArrayList<String>();
            row.add(Integer.toString(iSite+1)); 
            row.add(this.MAJOR_DELIMITER);
            for (int iGene = 0; iGene < numberOfGenes; iGene++) {

                for (int iSiteClass = 0; iSiteClass < numberOfSiteClasses; iSiteClass++) {
                    double prob = this.siteClassProbs[iSite][iGene][iSiteClass];
                    
                    double valueToPrint;
                    if (this.roundValues) {
                        valueToPrint = MatrixPrinter.roundDouble(prob, this.decimalPlaces);
                    }else {
                        valueToPrint = prob;
                    }

                    String marker = "";
                    if (this.addThresholdMarker && valueToPrint >= this.PROB_THRESH) {
                        marker += this.THRESH_MARKER;
                    }
                    
                    row.add( Double.toString(valueToPrint) + marker );
                }
                row.add(this.MAJOR_DELIMITER);
            }
            this.table.add(row);
        }
        
    }
    
    
    public static void main(String[] args){
        System.out.println("testing ProbSiteClassTable \n");
    
        double[][][] data = new double[][][]{ // 5 sites and 3 genes
            {{ 0.123456, 0.7891011, 0.94999999999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }},
            {{ 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }},
            {{ 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }},
            {{ 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }},
            {{ 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }},
        };
        int n = 3; // there are n site classes represented
        
        boolean addThreshMarker = false;
        boolean roundValues = false;
        ProbSiteClassTable table = new ProbSiteClassTable(data, n, addThreshMarker, roundValues);
        table.print(Constants.CLASSES);
        
        String dest = "/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/can/netbeans/probs.tsv";
        table.write(dest);
    }
    
    
}
