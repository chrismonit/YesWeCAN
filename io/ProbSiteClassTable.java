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
    
    public ProbSiteClassTable(double[][][] siteClassProbs, int numberOfSiteClasses){
        super(numberOfSiteClasses, siteClassProbs[0].length);
        this.siteClassProbs = siteClassProbs;
        makeTable();
    }
    
    public ProbSiteClassTable(
            double[][][] siteClassProbs, int numberOfSiteClasses,
            String majorDelimiter, String minorDelimiter
    ){
        this(siteClassProbs, numberOfSiteClasses);
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
                    double rounded = MatrixPrinter.roundDouble(this.siteClassProbs[iSite][iGene][iSiteClass], this.decimalPlaces);
                    row.add( Double.toString(rounded) );
                }
                row.add(this.MAJOR_DELIMITER);
            }
            this.table.add(row);
        }
        
    }
    
    
    public static void main(String[] args){
        System.out.println("testing ProbSiteClassTable \n");
    
        double[][][] data = new double[][][]{ // 5 sites and 3 genes
            {{ 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }},
            {{ 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }},
            {{ 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }},
            {{ 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }},
            {{ 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }, { 0.123456, 0.7891011, 0.08744289999999999 }},
        };
        int n = 3; // there are n site classes represented
        
        ProbSiteClassTable table = new ProbSiteClassTable(data, n);
        table.print(Constants.CLASSES);
        
        String dest = "/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/can/netbeans/probs.tsv";
        table.write(dest);
    }
    
    
}
