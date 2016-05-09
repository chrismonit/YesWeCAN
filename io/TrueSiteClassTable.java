/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.io;

import yeswecan.Constants;

/**
 *
 * @author cmonit1
 */
public class TrueSiteClassTable extends SiteClassTable {
    
    private String falseString = "F"; // default values
    private String trueString = "T";
    
    private int[][] siteClasses;
    private int numberOfSiteClasses; // put in super class?
    private int numberOfGenes;
    
    public TrueSiteClassTable(int[][] siteClasses, int numberOfSiteClasses){
        super();
        this.siteClasses = siteClasses;
        this.numberOfSiteClasses = numberOfSiteClasses;
        this.numberOfGenes = siteClasses[0].length;
        makeTable();
    }
    
    public TrueSiteClassTable(
            int[][] siteClasses, int numberOfSiteClasses, 
            String trueString, String falseString,
            String majorDelimiter, String minorDelimiter
    ){
        this(siteClasses, numberOfSiteClasses);
        this.falseString = falseString;
        this.trueString = trueString;
        this.majorDelimiter = majorDelimiter;
        this.minorDelimiter = minorDelimiter;
    }
    
    
    
    @Override
    protected void makeTable(){
        super.prependHeader(this.numberOfSiteClasses, this.numberOfGenes);
        
        // TODO parse data to string lists etc
        
    }
    
    
    public static void main(String[] args){
        System.out.println("hello");
    
        int[][] data = new int[][]{ // 5 sites and 3 genes
            {1, 0, 1},
            {1, 0, 1},
            {1, 0, 1},
            {1, 0, 1},
            {1, 0, 1},
        };
        int n = 3; // there are 2 site classes represented
        
        TrueSiteClassTable table = new TrueSiteClassTable(data, n);
        table.print(Constants.CLASSES);
        
        String dest = "/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/can/netbeans/true.csv";
        table.write(dest);
    }
    
    
}
