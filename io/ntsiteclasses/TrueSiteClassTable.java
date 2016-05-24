/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.io.ntsiteclasses;

import java.util.ArrayList;
import java.util.List;
import yeswecan.Constants;

/**
 *
 * @author cmonit1
 */
public class TrueSiteClassTable extends SiteClassTable {
    
    private String FALSE_STRING = Constants.FALSE_STRING; // default values
    private String TRUE_STRING = Constants.TRUE_STRING;
    
    private int[][] siteClasses;

    
    public TrueSiteClassTable(int[][] siteClasses, int numberOfSiteClasses){
        super(numberOfSiteClasses, siteClasses[0].length);
        this.siteClasses = siteClasses;
        makeTable();
    }
    
    public TrueSiteClassTable(
            int[][] siteClasses, int numberOfSiteClasses, 
            String majorDelimiter, String minorDelimiter
    ){
        this(siteClasses, numberOfSiteClasses);
        this.MAJOR_DELIMITER = majorDelimiter;
        this.MINOR_DELIMITER = minorDelimiter;
    }
    
    
    
    @Override
    protected void makeTable(){
        super.prependHeader();
        
        for (int iSite = 0; iSite < this.siteClasses.length; iSite++) {
            
            List<String> row = new ArrayList<String>();
            row.add(Integer.toString(iSite+1)); 
            row.add(this.MAJOR_DELIMITER);
            for (int iGene = 1; iGene < this.numberOfGenes; iGene++) { // start at 1 to avoid noncoding

                for (int iSiteClass = 0; iSiteClass < numberOfSiteClasses; iSiteClass++) {
                    String fieldText;
                    if (iSiteClass == this.siteClasses[iSite][iGene]) {
                        fieldText = this.TRUE_STRING;
                    }else {
                        fieldText = this.FALSE_STRING;
                    }
                    
                    row.add( fieldText );
                }
                row.add(this.MAJOR_DELIMITER);
            }
            this.table.add(row);
        }
        
    }
    
    
    
    
    public static void main(String[] args){
        System.out.println("testing \n");
    
        int[][] data = new int[][]{ // 5 sites and 3 genes
            {2, 0, 1},
            {1, 0, 2},
            {1, 0, 1},
            {1, 0, 2},
            {1, 0, 0},
        };
        int n = 3; // there are 2 site classes represented
        
        TrueSiteClassTable table = new TrueSiteClassTable(data, n);
        table.print(Constants.CLASSES);
        
        String dest = "/Users/cmonit1/Desktop/overlapping_ORF/CAN_model/YesWeCAN/test/can/netbeans/true.csv";
        table.write(dest);
    }
    
    
}
