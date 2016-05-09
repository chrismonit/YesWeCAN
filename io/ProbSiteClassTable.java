/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.io;

/**
 *
 * @author cmonit1
 */
public class ProbSiteClassTable extends SiteClassTable {
    
    private double[][] siteClassProbs;
    
    public ProbSiteClassTable(double[][] siteClassProbs, int numberOfSiteClasses){
        super();
        this.siteClassProbs = siteClassProbs;
        this.numberOfSiteClasses = numberOfSiteClasses;
        this.numberOfGenes = siteClassProbs[0].length;
        // makeTable();
    }
    
    public ProbSiteClassTable(
            double[][] siteClassProbs, int numberOfSiteClasses,
            String majorDelimiter, String minorDelimiter
    ){
        this(siteClassProbs, numberOfSiteClasses);
        this.MAJOR_DELIMITER = majorDelimiter;
        this.MINOR_DELIMITER = minorDelimiter;
    }

    
    
    
    
}
