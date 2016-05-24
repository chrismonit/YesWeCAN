/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yeswecan.io.ntsiteclasses;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import yeswecan.Constants;

/**
 *
 * @author cmonit1
 */
public abstract class SiteClassTable {
    
    protected List<List<String>> table;
    protected String MAJOR_DELIMITER = Constants.MAJOR_DELIM;
    protected String MINOR_DELIMITER = Constants.DEL;
    protected String SITES_HEADER = Constants.SITES;
    
    protected int numberOfSiteClasses; // put in super class?
    protected int numberOfGenes;
    
    public SiteClassTable(int numberOfSiteClasses, int numberOfGenes){
        this.table = new ArrayList<List<String>>();
        this.numberOfSiteClasses = numberOfSiteClasses;
        this.numberOfGenes = numberOfGenes;
    }
    
    
    public SiteClassTable(List<List<String>> table){
        this.table = table;
    }
    
    
    public void addDelimiterColumn(String delimiter){
        for (List<String> stringList : this.table){
            stringList.add(delimiter);
        }        
    }
    // could have equivalent for adding rows
    
    public void write(String filePath){
        
        try{
            FileWriter fileWriter = new FileWriter(filePath);
            BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
            
            for (List<String> stringList : this.table){
                bufferedWriter.write(String.join(this.MINOR_DELIMITER, stringList));
                bufferedWriter.newLine();
            }
            
            bufferedWriter.close();
        
        }catch (IOException e){
            e.printStackTrace();
        }
        
    }
    
    public void print(String lineMarker){
        for (List<String> stringList : this.table){
            System.out.println(lineMarker+this.MINOR_DELIMITER+String.join(this.MINOR_DELIMITER, stringList));
            //System.out.println(String.join(this.minorDelimiter, stringList));
        }
    }
    
    protected abstract void makeTable();
    
    protected void prependHeader(){
        List<String> header = new ArrayList<String>();
        header.add(this.SITES_HEADER); // blank (could name column SITES or something)
        header.add(this.MAJOR_DELIMITER);
        for (int iGene = 1; iGene < this.numberOfGenes; iGene++) { // start at 1 to ignore noncoding
            
            for (int iSiteClass = 0; iSiteClass < this.numberOfSiteClasses; iSiteClass++) {
                header.add(Integer.toString(iGene));
            }
            header.add(this.MAJOR_DELIMITER);
        }
        
        this.table.add(0, header);
    }
    
    
}
