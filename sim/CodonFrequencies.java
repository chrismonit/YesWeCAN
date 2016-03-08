/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.sim;

import com.opencsv.CSVReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import yeswecan.utils.MatrixPrinter;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class CodonFrequencies {
    
    
    // take CSV instance
    // have get methods returning freq for given codon
    
    private double[][] frequencies;
    
    // for no gene (neutral) case. all codon frequencies = 1/64
    public CodonFrequencies(){
        this.frequencies = new double[16][4];
        for (int i = 0; i < this.frequencies.length; i++) {
            for (int j = 0; j < this.frequencies[0].length; j++) {
                this.frequencies[i][j] = 1.0/64.0;
            }
        }
    }
    
    
    public CodonFrequencies(String frequenciesFilePath){
        this.frequencies = new double[16][4];
        List<String[]> textData;

        try{
            CSVReader reader = new CSVReader(new FileReader(frequenciesFilePath));
            textData = reader.readAll();
        }catch(IOException e){
            System.out.println("Failed to read codon frequencies CSV");
            textData = null;
            e.printStackTrace();
        }

        try{
            for (int i = 0; i < 16; i++) {
                for (int j = 0; j < 4; j++) {
                    this.frequencies[i][j] = Double.parseDouble(textData.get(i)[j]);
                }
            }
        }catch(ArrayIndexOutOfBoundsException e){
            System.out.println("Codon frequencies file not formatted as expected: must be in 16x4 TCAG");
            e.printStackTrace();
            System.exit(1);
        }
        catch(NumberFormatException e){
            System.out.println("Cannot read codon frequency element as number. Please check codon frequency file");
            e.printStackTrace();
            System.exit(1);
        }
        
        //MatrixPrinter.PrintMatrix(frequencies, "frequencies");
    }
    
    /**
     * Assuming codons are presented out as TCAG in table,
     * and nucleotides {T, C, A, G} = { 0, 1, 2, 3 }
     * and a codon is made up of nucleotides { x, y, z }
     * then
     * i = 4x + z
     * j = y
     * and i and j can be used to access the codon table
     * 
     */
    public double getFrequency(int[] codon){
        int i =  4*codon[0] + codon[2];
        int j = codon[1];
        return this.frequencies[i][j];
    }
    

    
    
    public static void main(String[] args){
        System.out.println("Testing CodonFreuqneices");
        
        CodonFrequencies cf = new CodonFrequencies(args[0]);
        System.out.println(cf.getFrequency(new int[]{1,0,2}));
        System.out.println("End main");
    }
    
    
}
