/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.utils;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class ArrayPrinter {
    
    public static String toString(double[] array, String delim){
        String str = "";
        for (int i = 0; i < array.length; i++) {
            str +=  Double.toString(array[i]) + delim;
        }
        return str;
    }
    
    public static void print(double[] array, String delim){
        System.out.println(toString(array,delim));
    }
 
      public static String toString(int[] array, String delim){
        String str = "";
        for (int i = 0; i < array.length; i++) {
            str +=  Integer.toString(array[i]) + delim;
        }
        return str;
    }
    
    public static void print(int[] array, String delim){
        System.out.println(toString(array,delim));
    }
    
    
}
