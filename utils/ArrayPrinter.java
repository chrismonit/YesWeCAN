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
    
    public static void print(double[] array, String delim){
        String str = "";
        for (int i = 0; i < 4; i++) {
            str += delim + Double.toString(array[i]);
        }
        System.out.println(str);
            
    
    }
    
}
