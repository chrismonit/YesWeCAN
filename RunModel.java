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
public abstract class RunModel {
    
    abstract String getHeader();
    
    abstract double[] fit();
    
    abstract double[] calculatelnL();
    
    
}
