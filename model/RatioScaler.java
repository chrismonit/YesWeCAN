/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public interface RatioScaler {
    
    public double get(double ratio, int siteType, int frame);
    
}
