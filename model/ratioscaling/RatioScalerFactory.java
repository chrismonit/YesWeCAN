/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package yeswecan.model.ratioscaling;

import yeswecan.model.ratioscaling.RatioScaler;
import yeswecan.model.ratioscaling.ProportionScaler;

/**
 *
 * @author Christopher Monit <c.monit.12@ucl.ac.uk>
 */
public class RatioScalerFactory {
     public static RatioScaler getRatioScaler(){
        return new ProportionScaler();
     }
}
