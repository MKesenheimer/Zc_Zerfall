/*
 Copyright (C) 2013-1024 Matthias Kesenheimer
 m.kesenheimer@gmx.net
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 Dieses Programm ist Freie Software: Sie können es unter den Bedingungen
 der GNU General Public License, wie von der Free Software Foundation,
 Version 3 der Lizenz oder (nach Ihrer Wahl) jeder neueren
 veröffentlichten Version, weiterverbreiten und/oder modifizieren.
 
 Dieses Programm wird in der Hoffnung, dass es nützlich sein wird, aber
 OHNE JEDE GEWÄHRLEISTUNG, bereitgestellt; sogar ohne die implizite
 Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.
 Siehe die GNU General Public License für weitere Details.
 
 Sie sollten eine Kopie der GNU General Public License zusammen mit diesem
 Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.
 */

//hier quadriertes Matrixelement in Abhängigkeit von eps, lam, s12, s23 und den Formfaktoren eingeben:

double getMinv(double eps, double lam, double s12, double s31, double F1, double F23, double F4, double F5, double F6) {

      
      double s23 = Power(M(eps),2) + Power(m1,2) + Power(m2,2) + Power(m3,2) - s12 - s31;
      double minv;
      /*********************************************************************/
      /*                                                                   */
      /*        Photonmasse gleich Null                                    */
      /*                                                                   */
      /*********************************************************************/
      
      if (mGam == 0.) {
      minv =
      
      //------------------------------------------------------------------------------------------------------------------------------------------------
      (Power(F23,2)*Power(m1,2)*Power(-Power(m1,2) + s12 + s23,2) +
       4*Power(F1,2)*(-5*Power(m2,2) + Power(-Power(m1,2) + s12 + s23,2)/(2.*Power(M(eps),2))) -
       4*Power(F23,2)*Power(m1,2)*Power(m2,2)*Power(M(eps),2) - Power(F23,2)*Power(Power(m1,2) + Power(m2,2) - s12,2)*Power(M(eps),2) -
       Power(F23,2)*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) - s12 - s23)*(Power(m1,2) - s23 + Power(M(eps),2)) +
       Power(F5,2)*Power(m1,2)*(-8*Power(m2,4) + 12*Power(m2,2)*(-Power(m1,2) + s12 + s23) - 2*Power(-Power(m1,2) + s12 + s23,2) -
                                (Power(m2,2)*Power(-Power(m1,2) + s12 + s23,2))/Power(M(eps),2) - 4*Power(m2,2)*Power(M(eps),2)) -
       Power(F4,2)*((Power(m1,2) - s12 - s23)*(Power(m1,2)*Power(m2,2) + (-Power(m2,2) + s12)*s23) +
                    (2*Power(m1,4) + Power(m1,2)*(7*Power(m2,2) - 4*s12 - s23) + (Power(m2,2) - s12)*(Power(m2,2) - 2*s12 - s23))*Power(M(eps),2)) -
       2*F23*F6*((Power(m1,2) - s12 - s23)*(Power(m1,2)*Power(m2,2) + (-Power(m2,2) + s12)*s23) +
                 (2*Power(m1,4) + Power(m1,2)*(7*Power(m2,2) - 4*s12 - s23) + (Power(m2,2) - s12)*(Power(m2,2) - 2*s12 - s23))*Power(M(eps),2)) -
       Power(F6,2)*((Power(m1,2) - s12 - s23)*(Power(m1,2)*Power(m2,2) + (-Power(m2,2) + s12)*s23) +
                    (2*Power(m1,4) + Power(m1,2)*(7*Power(m2,2) - 4*s12 - s23) + (Power(m2,2) - s12)*(Power(m2,2) - 2*s12 - s23))*Power(M(eps),2)) -
       (F5*(F23 + F6)*(Power(m1,2) - s23 + Power(M(eps),2))*
        (Power(m2,2)*(Power(m1,2) - s23)*(Power(m1,2) - s12 - s23) +
         (2*Power(m1,4) + 4*Power(m2,4) + 2*s12*(s12 + s23) - Power(m2,2)*(7*s12 + 5*s23) + Power(m1,2)*(9*Power(m2,2) - 2*(2*s12 + s23)))*
         Power(M(eps),2) + 2*Power(m2,2)*Power(M(eps),4)))/Power(M(eps),2) +
       2*F4*(F23*Power(m1,2)*Power(-Power(m1,2) + s12 + s23,2) - 4*F23*Power(m1,2)*Power(m2,2)*Power(M(eps),2) -
             F23*Power(Power(m1,2) + Power(m2,2) - s12,2)*Power(M(eps),2) -
             F23*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) - s12 - s23)*(Power(m1,2) - s23 + Power(M(eps),2)) +
             F6*(-((Power(m1,2) - s12 - s23)*(Power(m1,2)*Power(m2,2) + (-Power(m2,2) + s12)*s23)) -
                 (2*Power(m1,4) + Power(m1,2)*(7*Power(m2,2) - 4*s12 - s23) + (Power(m2,2) - s12)*(Power(m2,2) - 2*s12 - s23))*Power(M(eps),2)) -
             (F5*(Power(m1,2) - s23 + Power(M(eps),2))*(Power(m2,2)*(Power(m1,2) - s23)*(Power(m1,2) - s12 - s23) +
                                                        (2*Power(m1,4) + 4*Power(m2,4) + 2*s12*(s12 + s23) - Power(m2,2)*(7*s12 + 5*s23) +
                                                         Power(m1,2)*(9*Power(m2,2) - 2*(2*s12 + s23)))*Power(M(eps),2) + 2*Power(m2,2)*Power(M(eps),4)))/(2.*Power(M(eps),2))) +
       4*F1*(-(F23*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) - s12 - s23)) - F23*Power(m2,2)*(Power(m1,2) - s23 + Power(M(eps),2)) -
             (F5*(Power(m2,2)*(Power(m1,2) - s23)*(Power(m1,2) - s12 - s23) +
                  (2*Power(m1,4) + 4*Power(m2,4) + 2*s12*(s12 + s23) - Power(m2,2)*(7*s12 + 5*s23) +
                   Power(m1,2)*(9*Power(m2,2) - 2*(2*s12 + s23)))*Power(M(eps),2) + 2*Power(m2,2)*Power(M(eps),4)))/(2.*Power(M(eps),2)) +
             F4*(-((Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) - s12 - s23)) - Power(m2,2)*(Power(m1,2) - s23 + Power(M(eps),2))) +
             F6*(-((Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) - s12 - s23)) - Power(m2,2)*(Power(m1,2) - s23 + Power(M(eps),2)))))/4.
      //------------------------------------------------------------------------------------------------------------------------------------------------
      ;
      }
      
      
      /*********************************************************************/
      /*                                                                   */
      /*        Photonmasse ungleich Null                                  */
      /*                                                                   */
      /*********************************************************************/

      if (mGam != 0.) {
      minv =
      
      //------------------------------------------------------------------------------------------------------------------------------------------------
      (-16*Power(F23,2)*Power(m1,2)*Power(m2,2)*Power(M(eps),2) - 4*Power(F23,2)*Power(Power(m1,2) + Power(m2,2) - s12,2)*Power(M(eps),2) +
       (4*Power(F23,2)*Power(m1,2)*Power(m2,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2))/Power(m3,2) +
       (Power(F23,2)*Power(Power(m1,2) + Power(m2,2) - s12,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2))/Power(m3,2) +
       4*Power(F23,2)*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) - s23 + Power(M(eps),2))*(Power(m2,2) - s31 + Power(M(eps),2)) -
       (Power(F23,2)*(Power(m1,2) + Power(m2,2) - s12)*Power(Power(m3,2) - s12 + Power(M(eps),2),2)*(Power(m1,2) - s23 + Power(M(eps),2))*
        (Power(m2,2) - s31 + Power(M(eps),2)))/(Power(m3,2)*Power(M(eps),2)) +
       4*Power(F23,2)*Power(m1,2)*Power(Power(m2,2) - s31 + Power(M(eps),2),2) -
       (Power(F23,2)*Power(m1,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2)*Power(Power(m2,2) - s31 + Power(M(eps),2),2))/
       (Power(m3,2)*Power(M(eps),2)) + (Power(F5,2)*(Power(m1,4) + Power(Power(m3,2) - s31,2) - 2*Power(m1,2)*(Power(m3,2) + s31))*
                                        (Power(Power(m2,3) - m2*s31,2) + 2*s31*(3*Power(m2,2) + s31)*Power(M(eps),2) - (3*Power(m2,2) + 4*s31)*Power(M(eps),4) +
                                         2*Power(M(eps),6)))/(Power(m3,2)*Power(M(eps),2)) +
       4*Power(F1,2)*(-16*Power(m2,2) + Power(Power(m2,2) + Power(m3,2) - s23,2)/Power(m3,2) -
                      ((Power(m2,2) + Power(m3,2) - s23)*(Power(m3,2) - s12 + Power(M(eps),2))*(Power(m2,2) - s31 + Power(M(eps),2)))/
                      (Power(m3,2)*Power(M(eps),2)) + Power(Power(m2,2) - s31 + Power(M(eps),2),2)/Power(M(eps),2)) +
       Power(F4,2)*(-16*Power(m1,2)*Power(m2,2)*Power(M(eps),2) - 4*Power(Power(m1,2) + Power(m2,2) - s12,2)*Power(M(eps),2) +
         (4*Power(m1,2)*Power(m2,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2))/Power(m3,2) +
         (Power(Power(m1,2) + Power(m2,2) - s12,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2))/Power(m3,2) +
         4*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) - s23 + Power(M(eps),2))*(Power(m2,2) - s31 + Power(M(eps),2)) -
         ((Power(m1,2) + Power(m2,2) - s12)*Power(Power(m3,2) - s12 + Power(M(eps),2),2)*(Power(m1,2) - s23 + Power(M(eps),2))*
          (Power(m2,2) - s31 + Power(M(eps),2)))/(Power(m3,2)*Power(M(eps),2)) +
         4*Power(m1,2)*Power(Power(m2,2) - s31 + Power(M(eps),2),2) -
         (Power(m1,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2)*Power(Power(m2,2) - s31 + Power(M(eps),2),2))/
         (Power(m3,2)*Power(M(eps),2))) + 2*F23*F6*(-16*Power(m1,2)*Power(m2,2)*Power(M(eps),2) -
         4*Power(Power(m1,2) + Power(m2,2) - s12,2)*Power(M(eps),2) +
         (4*Power(m1,2)*Power(m2,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2))/Power(m3,2) +
         (Power(Power(m1,2) + Power(m2,2) - s12,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2))/Power(m3,2) +
         4*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) - s23 + Power(M(eps),2))*(Power(m2,2) - s31 + Power(M(eps),2)) -
         ((Power(m1,2) + Power(m2,2) - s12)*Power(Power(m3,2) - s12 + Power(M(eps),2),2)*(Power(m1,2) - s23 + Power(M(eps),2))*
          (Power(m2,2) - s31 + Power(M(eps),2)))/(Power(m3,2)*Power(M(eps),2)) +
         4*Power(m1,2)*Power(Power(m2,2) - s31 + Power(M(eps),2),2) -
         (Power(m1,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2)*Power(Power(m2,2) - s31 + Power(M(eps),2),2))/
         (Power(m3,2)*Power(M(eps),2))) + Power(F6,2)*(-16*Power(m1,2)*Power(m2,2)*Power(M(eps),2) -
         4*Power(Power(m1,2) + Power(m2,2) - s12,2)*Power(M(eps),2) +
         (4*Power(m1,2)*Power(m2,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2))/Power(m3,2) +
         (Power(Power(m1,2) + Power(m2,2) - s12,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2))/Power(m3,2) +
         4*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) - s23 + Power(M(eps),2))*(Power(m2,2) - s31 + Power(M(eps),2)) -
         ((Power(m1,2) + Power(m2,2) - s12)*Power(Power(m3,2) - s12 + Power(M(eps),2),2)*(Power(m1,2) - s23 + Power(M(eps),2))*
          (Power(m2,2) - s31 + Power(M(eps),2)))/(Power(m3,2)*Power(M(eps),2)) +
         4*Power(m1,2)*Power(Power(m2,2) - s31 + Power(M(eps),2),2) -
         (Power(m1,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2)*Power(Power(m2,2) - s31 + Power(M(eps),2),2))/
         (Power(m3,2)*Power(M(eps),2))) + 2*F5*((-4*F23*Power(m2,2)*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) + Power(m3,2) - s31)*
         (Power(m3,2) - s12 + Power(M(eps),2)))/Power(m3,2) -
        8*F23*Power(m2,2)*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) - s23 + Power(M(eps),2)) -
        (2*F23*Power(m2,2)*(Power(m1,2) + Power(m3,2) - s31)*(Power(m3,2) - s12 + Power(M(eps),2))*(Power(m1,2) - s23 + Power(M(eps),2)))/
        Power(m3,2) - 4*F23*Power(m2,2)*Power(Power(m1,2) - s23 + Power(M(eps),2),2) +
        (2*F23*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) + Power(m3,2) - s31)*(Power(m3,2) - s12 + Power(M(eps),2))*
         (Power(m2,2) - s31 + Power(M(eps),2)))/Power(m3,2) +
        4*F23*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) - s23 + Power(M(eps),2))*(Power(m2,2) - s31 + Power(M(eps),2)) +
        (F23*Power(m2,2)*(Power(m1,2) + Power(m3,2) - s31)*(Power(m3,2) - s12 + Power(M(eps),2))*(Power(m1,2) - s23 + Power(M(eps),2))*
         (Power(m2,2) - s31 + Power(M(eps),2)))/(Power(m3,2)*Power(M(eps),2)) +
        (2*F23*Power(m2,2)*Power(Power(m1,2) - s23 + Power(M(eps),2),2)*(Power(m2,2) - s31 + Power(M(eps),2)))/Power(M(eps),2) +
        F6*((-4*Power(m2,2)*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) + Power(m3,2) - s31)*(Power(m3,2) - s12 + Power(M(eps),2)))/
            Power(m3,2) - 8*Power(m2,2)*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) - s23 + Power(M(eps),2)) -
            (2*Power(m2,2)*(Power(m1,2) + Power(m3,2) - s31)*(Power(m3,2) - s12 + Power(M(eps),2))*(Power(m1,2) - s23 + Power(M(eps),2)))/
            Power(m3,2) - 4*Power(m2,2)*Power(Power(m1,2) - s23 + Power(M(eps),2),2) +
            (2*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) + Power(m3,2) - s31)*(Power(m3,2) - s12 + Power(M(eps),2))*
             (Power(m2,2) - s31 + Power(M(eps),2)))/Power(m3,2) +
            4*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) - s23 + Power(M(eps),2))*(Power(m2,2) - s31 + Power(M(eps),2)) +
            (Power(m2,2)*(Power(m1,2) + Power(m3,2) - s31)*(Power(m3,2) - s12 + Power(M(eps),2))*(Power(m1,2) - s23 + Power(M(eps),2))*
             (Power(m2,2) - s31 + Power(M(eps),2)))/(Power(m3,2)*Power(M(eps),2)) +
            (2*Power(m2,2)*Power(Power(m1,2) - s23 + Power(M(eps),2),2)*(Power(m2,2) - s31 + Power(M(eps),2)))/Power(M(eps),2))) +
       2*F1*((2*F23*(Power(m1,2) + Power(m2,2) - s12)*(Power(m2,2) + Power(m3,2) - s23)*(Power(m3,2) - s12 + Power(M(eps),2)))/Power(m3,2) -
             (4*F23*Power(m2,2)*(Power(m1,2) + Power(m3,2) - s31)*(Power(m3,2) - s12 + Power(M(eps),2)))/Power(m3,2) -
             8*F23*Power(m2,2)*(Power(m1,2) - s23 + Power(M(eps),2)) +
             8*F23*(Power(m1,2) + Power(m2,2) - s12)*(Power(m2,2) - s31 + Power(M(eps),2)) +
             (F23*(Power(m2,2) + Power(m3,2) - s23)*(-Power(m1,2) + s23 - Power(M(eps),2))*(Power(m3,2) - s12 + Power(M(eps),2))*
              (Power(m2,2) - s31 + Power(M(eps),2)))/(Power(m3,2)*Power(M(eps),2)) -
             (F23*(Power(m1,2) + Power(m2,2) - s12)*Power(Power(m3,2) - s12 + Power(M(eps),2),2)*(Power(m2,2) - s31 + Power(M(eps),2)))/
             (Power(m3,2)*Power(M(eps),2)) + (F23*(Power(m1,2) + Power(m3,2) - s31)*(Power(m3,2) - s12 + Power(M(eps),2))*
                                              Power(Power(m2,2) - s31 + Power(M(eps),2),2))/(Power(m3,2)*Power(M(eps),2)) +
             (2*F5*(Power(m2,2)*(Power(m2,2) - s31)*(Power(m3,4) + Power(m1,2)*(3*Power(m3,2) - s12) + s12*s31 -
                                                     Power(m3,2)*(s12 + 2*s23 + s31)) + (Power(m2,4)*(-3*Power(m3,2) + s31) +
                2*s31*(-Power(m3,4) - s23*s31 + Power(m3,2)*(2*s12 + s23 + s31)) -
                Power(m1,2)*(Power(m2,4) + 6*Power(m3,2)*s31 - 2*s23*s31 + Power(m2,2)*(9*Power(m3,2) - s12 - 2*s23 + 3*s31)) -
                Power(m2,2)*(3*Power(m3,4) + (s12 + 2*s23 - 3*s31)*s31 + Power(m3,2)*(-5*s12 - 4*s23 + 6*s31)))*Power(M(eps),2) +
                    (Power(m1,2)*(Power(m2,2) + 6*Power(m3,2) - 2*s23) + Power(m2,2)*(3*Power(m3,2) - s31) +
                     2*(Power(m3,4) + s23*s31 - Power(m3,2)*(2*s12 + s23 + s31)))*Power(M(eps),4)))/(Power(m3,2)*Power(M(eps),2)) +
             F4*((2*(Power(m1,2) + Power(m2,2) - s12)*(Power(m2,2) + Power(m3,2) - s23)*(Power(m3,2) - s12 + Power(M(eps),2)))/Power(m3,2) -
                 (4*Power(m2,2)*(Power(m1,2) + Power(m3,2) - s31)*(Power(m3,2) - s12 + Power(M(eps),2)))/Power(m3,2) -
                 8*Power(m2,2)*(Power(m1,2) - s23 + Power(M(eps),2)) + 
                 8*(Power(m1,2) + Power(m2,2) - s12)*(Power(m2,2) - s31 + Power(M(eps),2)) + 
                 ((Power(m2,2) + Power(m3,2) - s23)*(-Power(m1,2) + s23 - Power(M(eps),2))*(Power(m3,2) - s12 + Power(M(eps),2))*
                  (Power(m2,2) - s31 + Power(M(eps),2)))/(Power(m3,2)*Power(M(eps),2)) - 
                 ((Power(m1,2) + Power(m2,2) - s12)*Power(Power(m3,2) - s12 + Power(M(eps),2),2)*(Power(m2,2) - s31 + Power(M(eps),2)))/
                 (Power(m3,2)*Power(M(eps),2)) + ((Power(m1,2) + Power(m3,2) - s31)*(Power(m3,2) - s12 + Power(M(eps),2))*
                                                  Power(Power(m2,2) - s31 + Power(M(eps),2),2))/(Power(m3,2)*Power(M(eps),2))) + 
             F6*((2*(Power(m1,2) + Power(m2,2) - s12)*(Power(m2,2) + Power(m3,2) - s23)*(Power(m3,2) - s12 + Power(M(eps),2)))/Power(m3,2) - 
                 (4*Power(m2,2)*(Power(m1,2) + Power(m3,2) - s31)*(Power(m3,2) - s12 + Power(M(eps),2)))/Power(m3,2) - 
                 8*Power(m2,2)*(Power(m1,2) - s23 + Power(M(eps),2)) + 
                 8*(Power(m1,2) + Power(m2,2) - s12)*(Power(m2,2) - s31 + Power(M(eps),2)) + 
                 ((Power(m2,2) + Power(m3,2) - s23)*(-Power(m1,2) + s23 - Power(M(eps),2))*(Power(m3,2) - s12 + Power(M(eps),2))*
                  (Power(m2,2) - s31 + Power(M(eps),2)))/(Power(m3,2)*Power(M(eps),2)) - 
                 ((Power(m1,2) + Power(m2,2) - s12)*Power(Power(m3,2) - s12 + Power(M(eps),2),2)*(Power(m2,2) - s31 + Power(M(eps),2)))/
                 (Power(m3,2)*Power(M(eps),2)) + ((Power(m1,2) + Power(m3,2) - s31)*(Power(m3,2) - s12 + Power(M(eps),2))*
                                                  Power(Power(m2,2) - s31 + Power(M(eps),2),2))/(Power(m3,2)*Power(M(eps),2)))) + 
       2*F4*(-16*F23*Power(m1,2)*Power(m2,2)*Power(M(eps),2) - 4*F23*Power(Power(m1,2) + Power(m2,2) - s12,2)*Power(M(eps),2) + 
             (4*F23*Power(m1,2)*Power(m2,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2))/Power(m3,2) + 
             (F23*Power(Power(m1,2) + Power(m2,2) - s12,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2))/Power(m3,2) + 
             4*F23*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) - s23 + Power(M(eps),2))*(Power(m2,2) - s31 + Power(M(eps),2)) - 
             (F23*(Power(m1,2) + Power(m2,2) - s12)*Power(Power(m3,2) - s12 + Power(M(eps),2),2)*(Power(m1,2) - s23 + Power(M(eps),2))*
              (Power(m2,2) - s31 + Power(M(eps),2)))/(Power(m3,2)*Power(M(eps),2)) + 
             4*F23*Power(m1,2)*Power(Power(m2,2) - s31 + Power(M(eps),2),2) - 
             (F23*Power(m1,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2)*Power(Power(m2,2) - s31 + Power(M(eps),2),2))/
             (Power(m3,2)*Power(M(eps),2)) + F5*((-4*Power(m2,2)*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) + Power(m3,2) - s31)*
              (Power(m3,2) - s12 + Power(M(eps),2)))/Power(m3,2) - 
             8*Power(m2,2)*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) - s23 + Power(M(eps),2)) - 
             (2*Power(m2,2)*(Power(m1,2) + Power(m3,2) - s31)*(Power(m3,2) - s12 + Power(M(eps),2))*(Power(m1,2) - s23 + Power(M(eps),2)))/
             Power(m3,2) - 4*Power(m2,2)*Power(Power(m1,2) - s23 + Power(M(eps),2),2) + 
             (2*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) + Power(m3,2) - s31)*(Power(m3,2) - s12 + Power(M(eps),2))*
              (Power(m2,2) - s31 + Power(M(eps),2)))/Power(m3,2) + 
             4*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) - s23 + Power(M(eps),2))*(Power(m2,2) - s31 + Power(M(eps),2)) + 
             (Power(m2,2)*(Power(m1,2) + Power(m3,2) - s31)*(Power(m3,2) - s12 + Power(M(eps),2))*(Power(m1,2) - s23 + Power(M(eps),2))*
              (Power(m2,2) - s31 + Power(M(eps),2)))/(Power(m3,2)*Power(M(eps),2)) + 
             (2*Power(m2,2)*Power(Power(m1,2) - s23 + Power(M(eps),2),2)*(Power(m2,2) - s31 + Power(M(eps),2)))/Power(M(eps),2)) + 
             F6*(-16*Power(m1,2)*Power(m2,2)*Power(M(eps),2) - 4*Power(Power(m1,2) + Power(m2,2) - s12,2)*Power(M(eps),2) + 
                 (4*Power(m1,2)*Power(m2,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2))/Power(m3,2) + 
                 (Power(Power(m1,2) + Power(m2,2) - s12,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2))/Power(m3,2) + 
                 4*(Power(m1,2) + Power(m2,2) - s12)*(Power(m1,2) - s23 + Power(M(eps),2))*(Power(m2,2) - s31 + Power(M(eps),2)) - 
                 ((Power(m1,2) + Power(m2,2) - s12)*Power(Power(m3,2) - s12 + Power(M(eps),2),2)*(Power(m1,2) - s23 + Power(M(eps),2))*
                  (Power(m2,2) - s31 + Power(M(eps),2)))/(Power(m3,2)*Power(M(eps),2)) + 
                 4*Power(m1,2)*Power(Power(m2,2) - s31 + Power(M(eps),2),2) - 
                 (Power(m1,2)*Power(Power(m3,2) - s12 + Power(M(eps),2),2)*Power(Power(m2,2) - s31 + Power(M(eps),2),2))/
                 (Power(m3,2)*Power(M(eps),2)))))/16.
      //------------------------------------------------------------------------------------------------------------------------------------------------
      ;
      }
      
      //DEBUG
      //minv = 1;
      
      return minv;
}
