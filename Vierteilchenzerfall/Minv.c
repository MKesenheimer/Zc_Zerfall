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

double getMinv(double eps, double lam, double s12, double s34, double xsi, double tht, double phi, double F1, double F23, double F4, double F5, double F6) {

      double minv;
      
      //Funktionen:
      double sig12 = ksqrt(s12,Power(m1,2),Power(m2,2))/s12;
      double sig34 = ksqrt(s34,Power(m3,2),Power(m4,2))/s34;
      double chi = ksqrt(Power(M(eps),2),s12,s34)/2.;
      double KK = s12;
      double QQ = s34;
      double KQ = (Power(M(eps),2) - s12 - s34)/2.;
      double KL = Power(m1,2) - Power(m2,2);
      double QR = Power(m3,2) - Power(m4,2);
      double LQ = KL*KQ/KK+chi*sig12*Cos(xsi);
      double KR = QR*KQ/QQ+chi*sig34*Cos(tht);
      double LR = (KL*chi*sig34*Cos(tht))/KK + (QR*chi*sig12*Cos(xsi))/QQ +
                  KQ*sig12*sig34*Cos(tht)*Cos(xsi) + (KL*KQ*QR)/(KK*QQ) -
                  Sqrt(s12*s34)*sig12*sig34*Cos(phi)*Sin(tht)*Sin(xsi);
      double s134 = KQ + LQ + Power(m1,2) + s34;
      double F2346 = F23 + F4 + F6; //Zusammengefasster Formfaktor, da die Diagramme F23, F4 und F6 die gleiche Struktur besitzen
      
      //Matrixelement
      minv =
      //----------------------------------------------------------------------------------------------------------------------------------------------------------------
      
      (Power(q,2)*(-8*Power(F2346,2)*(KK - Power(m1,2) - Power(m2,2))*(KK + KQ + LQ + Power(m1,2) - Power(m2,2))*
                   (KK + KQ - LQ - Power(m1,2) + Power(m2,2))*Power(ml,2) +
                   8*Power(F2346,2)*Power(m1,2)*Power(KK + KQ - LQ - Power(m1,2) + Power(m2,2),2)*Power(ml,2) +
                   4*Power(F2346,2)*(KK - Power(m1,2) - Power(m2,2))*(KK + KQ + LQ + Power(m1,2) - Power(m2,2))*
                   (KK + KQ - LQ - Power(m1,2) + Power(m2,2))*(Power(m3,2) + Power(m4,2) - QQ) -
                   4*Power(F2346,2)*Power(m1,2)*Power(KK + KQ - LQ - Power(m1,2) + Power(m2,2),2)*
                   (Power(m3,2) + Power(m4,2) - QQ) - 32*Power(F2346,2)*Power(m1,2)*Power(m2,2)*Power(ml,2)*
                   (KK + 2*KQ + QQ) - 8*Power(F2346,2)*Power(-KK + Power(m1,2) + Power(m2,2),2)*Power(ml,2)*
                   (KK + 2*KQ + QQ) - 16*Power(F2346,2)*Power(m1,2)*Power(m2,2)*(KK + 2*KQ + QQ)*
                   (-Power(m3,2) - Power(m4,2) + QQ) - 4*Power(F2346,2)*Power(-KK + Power(m1,2) + Power(m2,2),2)*
                   (KK + 2*KQ + QQ)*(-Power(m3,2) - Power(m4,2) + QQ) +
                   16*Power(F2346,2)*Power(m1,2)*Power(m2,2)*(KQ + KR + Power(m3,2) - Power(m4,2) + QQ)*
                   (KQ - KR - Power(m3,2) + Power(m4,2) + QQ) +
                   4*Power(F2346,2)*Power(-KK + Power(m1,2) + Power(m2,2),2)*(KQ + KR + Power(m3,2) - Power(m4,2) + QQ)*
                   (KQ - KR - Power(m3,2) + Power(m4,2) + QQ) +
                   (4*Power(F2346,2)*(KK - Power(m1,2) - Power(m2,2))*(KK + KQ + LQ + Power(m1,2) - Power(m2,2))*
                    (KK + KQ - LQ - Power(m1,2) + Power(m2,2))*(KQ + KR + Power(m3,2) - Power(m4,2) + QQ)*
                    (KQ - KR - Power(m3,2) + Power(m4,2) + QQ))/(KK + 2*KQ + QQ) -
                   (4*Power(F2346,2)*Power(m1,2)*Power(KK + KQ - LQ - Power(m1,2) + Power(m2,2),2)*
                    (KQ + KR + Power(m3,2) - Power(m4,2) + QQ)*(KQ - KR - Power(m3,2) + Power(m4,2) + QQ))/
                   (KK + 2*KQ + QQ) + 4*Power(F1,2)*((KQ + KR - LQ - LR)*(KQ - KR - LQ + LR) -
                                                     40*Power(m2,2)*Power(ml,2) + 12*Power(m2,2)*(Power(m3,2) + Power(m4,2) - QQ) +
                                                     (4*Power(KK + KQ - LQ - Power(m1,2) + Power(m2,2),2)*Power(ml,2))/(KK + 2*KQ + QQ) +
                                                     ((KQ - KR - LQ + LR)*(KK + KQ - LQ - Power(m1,2) + Power(m2,2))*
                                                      (KQ + KR + Power(m3,2) - Power(m4,2) + QQ))/(KK + 2*KQ + QQ) +
                                                     ((KQ + KR - LQ - LR)*(KK + KQ - LQ - Power(m1,2) + Power(m2,2))*
                                                      (KQ - KR - Power(m3,2) + Power(m4,2) + QQ))/(KK + 2*KQ + QQ)) +
                   2*F1*(16*F2346*Power(m2,2)*(-KK - KQ - LQ - Power(m1,2) + Power(m2,2))*Power(ml,2) -
                         16*F2346*(KK - Power(m1,2) - Power(m2,2))*(KK + KQ - LQ - Power(m1,2) + Power(m2,2))*Power(ml,2) -
                         8*F2346*Power(m2,2)*(-KK - KQ - LQ - Power(m1,2) + Power(m2,2))*(Power(m3,2) + Power(m4,2) - QQ) +
                         8*F2346*(KK - Power(m1,2) - Power(m2,2))*(KK + KQ - LQ - Power(m1,2) + Power(m2,2))*
                         (Power(m3,2) + Power(m4,2) - QQ) +
                         4*F2346*(KQ - KR + LQ - LR)*Power(m2,2)*(KQ + KR + Power(m3,2) - Power(m4,2) + QQ) +
                         2*F2346*(KQ - KR - LQ + LR)*(KK - Power(m1,2) - Power(m2,2))*
                         (KQ + KR + Power(m3,2) - Power(m4,2) + QQ) +
                         (F2346*(KQ - KR - LQ + LR)*(KK + KQ + LQ + Power(m1,2) - Power(m2,2))*
                          (KK + KQ - LQ - Power(m1,2) + Power(m2,2))*(KQ + KR + Power(m3,2) - Power(m4,2) + QQ))/
                         (KK + 2*KQ + QQ) - (F2346*(KQ - KR + LQ - LR)*Power(KK + KQ - LQ - Power(m1,2) + Power(m2,2),2)*
                                             (KQ + KR + Power(m3,2) - Power(m4,2) + QQ))/(KK + 2*KQ + QQ) +
                         4*F2346*(KQ + KR + LQ + LR)*Power(m2,2)*(KQ - KR - Power(m3,2) + Power(m4,2) + QQ) +
                         2*F2346*(KQ + KR - LQ - LR)*(KK - Power(m1,2) - Power(m2,2))*
                         (KQ - KR - Power(m3,2) + Power(m4,2) + QQ) +
                         (F2346*(KQ + KR - LQ - LR)*(KK + KQ + LQ + Power(m1,2) - Power(m2,2))*
                          (KK + KQ - LQ - Power(m1,2) + Power(m2,2))*(KQ - KR - Power(m3,2) + Power(m4,2) + QQ))/
                         (KK + 2*KQ + QQ) - (F2346*(KQ + KR + LQ + LR)*Power(KK + KQ - LQ - Power(m1,2) + Power(m2,2),2)*
                                             (KQ - KR - Power(m3,2) + Power(m4,2) + QQ))/(KK + 2*KQ + QQ) +
                         (4*F2346*(KK - Power(m1,2) - Power(m2,2))*(KK + KQ - LQ - Power(m1,2) + Power(m2,2))*
                          (KQ + KR + Power(m3,2) - Power(m4,2) + QQ)*(KQ - KR - Power(m3,2) + Power(m4,2) + QQ))/
                         (KK + 2*KQ + QQ) + (2*F5*(-2*(KQ + KR - LQ - LR)*(KQ - KR + LQ - LR)*Power(m2,2) -
                                                   2*(KQ - KR - LQ + LR)*(KQ + KR + LQ + LR)*Power(m2,2) +
                                                   (KQ + KR - LQ - LR)*(KQ - KR + LQ - LR)*(KK + KQ - LQ - Power(m1,2) + Power(m2,2)) +
                                                   (KQ - KR - LQ + LR)*(KQ + KR + LQ + LR)*(KK + KQ - LQ - Power(m1,2) + Power(m2,2)) +
                                                   8*Power(m2,2)*(-KK - KQ - LQ - Power(m1,2) + Power(m2,2))*Power(ml,2) -
                                                   8*(KK - Power(m1,2) - Power(m2,2))*(KK + KQ - LQ - Power(m1,2) + Power(m2,2))*Power(ml,2) -
                                                   16*Power(m2,2)*(-KK + Power(m1,2) + Power(m2,2))*Power(ml,2) -
                                                   4*Power(m2,2)*(-KK - KQ - LQ - Power(m1,2) + Power(m2,2))*(Power(m3,2) + Power(m4,2) - QQ) +
                                                   4*(KK - Power(m1,2) - Power(m2,2))*(KK + KQ - LQ - Power(m1,2) + Power(m2,2))*
                                                   (Power(m3,2) + Power(m4,2) - QQ) +
                                                   8*Power(m2,2)*(-KK + Power(m1,2) + Power(m2,2))*(Power(m3,2) + Power(m4,2) - QQ) -
                                                   (4*Power(m2,2)*(-KK - KQ - LQ - Power(m1,2) + Power(m2,2))*
                                                    (KK + KQ - LQ - Power(m1,2) + Power(m2,2))*Power(ml,2))/(KK + 2*KQ + QQ) +
                                                   (2*Power(m2,2)*(-KK - KQ - LQ - Power(m1,2) + Power(m2,2))*
                                                    (KK + KQ - LQ - Power(m1,2) + Power(m2,2))*(Power(m3,2) + Power(m4,2) - QQ))/
                                                   (KK + 2*KQ + QQ) + 2*(KQ - KR + LQ - LR)*Power(m2,2)*
                                                   (KQ + KR + Power(m3,2) - Power(m4,2) + QQ) -
                                                   ((KQ - KR + LQ - LR)*Power(m2,2)*(KK + KQ - LQ - Power(m1,2) + Power(m2,2))*
                                                    (KQ + KR + Power(m3,2) - Power(m4,2) + QQ))/(KK + 2*KQ + QQ) +
                                                   2*(KQ + KR + LQ + LR)*Power(m2,2)*(KQ - KR - Power(m3,2) + Power(m4,2) + QQ) -
                                                   ((KQ + KR + LQ + LR)*Power(m2,2)*(KK + KQ - LQ - Power(m1,2) + Power(m2,2))*
                                                    (KQ - KR - Power(m3,2) + Power(m4,2) + QQ))/(KK + 2*KQ + QQ)))/(Power(m1,2) - s134)) +
                   (Power(F5,2)*(Power(KQ,2) - Power(KR,2) + 2*KQ*LQ + Power(LQ,2) - 2*KR*LR - Power(LR,2) +
                                 4*Power(m1,2)*Power(m3,2) + 4*Power(m1,2)*Power(m4,2) - 8*Power(m1,2)*Power(ml,2) -
                                 4*Power(m1,2)*QQ)*(2*Power(KK,3) + 4*Power(KQ,3) + Power(LQ,2)*Power(m2,2) +
                                                    2*LQ*Power(m1,2)*Power(m2,2) + Power(m1,4)*Power(m2,2) - 2*LQ*Power(m2,4) -
                                                    2*Power(m1,2)*Power(m2,4) + Power(m2,6) + 2*Power(LQ,2)*QQ + 4*LQ*Power(m1,2)*QQ +
                                                    2*Power(m1,4)*QQ + 8*LQ*Power(m2,2)*QQ + 8*Power(m1,2)*Power(m2,2)*QQ - 2*Power(m2,4)*QQ +
                                                    4*Power(m2,2)*Power(QQ,2) + Power(KK,2)*(8*KQ - 4*LQ - 4*Power(m1,2) - 3*Power(m2,2) + 2*QQ) +
                                                    2*KK*(KQ - LQ - Power(m1,2))*(5*KQ - LQ - Power(m1,2) - 3*Power(m2,2) + 2*QQ) +
                                                    Power(KQ,2)*(-8*LQ - 8*Power(m1,2) + Power(m2,2) + 2*QQ) + 
                                                    2*KQ*(2*Power(LQ,2) + 2*Power(m1,4) - Power(m2,4) + Power(m1,2)*(7*Power(m2,2) - 2*QQ) + 
                                                          LQ*(4*Power(m1,2) + 7*Power(m2,2) - 2*QQ) + 4*Power(m2,2)*QQ)))/
                   ((KK + 2*KQ + QQ)*Power(Power(m1,2) - s134,2)) + 
                   (4*F2346*F5*(Power(KQ,2) - Power(KR,2) + KK*Power(m3,2) + LQ*Power(m3,2) - LR*Power(m3,2) + 
                                Power(m1,2)*Power(m3,2) - Power(m2,2)*Power(m3,2) + KK*Power(m4,2) + LQ*Power(m4,2) + 
                                LR*Power(m4,2) + Power(m1,2)*Power(m4,2) - Power(m2,2)*Power(m4,2) - 
                                KR*(LR + Power(m3,2) - Power(m4,2)) - 2*KK*Power(ml,2) - 2*LQ*Power(ml,2) - 
                                2*Power(m1,2)*Power(ml,2) + 2*Power(m2,2)*Power(ml,2) + 
                                KQ*(LQ + Power(m3,2) + Power(m4,2) - 2*Power(ml,2)) - KK*QQ - Power(m1,2)*QQ + Power(m2,2)*QQ)*
                    (2*Power(KK,3) + Power(LQ,2)*Power(m2,2) + 2*LQ*Power(m1,2)*Power(m2,2) + Power(m1,4)*Power(m2,2) - 
                     2*LQ*Power(m2,4) - 2*Power(m1,2)*Power(m2,4) + Power(m2,6) - 
                     Power(KQ,2)*(4*Power(m1,2) + Power(m2,2)) + 
                     2*KQ*(2*LQ*(Power(m1,2) + 2*Power(m2,2)) + Power(m1,2)*(2*Power(m1,2) + 6*Power(m2,2) - QQ)) + 
                     2*LQ*Power(m1,2)*QQ + 2*Power(m1,4)*QQ + 4*LQ*Power(m2,2)*QQ + 6*Power(m1,2)*Power(m2,2)*QQ + 
                     Power(KK,2)*(6*KQ - 2*LQ - 4*Power(m1,2) - 3*Power(m2,2) + 2*QQ) + 
                     2*KK*(2*Power(KQ,2) + Power(m1,4) + 3*Power(m1,2)*Power(m2,2) + 
                           LQ*(Power(m1,2) + 2*Power(m2,2) - QQ) - 2*Power(m1,2)*QQ - Power(m2,2)*QQ + 
                           KQ*(-2*LQ - 5*Power(m1,2) - 3*Power(m2,2) + QQ))))/((KK + 2*KQ + QQ)*(Power(m1,2) - s134))))/
      (8.*Power(s34,2))
      
      //----------------------------------------------------------------------------------------------------------------------------------------------------------------
      ;
      
      return minv;
}
