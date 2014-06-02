/* Kompilieren mit g++ -o PS_Zc PS_Zc.cpp -lcuba -lm -lstdc++ -lncurses -O3 */

/*
 * Programm um die einfach differentielle Zerfallsrate dG/ds34 in Abhängigkeit
 * der Bindungsenergie Epsilon und des Größenparameters Epsilon
 * zu berechen
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cuba.h"
#include <time.h>
#include <string.h>
#include <string>
#include <unistd.h>
#include <curses.h>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

//Für welches System wird kompiliert?
#define MAC_OSX
//#define LINUX

/*********************************************************************/
/*                                                                   */
/*        Variablen Definitionen, alle Größen in GeV                 */
/*                                                                   */
/*********************************************************************/
//Parameter
#define lambda 0.5
#define lambdaStep 0.05
#define epsilon 0.0005
#define epsilonStep 0.0005

//Massen
#define mPi 0.1396
#define mJ 3.0969
#define mGam 0.
#define me 511.0e-6
#define mmu 105.6e-3
#define mD 1.865
#define mDs 2.010
//Die Zc Masse ist von eps abhängig

//Zerfallsraten
#define Gam1 2.5283751359093973e-17 //Zerfallsrate des Pions
#define GamM 7e-3 //Zerfallsrate Z_c (Schätzwert)

//Zuordnungen
//#define M mZ
#define m1 mPi
#define m2 mJ
#define m3 me
#define m4 me
#define ml me
#define MD mD
#define MDs mDs

//Kopplungskonstanten
#define gJDD 7.44 /* für JPsi in der 1s Welle*/
#define gDsDPi 17.9
#define gDDsJPi (gJDD*gDsDPi/(2.0*sqrt(2.0)))
#define alpFein (1.0/137.036)
#define q (sqrt(4*M_PI*alpFein))

//Sonstige Konstanten
#define Pi M_PI

/*********************************************************************/
/*                                                                   */
/*                Parameter für die Integration                      */
/*                                                                   */
/*********************************************************************/
//Feynmanintegration
#define NDIMF 3 //Dimensionalität des Hyperkubus, in dem die Integration stattfindet
#define NCOMPF 5 //Anzahl der Komponenten der Funktion f (Vektor)
#define USERDATAF NULL
#define EPSRELF 1e-3 //Aus EPSREL und EPSABS wird das Maximum gesucht: max(EPSREL * Integral, EPSABS)
#define EPSABSF 1e-12 //1e-12
#define LASTF 4
#define VERBOSEF 0
#define SEEDF 0
#define MINEVALF 0
#define MAXEVALF 1000000 //erst bei 2.5M Evaluationen wird eine relative Genauigkeit von 1e-2 erreicht! Mit 10000000 ist man sicher.
#define NSTARTF 500 //1000
#define NINCREASEF 500 //500
#define NBATCHF 1000
#define GRIDNOF 0 //
#define STATEFILEF NULL
#define NNEWF 1000
#define FLATNESSF 25.
#define KEYF 0
#define KEY1F 47
#define KEY2F 1
#define KEY3F 1
#define MAXPASSF 5
#define BORDERF 0.
#define MAXCHISQF 10.
#define MINDEVIATIONF .25
#define NGIVENF 0
#define LDXGIVENF NDIMF
#define NEXTRAF 0

//Phasenraumintegration
#define NDIMPS 4
#define NCOMPPS 1
#define USERDATAPS NULL
#define EPSRELPS 1e-3 //1e-3
#define EPSABSPS 1e-12 //1e-12
#define LASTPS 4
#define VERBOSEPS 0
#define SEEDPS 4
#define MINEVALPS 0
#define MAXEVALPS 50000 //50000
#define NSTARTPS 500
#define NINCREASEPS 500
#define NBATCHPS 1000
#define GRIDNOPS 0
#define STATEFILEPS NULL
#define NNEWPS 1000
#define FLATNESSPS 25.
#define KEYPS 0
#define KEY1PS 47
#define KEY2PS 1
#define KEY3PS 1
#define MAXPASSPS 5
#define BORDERPS 0.
#define MAXCHISQPS 10.
#define MINDEVIATIONPS .25
#define NGIVENPS 0
#define LDXGIVENPS NDIMPS
#define NEXTRAPS 0

//Wie viele Werte berechnen?
#define EPSCOUNT 1 //10  //Epsilon nach unten
#define LAMCOUNT 1 //7  //Lambda nach rechts
//Wo mit Rechnen anfangen?
#define STRTEPS (1-1)
#define STRTLAM (1-1)
//Maximale Anzahl der gZc-Werte aus Mathematica
#define MAXEPSCOUNT 10
#define MAXLAMCOUNT 7
//Schrittweite und Variablen für die partielle Zerfallsbreite
#define MAXNS 100

//aktuelle Ergebnisse
static double total[3 * NCOMPF + 3];
static int isnanFail = 0;

/*********************************************************************/
/*                                                                   */
/*                  Wahl der Integrationsmethode                     */
/*                                                                   */
/*********************************************************************/
#define VEGAS 0
#define SUAVE 1
#define DIVONNE 2
#define CUHRE 3

#ifndef INTMETH
#define INTMETH VEGAS //Gute Ergebnisse mit Suave und Vegas
#endif

//Optimierungen
//#define CUBACORES 1 //Benutze nur einen CPU-Kern
//#define OPTIMIZED //Programm mit minimaler IO Ausgabe ausführen, damit um nahezu 20% schneller

/*********************************************************************/
/*                                                                   */
/*                        Helper Functions                           */
/*                                                                   */
/*********************************************************************/
//Prototypen
double M(double eps);

double Power(double base, double exponent) {
      return pow(base, exponent);
}

double Sqrt(double arg) {
      return sqrt(arg);
}

double Cos(double arg) {
      return cos(arg);
}

double Sin(double arg) {
      return sin(arg);
}

//Betragsfunktion
double absolut(double x) {
      return sqrt(pow(x,2));
}

#ifdef LINUX
string to_string(int number){
      string number_string = "";
      char ones_char;
      int ones = 0;
      while(true){
            ones = number % 10;
            switch(ones){
                  case 0: ones_char = '0'; break;
                  case 1: ones_char = '1'; break;
                  case 2: ones_char = '2'; break;
                  case 3: ones_char = '3'; break;
                  case 4: ones_char = '4'; break;
                  case 5: ones_char = '5'; break;
                  case 6: ones_char = '6'; break;
                  case 7: ones_char = '7'; break;
                  case 8: ones_char = '8'; break;
                  case 9: ones_char = '9'; break;
            }
            number -= ones;
            number_string = ones_char + number_string;
            if(number == 0){
                  break;
            }
            number = number/10;
      }
      return number_string;
}
#endif

/*********************************************************************/
/*                                                                   */
/*                  Statusausgabe, Dateiausgabe                      */
/*                                                                   */
/*********************************************************************/
#define SCREENWIDTH 80
#define SCREENHEIGHT 40
static int pid = getpid();
static int childPID1 = 0;
static int childPID2 = 0;
//Dateiausgabe
ofstream oFF;
ofstream results;
ofstream dGamma;

//Integrationszyklen mitzählen um abschätzen zu können,
//wie lange die Rechnung dauert
static int cyclePS = 0;
static float percent = 0;
static int cycleTime = 0;

double getPercent(void){
      percent = (float)(cyclePS)/(EPSCOUNT*LAMCOUNT*MAXNS)*100;
      
      if (percent >= 100) {
            percent = 100;
      }
      return percent;
}

void getStatus(double s12, double s34, double xsi, double tht, double phi, double gZc, double g, double eps,
               double lam, double F1, double errF1, double F23, double errF23, double F4, double errF4,
               double F5, double errF5, double F6, double errF6, int nregionsf, int nevalf, int failf, double minv){
      clear();
      move(0,0);
      
      printw("%Master PID: %d\n",pid);
      
      //Werte ausgeben
      printw("cyclePS:\t%d\ns12:\t\t%f\t[GeV^2]\ns34:\t\t%f\t[GeV^2]\nxsi:\t\t%f\ntht:\t\t%f\nphi:\t\t%f\ngZc:\t\t%f\ng:\t\t%f\nmGam:\t\t%.0f\t\t[MeV]\n\nF1:\t\t%f +- %f\nF23:\t\t%f +- %f\n",
             cyclePS, s12, s34, xsi, tht, phi, gZc, g, mGam*1000, F1, errF1, F23, errF23);
      printw("F4:\t\t%f +- %f\nF5:\t\t%f +- %f\nF6:\t\t%f +- %f\n",
             F4, errF4, F5, errF5, F6, errF6);

      printw("nregionsf: %d\t nevalf: %d\t failf: %d\n\n",nregionsf, nevalf, failf);
      printw("Minv2:\t\t%f\n\n",minv);
      
      //Ladebalken
      move(20,0);
      printw("[");
      for (int n = 0; n <= (int)((SCREENWIDTH-10)*percent/100); n++) printw("#");
      move(20,SCREENWIDTH-9);
      printw("]%.1f%%", getPercent());
      refresh();
}

/*********************************************************************/
/*                                                                   */
/*        Parameter, die von anderen Größen abhängig sind            */
/*                                                                   */
/*********************************************************************/
//Källen Funktion
double k(double x, double y, double z){
      return (pow(x,2)+pow(y,2)+pow(z,2)-2*x*y-2*x*z-2*y*z);
}

//Wurzel aus Källen Funktion
double ksqrt(double x, double y, double z){
      
      //Numerischen Quatsch vermeiden:
      if (k(x,y,z) < 0) {
            return 0;
      }
      
      return (sqrt(pow(x,2)+pow(y,2)+pow(z,2)-2*x*y-2*x*z-2*y*z));
}


//Z_c Masse in Abhängigkeit von eps:
double M(double eps) {
      return (mD+mDs-eps);
}

//Kopplungskonstante gZc in Abhängigkeit von eps und lam aus Datei laden:
double getgZc(int neps, int nlam) {
      
      //Wert
      double gZc = 0;
      
      //welchen Wert laden? (beachte, dass beim Laden der Werte aus der Datei
      //die Werte nacheinander eingelesen werden. Die Werte sind in der Datei also
      //nach folgendem Schema "nummeriert":
      /*
       
        ----------------------
       |  1  2  3  4  5  6  7 |
       |  8  9 10 11 12 13 14 |
       | 15 16 17 18 19 20 21 |
       | ...                  |
       |                      |
        ----------------------
      
      */
      //Um also den richtigen Wert zu laden muss folgende Formel verwendet werden:
      int which = MAXLAMCOUNT * neps + nlam + 1;
      
      ifstream is ("gZc.dat", ifstream::binary);
      if (is) {
            // get length of file:
            is.seekg (0, is.end);
            int length = is.tellg();
            is.seekg (0, is.beg);
            
            char * buffer = new char [length];
            
            // read data as a block:
            is.read (buffer,length);
            is.close();
            
            // ...now buffer contains the entire file...
            
            char* c = buffer;
            for (int i = 0; i < which; i++) //Zeiger auf die Zahl bewegen, die gelesen werden soll
                  gZc = strtod(c, &c); //string to double conversion
            
            delete[] buffer;
      
      } else {
            gZc = 0;
      }
      
      return gZc;
}

//Auch die Kopplungskonstante g ist von eps abhängig
double g(double eps) {
      double geps = (sqrt(6)*gDDsJPi/((pow(M(eps),2)-pow(mJ,2))*(1+pow(mJ,2)/(2*pow(M(eps),2)))));
      return geps;
}

/*********************************************************************/
/*                                                                   */
/*                      Phasenraumgrenzen                            */
/*                                                                   */
/*********************************************************************/
#define s34m (pow(m3+m4,2))
#define s12m (pow(m1+m2,2))
#define xsim 0
#define thtm 0
#define phim (-M_PI)
#define xsip M_PI
#define thtp M_PI
#define phip M_PI

double s34p(double eps) {
      return (pow(M(eps)-m1-m2,2));
}

double s12p(double s34, double eps) {
      return (pow(M(eps)-sqrt(s34),2));
}

/*********************************************************************/
/*                                                                   */
/*                            Formfaktoren                           */
/*                                                                   */
/*********************************************************************/
static int Formfaktoren(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *userdata) {
      
      //Feynmanparameter
	#define u xx[0]
	#define v xx[1]
      #define w xx[2]
      #define t xx[2] //Feynmanparamter für die t-Integration von Diagramm 4
      //einzelne Formfaktoren (also Komponenten eines Vektors f)
	#define ff1 ff[0]
	#define ff23 ff[1]
      #define ff4 ff[2]
      #define ff5 ff[3]
      #define ff6 ff[4]
      
      //Umdefinition der Integrationsvariablen
      // => u,v,w -> 0 bis 1
      // => x,y,z -> 0 bis inf
      #define x (u/(1-u))
	#define y (v/(1-v))
	#define z (w/(1-w))
      
      //Jacobian x,y,z -> u,v,w
      //mit verschiedenen Dimensionalitäten
      #define dm1 (pow(u-1,-2))
      #define dm2 (pow((u-1)*(v-1),-2))
      #define dm3 (pow((u-1)*(v-1)*(w-1),-2))
      
      /*
      //Es ist auch eine andere Umdefinition möglich:
      // => u,v,w -> 0 bis 1
      // => x,y,z -> 0 bis inf
      #define x (tan(M_PI/2*u))
      #define y (tan(M_PI/2*v))
      #define z (tan(M_PI/2*w))
      
      //Jacobian x,y,z -> u,v,w
      //mit verschiedenen Dimensionalitäten
      #define dm1 (M_PI/(2*pow(cos(M_PI/2*u),2)))
      #define dm2 (M_PI/(2*pow(cos(M_PI/2*u),2))*M_PI/(2*pow(cos(M_PI/2*v),2)))
      #define dm3 (M_PI/(2*pow(cos(M_PI/2*u),2))*M_PI/(2*pow(cos(M_PI/2*v),2))*M_PI/(2*pow(cos(M_PI/2*v),2)))
      */
      
	//Modellparameter
      double eps = *((double*)userdata);
	double lam = *((double*)userdata+1);
      double s   = pow(lam,-2);
      //Phasenraumvariable
      double s12 = *((double*)userdata+2);
      double s34 = *((double*)userdata+3);
      //Kopplungskonstante laden:
      double gZc = *((double*)userdata+4);
      
	//Notiz: Aufpassen!
      // #defines immer nur mit Klammen setzen!
      
      #define cc (gZc*M(eps)*g(eps)*q/(16*Power(M_PI,2)))
      
      //Definitionen für F1
      #define r1_2 ((Power(s+2*y,2)*Power(M(eps),2))/4.)
      #define a1 (s+x+y)
      #define M1_2 (Power(MDs,2)*x+Power(MD,2)*y-((s+4*y)*Power(M(eps),2))/4.)
      ff1 = dm2*2*cc*exp(-r1_2/a1-M1_2)/(Power(a1,2));
      
      //Definitionen für F23
      #define r23_2 ((2*y*(-(s*s12)+s*s34+2*s34*y-2*s12*z+2*s34*z)+(s+2*z)*(s+2*(y+z))*Power(M(eps),2))/4.)
      #define a23 (s+x+y+z)
      #define M231_2 (-(s34*y)+Power(mD,2)*(x+y)+Power(mDs,2)*z-((s+4*z)*Power(M(eps),2))/4.)
      #define M232_2 (-(s34*y)+Power(mDs,2)*(x+y)+Power(mD,2)*z-((s+4*z)*Power(M(eps),2))/4.)
      ff23 = dm3*cc*(s+2*z)*(exp(-r23_2/a23-M231_2)+exp(-r23_2/a23-M232_2))/(Power(a23,3));
      
      //Definitionen für F4
      #define r4_2 (((s*(-1+2*t)-3*x+y)*(Power(m3,2)*(s*(-1+2*t)-x-y)+2*s12*(-x+y))+2*(s*(-1+2*t)-x-y)*(x-y)*Power(M(eps),2))/16.)
      #define a4 (s+x+y)
      #define M41_2 ((Power(m3,2)*(-s-3*x+y)+2*(8*Power(mD,2)*x+8*Power(mDs,2)*y-s12*(3*x+y))+2*(x-y)*Power(M(eps),2))/16.)
      #define M42_2 ((Power(m3,2)*(-s-3*x+y)+2*(8*Power(mDs,2)*x+8*Power(mD,2)*y-s12*(3*x+y))+2*(x-y)*Power(M(eps),2))/16.)
      ff4 = -dm2*cc*(s*(x-y))/2*(exp(-r4_2/a4-M41_2)+exp(-r4_2/a4-M42_2))/(Power(a4,3));
      
      //Definitionen für F5
      #define r5_2 ((Power(s+2*y,2)*Power(M(eps),2))/4.)
      #define a5 (s+x+y)
      #define M5_2 (Power(MDs,2)*x+Power(MD,2)*y-((s+4*y)*Power(M(eps),2))/4.)
      ff5 = dm2*4*cc*exp(-r5_2/a5-M5_2)/(Power(a5,2)); //(Power(m1,2)-s34) //Der Propagator steht für diese Rechnung im Matrixelement
      
      //Definitionen für F6
      #define r6_2 ((-((s+2*y)*(s*(s12-2*s34)+2*s12*x-2*s34*(x+y)))+2*(s+2*x)*(s+x+y)*Power(M(eps),2))/4.)
      #define a6 (s+x+y)
      #define M6_2 ((s*s12-2*s*s34+4*Power(MDs,2)*x+4*Power(MD,2)*y-4*s34*y-2*(s+2*x)*Power(M(eps),2))/4.)
      ff6 = dm2*4*cc*exp(-r6_2/a6-M6_2)/(Power(a6,2)*(Power(M(eps),2)-s12));
      
	return 0;
}

/*********************************************************************/
/*                                                                   */
/*                      Feynmanintegration                           */
/*                                                                   */
/*********************************************************************/
double *F(double eps, double lam, double s12, double s34, double xsi, double tht, double phi, double gZc) {
      //Rückgabewerte für die Feynmanintegration
	int nregionsf, nevalf, failf;
	double FF[NCOMPF], errorf[NCOMPF], probf[NCOMPF];
	
	double dataf[] = {eps, lam, s12, s34, gZc};
	
#if INTMETH == CUHRE
      Cuhre(NDIMF, NCOMPF, Formfaktoren, dataf,
            EPSRELF, EPSABSF, VERBOSEF | LASTF,
            MINEVALF, MAXEVALF, KEYF,
            STATEFILEF,
            &nregionsf, &nevalf, &failf, FF, errorf, probf);
#endif
#if INTMETH == SUAVE
      Suave(NDIMF, NCOMPF, Formfaktoren, dataf,
            EPSRELF, EPSABSF, VERBOSEF | LASTF, SEEDF,
            MINEVALF, MAXEVALF, NNEWF, FLATNESSF,
            STATEFILEF,
            &nregionsf, &nevalf, &failf, FF, errorf, probf);
#endif
#if INTMETH == VEGAS
      Vegas(NDIMF, NCOMPF, Formfaktoren, dataf,
            EPSRELF, EPSABSF, VERBOSEF | LASTF, SEEDF,
            MINEVALF, MAXEVALF, NSTARTF, NINCREASEF, NBATCHF,
            GRIDNOF, STATEFILEF,
            &nevalf, &failf, FF, errorf, probf);
#endif

#if INTMETH == DIVONNE
      Divonne(NDIMF, NCOMPF, Formfaktoren, dataf,
              EPSRELF, EPSABSF, VERBOSEF | LASTF, SEEDF,
              MINEVALF, MAXEVALF, KEY1F, KEY2F, KEY3F, MAXPASSF,
              BORDERF, MAXCHISQF, MINDEVIATIONF,
              NGIVENF, LDXGIVENF, NULL, NEXTRAF, NULL,
              STATEFILEF,
              &nregionsf, &nevalf, &failf, FF, errorf, probf);
#endif
      
      //Speichern der Ergebnisse in einem Vektorn total[]
      for (int i = 0; i < NCOMPF; i++) {
            total[i] = FF[i];
            total[i + NCOMPF] = errorf[i];
            total[i + 2 * NCOMPF] = probf[i];
      }
      total[3 * NCOMPF]     = (double)nregionsf;
      total[1 + 3 * NCOMPF] = (double)nevalf;
      total[2 + 3 * NCOMPF] = (double)failf;
      
      return total;
}

#include "Minv.c"

/*********************************************************************/
/*                                                                   */
/*               Betragsquadriertes Matrixelement                    */
/*                                                                   */
/*********************************************************************/
double Minv(double eps, double lam, double s12, double s34, double xsi, double tht, double phi, double gZc) {
      double F1, F23, F4, F5, F6, minv;
      double errF1, errF23, errF4, errF5, errF6;
      int nregionsf, nevalf, failf;
      
      double erg[3 * NCOMPF + 3];
      
      //Feynmanintegration der Formfaktoren durchführen
      double *ptr = F(eps, lam, s12, s34, xsi, tht, phi, gZc);
      
      //Ergebnisse in einem Vektor erg[] speichern
      for (int i = 0; i < NCOMPF; i++) {
            erg[i] = *(ptr + i);
            erg[i + NCOMPF] = *(ptr + i + NCOMPF);
            erg[i + 2 * NCOMPF] = *(ptr + i + 2 * NCOMPF);
      }
      nregionsf = (int)(*(ptr + 3 * NCOMPF));
      nevalf = (int)(*(ptr + 3 * NCOMPF + 1));
      failf = (int)(*(ptr + 3 * NCOMPF + 2));
      
      //Formfaktoren
      F1  = erg[0];
      F23  = erg[1];
      F4  = erg[2];
      F5  = erg[3];
      F6  = erg[4];
      
      //Fehler der Formfaktoren
      errF1  = absolut(erg[5]);
      errF23  = absolut(erg[6]);
      errF4  = absolut(erg[7]);
      errF5  = absolut(erg[8]);
      errF6  = absolut(erg[9]);
      
      //Matrixelement laden
      minv = getMinv(eps, lam, s12, s34, xsi, tht, phi, F1, F23, F4, F5, F6);
      
      
      //Zwischenwerte auf Bildschirm ausgeben und in Dateien schreiben
#ifndef DEBUG
 #ifndef OPTIMIZED
      getStatus(s12, s34, xsi, tht, phi, gZc, g(eps), eps, lam, F1, errF1, F23, errF23, F4, errF4,
                F5, errF5, F6, errF6, nregionsf, nevalf, failf, minv);
 #endif
#endif
      
#ifdef OPTIMIZED
      if((cyclePS % 100) == 1)
            printf("%.0f\n",getPercent());
#endif
      
      move(0,0);
      refresh();
      
      if(!isnan(minv)) {
          return minv;
      }
    
      isnanFail++; //Falls das Matrixelement ein NaN ist.
      return 0;
}

/*********************************************************************/
/*                                                                   */
/*                  Definition des Phasenraums                       */
/*                                                                   */
/*********************************************************************/
static int Phasenraum(const int *ndim, const double xx[],
                        const int *ncomp, double ff[], void *userdata) {
      
      //Phasenraumvariable
      #define qxsi xx[0]
      #define qtht xx[1]
      #define qphi xx[2]
      #define q12 xx[3]
      #define PS ff[0]
      
	//Modellparameter
      double eps = *((double*)userdata);
	double lam = *((double*)userdata+1);
	double gZc = *((double*)userdata+2);
      double s34 = *((double*)userdata+3);
      
      //Konstante für Phasenraumintegration
      #define c (1.0/(3.0*pow(2.0,15.0)*pow(M(eps),3.0)*pow(M_PI,6.0)))
      
      //Umdefinition der Variaben [0, 1] -> [s34m, s34p]:
      double s12 = s12m - q12*(s12m-s12p(s34,eps));
      double ds12 = (s12p(s34,eps)-s12m);
      
      //Umdefinition der Winkel
      double xsi = xsim - qxsi*(xsim - xsip);
      double dxsi = (xsip - xsim);
      double tht = thtm - qtht*(thtm - thtp);
      double dtht = (thtp - thtm);
      double phi = phim - qphi*(phim - phip);
      double dphi = (phip - phim);
      
      //Jacobian der Phasenraumintegration und sonstige Vorfaktoren
      double jac = 1/(s34*s12)*ksqrt(pow(M(eps),2),s12,s34)*ksqrt(s12,pow(m1,2),pow(m2,2))*ksqrt(s34,pow(m3,2),pow(m4,2))*ds12*dxsi*dtht*dphi*sin(tht)*sin(xsi);
      
      //Phasenraumintegrant mit angepassten Grenzen 0 bis 1
      PS = c*Minv(eps, lam, s12, s34, xsi, tht, phi, gZc)*jac;
      
      refresh();
      
	return 0;
}

/*********************************************************************/
/*                                                                   */
/*        Winkelintegration und Integration von s12                  */
/*                                                                   */
/*********************************************************************/
void integrate(double s34, double eps, double lam, double gZc) {
      //Wie oft wurde bereits gerechnet?
      cyclePS++;
      
      //Rückgabewerte für die Phasenraumintegration
	int nregionsps, nevalps, failps;
	double integralps[NCOMPPS], errorps[NCOMPPS], probps[NCOMPPS];
      double dataps[] = {eps, lam, gZc, s34};
      
      //Phasenraumintegration
#if INTMETH == CUHRE
      Cuhre(NDIMPS, NCOMPPS, Phasenraum, dataps,
            EPSRELPS, EPSABSPS, VERBOSEPS | LASTPS,
            MINEVALPS, MAXEVALPS, KEYPS,
            STATEFILEPS,
            &nregionsps, &nevalps, &failps, integralps, errorps, probps);
#endif
#if INTMETH == SUAVE
      Suave(NDIMPS, NCOMPPS, Phasenraum, dataps,
            EPSRELPS, EPSABSPS, VERBOSEPS | LASTPS, SEEDPS,
            MINEVALPS, MAXEVALPS, NNEWPS, FLATNESSPS,
            STATEFILEPS,
            &nregionsps, &nevalps, &failps, integralps, errorps, probps);
#endif
#if INTMETH == DIVONNE
      Divonne(NDIMPS, NCOMPPS, Phasenraum, dataps,
              EPSRELPS, EPSABSPS, VERBOSEPS | LASTPS, SEEDPS,
              MINEVALPS, MAXEVALPS, KEY1PS, KEY2PS, KEY3PS, MAXPASSPS,
              BORDERPS, MAXCHISQPS, MINDEVIATIONPS,
              NGIVENPS, LDXGIVENPS, NULL, NEXTRAPS, NULL,
              STATEFILEPS,
              &nregionsps, &nevalps, &failps, integralps, errorps, probps);
#endif
#if INTMETH == VEGAS
      Vegas(NDIMPS, NCOMPPS, Phasenraum, dataps,
            EPSRELPS, EPSABSPS, VERBOSEPS | LASTPS, SEEDPS,
            MINEVALPS, MAXEVALPS, NSTARTPS, NINCREASEPS, NBATCHPS,
            GRIDNOPS, STATEFILEPS,
            &nevalps, &failps, integralps, errorps, probps);
#endif
      
#ifndef OPTIMIZED
      string filename = "dGamma2_eps";
      filename.append(to_string((int)(eps*1000000))); //Epsilon in keV an den Dateiname anhängen
      filename.append("keV_lam");
      filename.append(to_string((int)(lam*1000))); //Lambda in MeV
      filename.append("MeV.dat");
      dGamma.open(filename.c_str(), ios::out | ios::app);
      dGamma << s34 << "\t" << integralps[0] << "\n";
      dGamma.close();
#endif
      
}


/*********************************************************************/
/*                                                                   */
/*                            Main()                                 */
/*                                                                   */
/*********************************************************************/
int main() {
      time_t start, end;
      start = time (NULL);

#ifndef OPTIMIZED
      initscr();
      //wresize(stdscr,SCREENHEIGHT,SCREENWIDTH); //Resize the terminal
#endif

	for (int neps = STRTEPS; neps < EPSCOUNT; neps++) {
		for (int nlam = STRTLAM; nlam < LAMCOUNT; nlam++) {
			
                  //Aktuelle Parameter berechnen
                  double eps = epsilon + epsilonStep * neps;
                  double lam = lambda + lambdaStep * nlam;
                  //Kopplungskonstante laden und für aktuelle Rechnung speichern
                  double gZc = getgZc(neps, nlam);
                  //aktuelle Bindungsenergie und Größe in Vektor schreiben, wird nachher der Integration übergeben
                  double dataps[] = {eps, lam, gZc};
                  
#ifndef OPTIMIZED
                  //Dateiausgabe vorbereiten
                  string filename = "dGamma2_eps";
                  filename.append(to_string((int)(eps*1000000))); //Epsilon in keV an den Dateiname anhängen
                  filename.append("keV_lam");
                  filename.append(to_string((int)(lam*1000))); //Lambda in MeV
                  filename.append("MeV.dat");
                  dGamma.open(filename.c_str(), ios::out); //alte Dateien löschen!
                  //Header
                  dGamma << "eps: " << eps*1000 << "MeV" << endl;
                  dGamma << "lam: " << lam*1000 << "MeV" << endl;
                  dGamma << "mGam: " << mGam*1000 << "MeV" << endl;
                  dGamma << "MAXEVALF: " << MAXEVALF << endl;
                  dGamma << "MAXEVALPS: " << MAXEVALPS << endl;
                  dGamma << endl;
                  //Tabelle vorbereiten
                  dGamma << "s34\t\tdGamma" << endl;
                  dGamma.close();
#endif

                  //Integrationsschleife für die Punkte (s12,s34)
                  for (int ns34 = 0; ns34 <= MAXNS; ns34++) {
                              double q34 = (double)ns34/(double)MAXNS; //läuft von [~0,1]
                              double s34 = s34m + (s34p(eps) - s34m)*q34; //Invariante Massen
                              integrate(s34, eps, lam, gZc);
                  }
		}
	}
#ifndef OPTIMIZED
      move(21,0);
#endif
      
#ifdef OPTIMIZED
      initscr();
      move(0,0);
#endif

      printw("All done. Good bye.\n");
      refresh();
      getch();
      
      endwin();
      
      return 0;
}

