#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define TH       0.000000001
#define NBR_ABS  8
#define NBR_FLUO 8

typedef struct {
   char  output[8]; 
   double Fe055, Am241, Co057, Cs137;
   double CdTe_absd, Si_absd;
   char  abs_mat[NBR_ABS][8];
   double abs_thk[NBR_ABS];
   char  fluo_mat[NBR_FLUO][8];
   double fluo_int[NBR_FLUO];
   char  det_mat[8];
   double det_d, det_FM, det_res0;
   double t_out, f_in;
   double bin_w;
} config_t;

typedef struct {
   double erg;
   double in;
} spec_t;

double binning;  //in eV; set by configuration file

int   e2i(double erg);
double i2e(int i);
double myabs(double a);
void  InitConfiguration(config_t* c); 
void  ReadConfiguration(config_t* c);
void  AddFe055(spec_t* LineSpec,  int length, double inten);
void  AddAm241(spec_t* LineSpec, int length, double inten);
void  AddCo057(spec_t* LineSpec,  int length, double inten);
void  AddCs137(spec_t* LineSpec, int length, double inten);
void  ApplyAbsorber(spec_t* LineSpec, int length, config_t* c);
void  AddFluorescence(spec_t* LineSpec, int length, config_t c);
void  CalcCompt(spec_t* LineSpec, spec_t* ComSpec, int length, config_t c);
double DsigIS_De(double Erg, double k0);
double DsigCS_De(double k, double k0);
double S(char* element, double theta, double k0);
double theta(double k, double k0);
double theta_bar(double k, double k0);
void  ApplyDetectorEfficiency(spec_t* LineSpec, int length, config_t c);
double Transmission(int i_E, double d, char* material);
void  ApplyDetRes(spec_t* LineSpec, spec_t* DetSpec, int length, config_t c);
void  AddSpectra(spec_t* DetSpec, spec_t* ComSpec, int length);
double Efficiency(int energy, double d, char* mat);
double Gauss(double x0, double x, double height, config_t c);
double myGauss(double x0, double x, double height, config_t c);
void  Normalize(spec_t* DetSpec, int length);
void  ApplyEscape(spec_t* LineSpec, int length, config_t c);
void  AddComptonization(spec_t* ComSpec, spec_t* DetSpec, int length);
void  ApplyPileUp(spec_t* DetSpec, spec_t* PUSpec, int length, config_t* c);
void  WriteOutput(spec_t* DetSpec, int length);

int main(int argc, char* argv[])
{
   printf("###########################################################\n");
   printf("###########################################################\n");
   printf("#                           SynSpec                       #\n");
   printf("###########################################################\n");
   printf("###############################################27.12.2016##\n");
   if (argc==1) {
      printf("SynSpec is a tool to create a syntetic X-ray spectrum.\n");
      printf("All configuration are set by the SynSpec_config.txt file.\n");
      printf("To plot sepctrum use: gnuplot SynSpec.plt Synspec_out.txt\n");
   }
   
   int i;
   config_t c;
   InitConfiguration(&c);
   ReadConfiguration(&c);//SynSpec_config.txt
   
   spec_t *LineSpec, *DetSpec, *ComSpec, *PUSpec;
   int length = 2.0E4;//max. energy = length * binning [eV]
   LineSpec = (spec_t*)malloc(length * sizeof(spec_t));
   DetSpec  = (spec_t*)malloc(length * sizeof(spec_t));
   ComSpec  = (spec_t*)malloc(length * sizeof(spec_t));
   PUSpec   = (spec_t*)malloc(length * sizeof(spec_t));
   if (LineSpec == NULL) {printf("ERROR: cannot alloc LineSpec!\n");}
   if (DetSpec  == NULL) {printf("ERROR: cannot alloc DetSpec!\n");}
   if (ComSpec  == NULL) {printf("ERROR: cannot alloc ComSpec!\n");}
   if (PUSpec   == NULL) {printf("ERROR: cannot alloc PUSpec!\n");}
 
   double energy;
   for (i=0; i<=length-1; ++i) {
      LineSpec[i].in  = 0;
      DetSpec[i].in  = 0;
      ComSpec[i].in  = 0;
      energy = i2e(i);
      LineSpec[i].erg = energy;
      DetSpec[i].erg  = energy;
      ComSpec[i].erg  = energy;
      PUSpec[i].erg   = energy;
   }

   if (c.Fe055>0) {AddFe055(LineSpec, length, c.Fe055);}
   if (c.Am241>0) {AddAm241(LineSpec, length, c.Am241);}
   if (c.Co057>0) {AddCo057(LineSpec, length, c.Co057);}
   if (c.Cs137>0) {AddCs137(LineSpec, length, c.Cs137);}

   LineSpec[e2i(  60000)].in += 0.09; //testing...
   LineSpec[e2i( 120000)].in += 1.0; //testing...

   AddFluorescence(LineSpec, length, c);

   ApplyAbsorber(LineSpec, length, &c);

   ApplyEfficiency_PE(LineSpec, length, c);


   //Comptonization is not jet added to LineSpec in order to keep
   //the spectrum non-continuous
   CalcCompt(LineSpec, ComSpec, length, c);
   printf("Comptonization calculated.\n");
   
   //ApplyEscape(LineSpec, length, c);

   //ApplyDetRes(LineSpec, DetSpec, length, c);
   AddSpectra(LineSpec, ComSpec, length);
   //printf("detector resolution applied.\n");

   //ApplyPileUp(LineSpec, PUSpec, length, &c);
   //AddSpectra(LineSpec, PUSpec, length);
   //Normalize(LineSpec, length);

   WriteOutput(LineSpec, length);

   free(LineSpec);
   free(DetSpec);
   free(ComSpec);
   free(PUSpec);
   return 0;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void  InitConfiguration(config_t* c)
{
   int i;

   strcpy(c->output, "txt"); 
   c->Fe055 = 0.0;
   c->Am241 = 0.0; 
   c->Co057 = 0.0; 
   c->Cs137 = 0.0; 
   for (i=0; i<NBR_ABS; ++i) {
      strcpy(c->abs_mat[i], "");
      c->abs_thk[i] = 0.0;
   }
   for (i=0; i<NBR_FLUO; ++i) {
      strcpy(c->fluo_mat[i], "");
      c->fluo_int[i] = 0.0;
   }
   strcpy(c->det_mat, "");
   c->det_d    = 0.0;
   c->det_FM   = 1.0;
   c->det_res0 = 100;
   c->t_out    = 100;
   c->f_in     =   1;
   c->bin_w    = 100;
}
//------------------------------------------------------------------------------
void ReadConfiguration(config_t* c)
{
   int i;
   FILE* fpin_config;
   fpin_config = fopen("SynSpec_config.txt", "r");
   if (fpin_config == NULL) {
      printf("ERROR: cannot open SynSpec_config.txt\n");
   }
   
   char  line[256];
   while (fgets(line, sizeof(line), fpin_config)) {
      if (line[0]!='#') {
         sscanf(line, " output    = %s", (*c).output);
         sscanf(line, " Fe-55     = %lf", &((*c).Fe055));
         sscanf(line, " Am-241    = %lf", &((*c).Am241));
         sscanf(line, " Co-57     = %lf", &((*c).Co057));
         sscanf(line, " Cs-137    = %lf", &((*c).Cs137));
         sscanf(line, " abs_mat0  = %s", (*c).abs_mat[0]);
         sscanf(line, " abs_thk0  = %lf", &((*c).abs_thk[0]));
         sscanf(line, " abs_mat1  = %s", (*c).abs_mat[1]);
         sscanf(line, " abs_thk1  = %lf", &((*c).abs_thk[1]));
         sscanf(line, " abs_mat2  = %s", (*c).abs_mat[2]);
         sscanf(line, " abs_thk2  = %lf", &((*c).abs_thk[2]));
         sscanf(line, " abs_mat3  = %s", (*c).abs_mat[3]);
         sscanf(line, " abs_thk3  = %lf", &((*c).abs_thk[3]));
         sscanf(line, " abs_mat4  = %s", (*c).abs_mat[4]);
         sscanf(line, " abs_thk4  = %lf", &((*c).abs_thk[4]));
         sscanf(line, " abs_mat5  = %s", (*c).abs_mat[5]);
         sscanf(line, " abs_thk5  = %lf", &((*c).abs_thk[5]));
         sscanf(line, " abs_mat6  = %s", (*c).abs_mat[6]);
         sscanf(line, " abs_thk6  = %lf", &((*c).abs_thk[6]));
         sscanf(line, " abs_mat7  = %s", (*c).abs_mat[7]);
         sscanf(line, " abs_thk7  = %lf", &((*c).abs_thk[7]));
         sscanf(line, " fluo_mat0 = %s", (*c).fluo_mat[0]);
         sscanf(line, " fluo_int0 = %lf", &((*c).fluo_int[0]));
         sscanf(line, " fluo_mat1 = %s", (*c).fluo_mat[1]);
         sscanf(line, " fluo_int1 = %lf", &((*c).fluo_int[1]));
         sscanf(line, " fluo_mat2 = %s", (*c).fluo_mat[2]);
         sscanf(line, " fluo_int2 = %lf", &((*c).fluo_int[2]));
         sscanf(line, " fluo_mat3 = %s", (*c).fluo_mat[3]);
         sscanf(line, " fluo_int3 = %lf", &((*c).fluo_int[3]));
         sscanf(line, " fluo_mat4 = %s", (*c).fluo_mat[4]);
         sscanf(line, " fluo_int4 = %lf", &((*c).fluo_int[4]));
         sscanf(line, " fluo_mat5 = %s", (*c).fluo_mat[5]);
         sscanf(line, " fluo_int5 = %lf", &((*c).fluo_int[5]));
         sscanf(line, " fluo_mat6 = %s", (*c).fluo_mat[6]);
         sscanf(line, " fluo_int6 = %lf", &((*c).fluo_int[6]));
         sscanf(line, " fluo_mat7 = %s", (*c).fluo_mat[7]);
         sscanf(line, " fluo_int7 = %lf", &((*c).fluo_int[7]));
         sscanf(line, " CdTe_absd = %lf", &((*c).CdTe_absd));
         sscanf(line, " Si_absd   = %lf", &((*c).Si_absd));
         sscanf(line, " det_mat   = %s", (*c).det_mat);
         sscanf(line, " det_d     = %lf", &((*c).det_d));
         sscanf(line, " det_FM    = %lf", &((*c).det_FM));
         sscanf(line, " det_res0  = %lf", &((*c).det_res0));
         sscanf(line, " t_out     = %lf", &((*c).t_out));
         sscanf(line, " f_in      = %lf", &((*c).f_in));
         sscanf(line, " bin_width = %lf", &((*c).bin_w));
      }
   }
   
   printf("## reading config file...\n");
   printf("## OUTPUT: %s\n", c->output);
   printf("## SOURCE:\n");
   if (c->Fe055>0)     {printf("   Fe-55   I=%5.3f\n", c->Fe055);} 
   if (c->Am241>0)     {printf("   Am-241  I=%5.3f\n", c->Am241);}
   if (c->Co057>0)     {printf("   Co-57   I=%5.3f\n", c->Co057);}
   if (c->Cs137>0)     {printf("   Cs-137  I=%5.3f\n", c->Cs137);}
   printf("## ABSORBER:\n");
   for (i=0; i<NBR_ABS; ++i) {
      if (c->abs_thk[i] > TH) { 
         printf("i=%d   %s    d=%5.3f mm\n", i, c->abs_mat[i], c->abs_thk[i]);
      }
   }
   printf("## DETECTOR:\n");
   printf("   %s     d=%5.3f mm, res=%5.3f * F + %5.3f eV FWHM\n", 
          c->det_mat, c->det_d, c->det_FM, c->det_res0);
   printf("   t_out = %5.3fms,  f_in = %5.3f Hz/ch\n", c->t_out, c->f_in);
   printf("## CALCULATION:\n");
   binning = c->bin_w;
   printf("   bin_w: %5.2f eV\n", binning);
   fclose(fpin_config);
}
//------------------------------------------------------------------------------
void ApplyAbsorber(spec_t* LineSpec, int length, config_t* c)
{
   int i, i_E; //index for LineSpec
   for (i=0; i<NBR_ABS; ++i) {
      if (c->abs_thk[i] > 0.0) {
         for (i_E=0; i_E<= length-1; ++i_E) {
            if(LineSpec[i_E].in > 0.0) {
               printf("i_E=%d  E = %f  -> %f\n", 
                       i_E, LineSpec[i_E].erg, LineSpec[i_E].in);
               LineSpec[i_E].in *= Transmission(i_E, c->abs_thk[i], 
                                               c->abs_mat[i]);
            }
         }
      }
   }
}
//------------------------------------------------------------------------------
void  AddFluorescence(spec_t* LineSpec, int length, config_t c)
{
   //to do: implement more than Ka1, adjust intensity according to 
   //fluorescence yield and according to the amount of radiation that
   //can produce the fluorescence.
   int i;
   for (i=0; i<NBR_FLUO; ++i) {
      if (!strcmp(c.fluo_mat[i], "Al")) {
         if (  1487<length-1) {LineSpec[  1487].in += c.fluo_int[i];}
      }
      if (!strcmp(c.fluo_mat[i], "Si")) {
         if (  1740<length-1) {LineSpec[  1740].in += c.fluo_int[i];}
      }
      if (!strcmp(c.fluo_mat[i], "Ti")) {
         if (  4511<length-1) {LineSpec[  4511].in += c.fluo_int[i];}
      }
      if (!strcmp(c.fluo_mat[i], "Fe")) {
         if (  6406<length-1) {LineSpec[  6406].in += c.fluo_int[i];}
      }
      if (!strcmp(c.fluo_mat[i], "Cu")) {
         if (  8048<length-1) {LineSpec[  8048].in += c.fluo_int[i];}
      }
      if (!strcmp(c.fluo_mat[i], "Ge")) {
         if (  9886<length-1) {LineSpec[  9886].in += c.fluo_int[i];}
      }
      if (!strcmp(c.fluo_mat[i], "Ag")) {
         if ( 22163<length-1) {LineSpec[ 22163].in += c.fluo_int[i];}
      }
      if (!strcmp(c.fluo_mat[i], "Cd")) {
         if ( 23174<length-1) {LineSpec[ 23174].in += c.fluo_int[i];}
      }
      if (!strcmp(c.fluo_mat[i], "Sn")) {
         if ( 25271<length-1) {LineSpec[ 25271].in += c.fluo_int[i];}
      }
      if (!strcmp(c.fluo_mat[i], "Te")) {
         if ( 27472<length-1) {LineSpec[ 27472].in += c.fluo_int[i];}
      }
      if (!strcmp(c.fluo_mat[i], "W")) {
         if ( 59318<length-1) {LineSpec[ 59318].in += c.fluo_int[i];}
      }
      if (!strcmp(c.fluo_mat[i], "Pt")) {
         if ( 66832<length-1) {LineSpec[ 66832].in += c.fluo_int[i];}
      }
      if (!strcmp(c.fluo_mat[i], "Au")) {
         if ( 68804<length-1) {LineSpec[ 68804].in += c.fluo_int[i];}
      }
      if (!strcmp(c.fluo_mat[i], "Pb")) {
         if ( 74976<length-1) {LineSpec[ 74969].in += c.fluo_int[i];}
      }
   }
}
//------------------------------------------------------------------------------
void CalcCompt(spec_t* LineSpec, spec_t* ComSpec, int length, config_t c)
{
   int   i_E, ii_E;
   double ComptEdge, LineIn;
   for (i_E=0; i_E<=length-1; ++i_E) {
      if (LineSpec[i_E].in > TH) {
         LineIn    = LineSpec[i_E].in;
         ComptEdge = i2e(i_E)*(1.0 - 1.0/(1.0 + (2.0*i2e(i_E)/511000.0)));
         for (ii_E=0; ii_E<=e2i(ComptEdge); ++ii_E) {
            ComSpec[ii_E].in += LineIn * DsigIS_De(i2e(ii_E), i2e(i_E));
         }
      }
   }
}
//------------------------------------------------------------------------------
double DsigIS_De(double Erg, double k0)
{
   double k = k0 - Erg; //energy of sc. photon
   return DsigCS_De(k, k0) * sin(theta(k, k0)) * theta_bar(k, k0) *
          S("CdTe", theta(k, k0), k0);
}
//------------------------------------------------------------------------------
double DsigCS_De(double k, double k0)
{
   const double r_0_2 = 0.07740787; //sq of classical el. radius
   const double kdk0  = k/k0;
   return  M_PI * r_0_2 * (kdk0)*(kdk0) * 
           (kdk0 + k0/k - (sin(theta(k, k0)))*(sin(theta(k, k0))));
}
//------------------------------------------------------------------------------
double S(char* element, double theta, double k0)
{
   const double k0_keV = k0 / 1000.0;
   const double t_2    = theta / 2.0;
   if (!strcmp(element, "Si"))
      return 14.0 * (1-exp(-0.058942 * pow(k0_keV*sin(t_2), 1.1289))) / 
                    (1-0.638*exp(-0.09432*k0_keV*sin(t_2)));
   if (!strcmp(element, "Ge"))
      return 32.0 * (1-exp(-0.039373 * pow(k0_keV * sin(t_2), 1.1002))) / 
                    (1-0.570*exp(-0.05860*k0_keV*sin(t_2)));
   if (!strcmp(element, "Sn") || !strcmp(element, "CdTe"))
      return 50.0 * (1-exp(-0.025066 * pow(k0_keV * sin(t_2), 1.0714))) / 
                    (1-0.702*exp(-0.03408*k0_keV*sin(t_2)));
   return -1; //indicates error
}
//------------------------------------------------------------------------------
double theta(double k, double k0)
{
   return acos(1.0 - (511000.0*k0 / k-511000.0) / k0);
}
//------------------------------------------------------------------------------
double theta_bar(double k, double k0)
{
   return 1.0 / (sqrt(1.0 - (1.0 + 511000.0 * (k-k0)/(k*k0)) *
                            (1.0 + 511000.0 * (k-k0)/(k*k0)))) * 511000.0/(k*k);
}
//------------------------------------------------------------------------------
void AddComptonization(spec_t* ComSpec, spec_t* DetSpec, int length)
{
   int i_D;
   for (i_D=0; i_D<=length-1; ++i_D) {
      DetSpec[i_D].in += ComSpec[i_D].in;
   }
}
//------------------------------------------------------------------------------
double Transmission(int i_E, double d, char* mat)
{
   char file_name[32];
   sprintf(file_name, "SynSpec_material_%s.txt", mat);
   FILE* fp_mat;
   fp_mat = fopen(file_name, "r");
   if (fp_mat==NULL) {
      printf("ERROR: cannot open material definition file %s\n", file_name);
   }

   char  line[256];
   double erg0=i2e(i_E)/1E6; //in MeV
   double rho; //density in g/cm^3
   double att; //in cm^2/g
   double du;  //dump unused informations
   double erg=0, last_erg=0;
   double last_att=0;
   while (fgets(line, sizeof(line), fp_mat)) {
      if (line[0]!='#') {
         sscanf(line, " density    = %lf", &rho);
         sscanf(line, " %lf %lf %lf %lf %lf %lf %lf %lf",
                &erg, &du, &du, &du, &du, &du, &att, &du);
      }
      if ((last_erg < erg0) && (erg0 < erg)) {//interpolate attenuation
         att = last_att + (att-last_att)/(erg-last_erg)*(erg0-last_erg);
         break;
      }
      last_att = att;
      last_erg = erg;
   }
 
   double dcm = d/10.0; //thickness in cm
   fclose(fp_mat);
   printf("i_E = %d, E0=%f, E=%f, att=%f, trans=%f\n", 
           i_E,      erg0,   erg, att,    exp(-att*rho*dcm));
   return exp(-att*rho*dcm);
}
//------------------------------------------------------------------------------
void ApplyEscape(spec_t* LineSpec, int length, config_t c)
{
   int   i_E, j;
   int   E_esc[8];      //escaping energies arranged in order of size
   double P_esc[8];      //escape probabilities
   E_esc[0] = length;   //do not apply escape if no material matches
   if (!strcmp(c.det_mat, "Si")) {
      E_esc[0] = 1740;  //Si-Ka
      P_esc[0] = 0.01;  //
      E_esc[1] = -1;    //indicates end of escape list
      printf("apply Si escape...");
   }
   if (!strcmp(c.det_mat, "Ge")) {
      E_esc[0] = 9886;  //Ge-Ka
      P_esc[0] = 0.01;  //
      E_esc[1] = -1;    //indicates end of escape list
      printf("apply Ge escape...");
   }
   if (!strcmp(c.det_mat, "CdTe")) { 
      E_esc[0] = 23174; //Cd-Ka
      P_esc[0] = 0.01;  //
      E_esc[1] = 27472; //Te-Ka
      P_esc[1] = 0.01;  //
      E_esc[2] = -1;    //indicates end of escape list
      printf("apply CdTe escape...");
   }
   for (i_E=e2i(E_esc[0]); i_E<=length-1; ++i_E) {
      if(LineSpec[i_E].in > TH) {
         for (j=0; j<=7; ++j) {
            if (E_esc[j]==-1) {break;}
            if (i2e(i_E)>E_esc[j]) {
               LineSpec[e2i(i2e(i_E)-E_esc[j])].in += P_esc[j]*LineSpec[i_E].in;
            }
         }
      }
   }
   printf("done %s\n", c.det_mat);
}
//------------------------------------------------------------------------------
void ApplyDetectorEfficiency(spec_t* LineSpec, int length, config_t c)
{
   int i_E;
   for (i_E=0; i_E<=length-1; ++i_E) {
      if(LineSpec[i_E].in > TH) {
         LineSpec[i_E].in *= Efficiency(i_E, c.det_d, c.det_mat); 
      }
   }
}
//------------------------------------------------------------------------------
void ApplyDetRes(spec_t* LineSpec, spec_t* DetSpec, int length, config_t c)
{
   int i_E; //index for LineSpec
   int i_D; //index for DetSpec
   for (i_E=0; i_E<=length-1; ++i_E) {
      //DetSpec[i_E].in = LineSpec[i_E].in;
      if(LineSpec[i_E].in > TH) {
         //printf("apply Gauss for E=%f\n", i2e(i_E));
         for (i_D=0; i_D<= length-1; ++i_D) {
            DetSpec[i_D].in += Gauss(i2e(i_E), i2e(i_D), LineSpec[i_E].in, c);
         }
      }
   } 
}
//------------------------------------------------------------------------------
void AddSpectra(spec_t* DetSpec, spec_t* ComSpec, int length)
{
   int i_E; //index for LineSpec
   for (i_E=0; i_E<=length-1; ++i_E) DetSpec[i_E].in += ComSpec[i_E].in;
}
//------------------------------------------------------------------------------
double Efficiency(int energy, double d, char* mat)
{
   return 1-Transmission(energy, d, mat);
}
//------------------------------------------------------------------------------
double Gauss(double x0, double x, double height, config_t c)
{
   double var, omega=1, F=1;
   if (!strcmp(c.det_mat, "Si")) {
      omega = 3.63;
      F = 0.10;
   } else {
      if (!strcmp(c.det_mat, "CdTe")) {
         omega = 4.43;
         F = 0.15;
      } else {
         printf("ERROR: unknown detector material: %s\n", c.det_mat);
      }
   }
   var = sqrt((c.det_FM*c.det_FM * F * x0 * omega) * 
              (c.det_FM*c.det_FM * F * x0 * omega) + 
              pow(c.det_res0/2.3548, 4));
   return height/sqrt(2*M_PI*var) * exp(-(x-x0)*(x-x0) / (2*var));
}
//------------------------------------------------------------------------------
double myGauss(double x0, double x, double height, config_t c)
{
   double var, omega=1, F=1;
   if (!strcmp(c.det_mat, "Si")) {
      omega = 3.63;
      F = 0.10;
   } else {
      if (!strcmp(c.det_mat, "CdTe")) {
         omega = 4.43;
         F = 0.15;
      } else {
         printf("ERROR: unknown detector material: %s\n", c.det_mat);
      }
   }
   var = sqrt((c.det_FM * F * x0 * omega) * (c.det_FM * F * x0 * omega) + 
              pow(c.det_res0/(2.3548*omega), 4));
   return height * exp(-(x-x0)*(x-x0) / (2*var));
}
//------------------------------------------------------------------------------
void Normalize(spec_t* DetSpec, int length)
{
   int i;
   double max_val=0;
   for (i=0; i<=length-1; ++i) {
      if (DetSpec[i].in > max_val) {
         max_val = DetSpec[i].in;
      }
   } 
   for (i=0; i<=length-1; ++i) {
      DetSpec[i].in = DetSpec[i].in/max_val;
   }
}
//------------------------------------------------------------------------------
void AddCs137(spec_t* LineSpec, int length, double inten)
{
   if (length-1 < e2i(661657)) {
      printf("ERROR: length of LineSpec is not long enough!\n");
   } else {
      LineSpec[e2i(  3954)].in += inten * 0.0090;
      LineSpec[e2i( 31817)].in += inten * 0.0195;
      LineSpec[e2i( 32194)].in += inten * 0.0359;
      LineSpec[e2i( 36305)].in += inten * 0.01055;
      LineSpec[e2i( 36379)].in += inten * 0.01055;
      LineSpec[e2i( 36654)].in += inten * 0.01055;
      LineSpec[e2i( 37258)].in += inten * 0.00266;
      LineSpec[e2i( 37312)].in += inten * 0.00266;
      LineSpec[e2i( 37425)].in += inten * 0.00266;
      LineSpec[e2i(283500)].in += inten * 0.0000058;
      LineSpec[e2i(661657)].in += inten * 0.8499;
   }
}
//------------------------------------------------------------------------------
void AddAm241(spec_t* LineSpec, int length, double inten)
{
   if (length-1 < e2i(59541)) {
      printf("ERROR: length of LineSpec is not long enough!\n");
   } else {
      LineSpec[e2i( 13940)].in += inten * 0.012;
      LineSpec[e2i( 16910)].in += inten * 0.006;
      LineSpec[e2i( 17750)].in += inten * 0.012;
      LineSpec[e2i( 20780)].in += inten * 0.0034;
      LineSpec[e2i( 21214)].in += inten * 0.0023;
      LineSpec[e2i( 26345)].in += inten * 0.0231;
      LineSpec[e2i( 32183)].in += inten * 0.0359;
      LineSpec[e2i( 33196)].in += inten * 0.001215;
      LineSpec[e2i( 59541)].in += inten * 0.3592;
   }
}
//------------------------------------------------------------------------------
void AddFe055(spec_t* LineSpec, int length, double inten)
{ 
   if (length-1 < e2i(6535)) {
      printf("ERROR: length of LineSpec is not long enough!\n");
   } else {
      LineSpec[e2i(  5888)].in += inten * 0.0845;
      LineSpec[e2i(  5899)].in += inten * 0.1657;
      LineSpec[e2i(  6490)].in += inten * 0.0340;
      LineSpec[e2i(  6535)].in += inten * 0.0340;
   }
}
//------------------------------------------------------------------------------
void AddCo057(spec_t* LineSpec, int length, double inten)
{
   if (length-1 < e2i(692010)) {
      printf("ERROR: length of LineSpec is not long enough!\n");
   } else {
      LineSpec[e2i(  6391)].in += inten * 0.168;
      LineSpec[e2i(  6404)].in += inten * 0.332;
      LineSpec[e2i(  7058)].in += inten * 0.071;
      LineSpec[e2i(  7108)].in += inten * 0.071;
      LineSpec[e2i( 14413)].in += inten * 0.0915;
      LineSpec[e2i(122061)].in += inten * 0.8551;
      LineSpec[e2i(136474)].in += inten * 0.1071;
      LineSpec[e2i(692010)].in += inten * 0.00159;
   }
}
//------------------------------------------------------------------------------
void ApplyPileUp(spec_t* DetSpec, spec_t* PUSpec, int length, config_t* c)
{
   int i_E, ii_E;                              //index for DetSpec
   double la, la_0 = c->t_out/1000 * c->f_in;   //lambda for intensity = 1.0
   printf("lambda_0 = %E\n", la_0);
   for (i_E=0; i_E<=length-1; ++i_E) {
      if (DetSpec[i_E].in > TH) {
         for (ii_E=0; ii_E<=length-1; ++ii_E) {
            if (DetSpec[ii_E].in > TH) {
               if ((i_E + ii_E) <= length-1) {
                  la = la_0 * DetSpec[ii_E].in;
                  PUSpec[i_E + ii_E].in += DetSpec[i_E].in * la * exp(-la);
               }
            }      
         }
      }
   }
}
//------------------------------------------------------------------------------
void WriteOutput(spec_t* DetSpec, int length)
{
   int i;
   FILE *fpout;
   fpout = fopen("SynSpec_out.txt", "w");
   if (fpout==NULL) {
      printf("ERROR: cannot open SynSpec_out.txt");
   } else {
      fprintf(fpout, "#energy[eV]  rel.intensity\n");
      for (i=0; i<=length-1; ++i) {
         fprintf(fpout, "%f    %E\n", DetSpec[i].erg/1000.0, DetSpec[i].in);
      }
   fclose(fpout);
   }
}
//------------------------------------------------------------------------------
int e2i(double erg)
{
   return (int)(erg/binning);
}
//------------------------------------------------------------------------------
double i2e(int i)
{
   return (i+0.5) * binning;
}
//------------------------------------------------------------------------------
double myabs(double a)
{
   if (a >= 0) {
      return a;
   } else {
      return -a;
   }
}
