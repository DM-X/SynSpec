#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "SynSpec.h"

#define TH       0.00000000001
#define NBR_FLUO 5
float binning;  //in eV; set by configuration file
void AddAm241(spec_t* LineSpec, int length, double inten);

int main(int argc, char* argv[])
{
   printf("###########################################################\n");
   printf("###########################################################\n");
   printf("##                          SynSpec                      ##\n");
   printf("###########################################################\n");
   printf("###############################################24.02.2018##\n");
   if (argc==1) {
      printf("SynSpec is a tool to create a syntetic X-ray spectrum.\n");
      printf("All configuration are set by the SynSpec_config.txt file.\n");
      printf("To plot sepctrum use: gnuplot SynSpec.plt Synspec.dat\n");
   }
   
   config_t c;
   MakeConfig("SynSpec_config.txt", &c);
   InitConfiguration(&c);
   ReadConfiguration("SynSpec_config.txt", &c);
   PrintConfiguration(&c);
   
   spec_t *LineSpec, *DetSpec, *ComSpec, *PUSpec;
   int length = ceilf((c.bin_e - c.bin_s) / c.bin_w -0.000001);
   MakeSpecs(&LineSpec, &DetSpec, &ComSpec, &PUSpec, length);
   InitSpecs(LineSpec, DetSpec, ComSpec, PUSpec, length);
 
   int i, udr_cnt=0, ovr_cnt=0;
   for (i=0; i<=c.nbr_src-1; ++i) {
      AddSource(c.src_name[i],c.src_int[i],LineSpec, length, 
                &udr_cnt, &ovr_cnt, c.use_oor);
   }
   //AddAm241(LineSpec, length, 1);

  // printf("line spec: %f\n", LineSpec[e2i(1332500)].in);
  // //LineSpec[e2i(  60000)].in += 0.09; //testing...

  AddFluorescence(LineSpec, length, c);

  ApplyAbsorber(LineSpec, length, &c);

  ApplyEfficiency_PE(LineSpec, length, c);

  // //Comptonization is not jet added to LineSpec in order to keep
  // //the spectrum non-continuous
  // CalcCompt(LineSpec, ComSpec, length, c);
  // CalcCBS(LineSpec, ComSpec, length, c);
  // //printf("Comptonization calculated.\n");
  // 
  //ApplyEscape(LineSpec, length, c);

  // ApplyTail(LineSpec, ComSpec, length, &c);
  // AddSpectra(LineSpec, ComSpec, length);

  // //ApplyPileUp(LineSpec, PUSpec, length, &c);
  // //AddSpectra(LineSpec, PUSpec, length);

   ApplyLT(LineSpec, length, &c);
   ApplyDetRes(LineSpec, DetSpec, length, c);
  // printf("detector resolution applied.\n");

   Normalize(DetSpec, length, c.norm);

   WriteOutput(DetSpec, length);
   //WriteOutput(LineSpec, length);

   FreeSpecs(&LineSpec, &DetSpec, &ComSpec, &PUSpec);
   FreeConfig(&c);
   return 0;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void  InitConfiguration(config_t* c)
{
   int i;
   strcpy(c->output, "txt"); 
   for (i=0; i<=c->nbr_src-1; ++i) {
      strcpy(c->src_name[i], "");
      c->src_int[i] = 0.0;
   }
   for (i=0; i<=c->nbr_abs-1; ++i) {
      strcpy(c->abs_mat[i], "");
      c->abs_thk[i] = 0.0;
   }
   for (i=0; i<=c->nbr_fluo-1; ++i) {
      strcpy(c->fluo_mat[i], "");
      c->fluo_int[i] = 0.0;
   }
   strcpy(c->det_mat, "");
   c->det_d    = 0.0;
   c->det_FM   = 1.0;
   c->det_res0 = 100;
   c->det_LT   = 0;
   c->t_out    = 100;
   c->f_in     =   1;
   c->bin_w    = 100;
   c->use_oor  = 0;
}
//------------------------------------------------------------------------------
void MakeConfig(char* file_name, config_t* c)
{
   FILE* fpin_config;
   fpin_config = fopen(file_name, "r");
   if (fpin_config == NULL) {
      printf("ERROR: cannot open %s\n", file_name);
   }

   char line[256];
   char key[32];
   char src_key[32], abs_key[32], fluo_key[32];
   int  i_src=0, i_abs=0, i_fluo=0;
   while (fgets(line, sizeof(line), fpin_config)) {
      if (line[0]!='#') {
         sscanf(line, "%s =", key);

         sprintf(src_key,  "src_name_%d", i_src);
         sprintf(abs_key,  "abs_mat_%d",  i_abs);
         sprintf(fluo_key, "fluo_mat_%d", i_fluo);
         if (!strcmp(key, src_key))  ++i_src;
         if (!strcmp(key, abs_key))  ++i_abs;
         if (!strcmp(key, fluo_key)) ++i_fluo;
      }
   }
   printf("i_src = %d\n", i_src);
   printf("i_abs = %d\n", i_abs);
   printf("i_flu = %d\n", i_fluo);

   c->nbr_src  = i_src;
   c->nbr_abs  = i_abs;
   c->nbr_fluo = i_fluo;

   int i;
   c->src_name = (char**)malloc(i_src * sizeof(char*));
   if (c->src_name==NULL) {printf("ERROR: alloc. src_name\n"); exit(1);}
   for (i=0; i<=i_src-1; ++i) {
      c->src_name[i] = (char*)malloc(8 * sizeof(char));
      if (c->src_name[i]==NULL) {printf("ERROR: alloc. src_name\n"); exit(1);}
   }
   c->src_int = (float*)malloc(i_src * sizeof(float));
   if (c->src_int==NULL) {printf("ERROR: alloc. src_int\n"); exit(1);}

   c->abs_mat = (char**)malloc(i_abs * sizeof(char*));
   if (c->abs_mat==NULL) {printf("ERROR: alloc. abs_mat\n"); exit(1);}
   for (i=0; i<=i_abs-1; ++i) {
      c->abs_mat[i] = (char*)malloc(8 * sizeof(char));
      if (c->abs_mat[i]==NULL) {printf("ERROR: alloc. abs_mat\n"); exit(1);}
   }
   c->abs_thk = (float*)malloc(i_abs * sizeof(float));
   if (c->abs_thk==NULL) {printf("ERROR: alloc. abs_thk\n"); exit(1);}
   
   c->fluo_mat = (char**)malloc(i_fluo * sizeof(char*));
   if (c->fluo_mat==NULL) {printf("ERROR: alloc. fluo_mat\n"); exit(1);}
   for (i=0; i<=i_fluo-1; ++i) {
      c->fluo_mat[i] = (char*)malloc(8 * sizeof(char));
      if (c->fluo_mat[i]==NULL) {printf("ERROR: alloc. fluo_mat\n"); exit(1);}
   }
   c->fluo_int = (float*)malloc(i_fluo * sizeof(float));
   if (c->fluo_int==NULL) {printf("ERROR: alloc. fluo_int\n"); exit(1);}

   fclose(fpin_config);
}
//------------------------------------------------------------------------------
void ReadConfiguration(char* file_name, config_t* c)
{
   int i_src_name=0, i_src_int=0;
   int i_abs_mat =0, i_abs_thk=0;
   int i_fluo_mat=0, i_fluo_int=0;
   FILE* fpin_config;
   fpin_config = fopen(file_name, "r");
   if (fpin_config==NULL) {printf("ERROR: opening %s\n", file_name); exit(1);}
   
   printf("## reading config file...");
   char  line[256];
   char  key[32], cmp_key[32];

   while (fgets(line, sizeof(line), fpin_config)) {
      if (line[0]!='#') {
         sscanf(line, " output    = %s", c->output);

         sscanf(line, "%s =", key);
         sprintf(cmp_key,  "src_name_%d", i_src_name);
         if (!strcmp(key, cmp_key)) {
            sscanf(line, "%*s = %s", c->src_name[i_src_name]);
            ++i_src_name;
         }
         sprintf(cmp_key,  "src_int_%d", i_src_int);
         if (!strcmp(key, cmp_key)) {
            sscanf(line, "%*s = %f", &(c->src_int[i_src_int]));
            ++i_src_int;
         }

         sprintf(cmp_key,  "abs_mat_%d",  i_abs_mat);
         if (!strcmp(key, cmp_key)) {
            sscanf(line, "%*s = %s", c->abs_mat[i_abs_mat]);
            ++i_abs_mat;
         }
         sprintf(cmp_key,  "abs_thk_%d",  i_abs_thk);
         if (!strcmp(key, cmp_key)) {
            sscanf(line, "%*s = %f", &(c->abs_thk[i_abs_thk]));
            ++i_abs_thk;
         }

         sprintf(cmp_key, "fluo_mat_%d", i_fluo_mat);
         if (!strcmp(key, cmp_key)) {
            sscanf(line, "%*s = %s", c->fluo_mat[i_fluo_mat]);
            ++i_fluo_mat;
         }
         sprintf(cmp_key, "fluo_int_%d", i_fluo_int);
         if (!strcmp(key, cmp_key)) {
            sscanf(line, "%*s = %f", &(c->fluo_int[i_fluo_int]));
            ++i_fluo_int;
         }

         sscanf(line, " CdTe_absd = %f", &((*c).CdTe_absd));
         sscanf(line, " Si_absd   = %f", &((*c).Si_absd));
         sscanf(line, " det_mat   = %s", (*c).det_mat);
         sscanf(line, " det_d     = %f", &((*c).det_d));
         sscanf(line, " det_FM    = %f", &((*c).det_FM));
         sscanf(line, " det_res0  = %f", &((*c).det_res0));
         sscanf(line, " det_LT    = %f", &((*c).det_LT)); 
         sscanf(line, " depl_volt = %f", &((*c).depl_volt));
         sscanf(line, " mu        = %f", &((*c).mu));
         sscanf(line, " tau       = %f", &((*c).tau));
         sscanf(line, " t_out     = %f", &((*c).t_out));
         sscanf(line, " f_in      = %f", &((*c).f_in));
         sscanf(line, " bin_start = %f", &((*c).bin_s));
         sscanf(line, " bin_end   = %f", &((*c).bin_e));
         sscanf(line, " bin_width = %f", &((*c).bin_w));
         sscanf(line, " norm      = %f", &((*c).norm));
         sscanf(line, " use_oor   = %d", &((*c).use_oor));
      }
   }
   c->bin_s  *= 1000.0; //keV -> eV
   c->bin_e  *= 1000.0;
   c->bin_w  *= 1000.0;
   c->det_LT *= 1000.0;

   if (i_src_name != i_src_int  || i_src_name != c->nbr_src ||
       i_abs_mat  != i_abs_thk  || i_abs_mat  != c->nbr_abs ||
       i_fluo_mat != i_fluo_int || i_fluo_mat != c->nbr_fluo) {
      printf("ERROR: src, abs, or fluo input is inconsistent!\n");
      exit(1);
   }

   binning = c->bin_w;
   printf("done\n");
   fclose(fpin_config);
}
//------------------------------------------------------------------------------
void PrintConfiguration(config_t* c)
{
   int i;
   printf("## SOURCES: %d\n", c->nbr_src);
   for (i=0; i<=c->nbr_src-1; ++i) {
      printf("   [%d]  %s, int = %f\n", i, c->src_name[i], c->src_int[i]);
   }
   printf("## ABSORBER: %d\n", c->nbr_abs);
   for (i=0; i<=c->nbr_abs-1; ++i) {
      printf("   [%d]  %s, thk = %f\n", i, c->abs_mat[i], c->abs_thk[i]);
   }
   printf("## FLUORESCENCE SOURCES: %d\n", c->nbr_fluo);
   for (i=0; i<=c->nbr_fluo-1; ++i) {
      printf("   [%d]  %s, int = %f\n", i, c->fluo_mat[i], c->fluo_int[i]);
   }
   printf("## OUTPUT: %s\n", c->output);
   printf("## DETECTOR:\n");
   printf("   %s     d=%5.3f mm, res = %5.3f * F + %5.3f eV FWHM, LT=%f keV\n", 
          c->det_mat, c->det_d, c->det_FM, c->det_res0, c->det_LT);
   printf("   depl = %5.1f V,  mu = %5.1f cm^2/(Vs),  tau =%5.1f us\n",
          c->depl_volt, c->mu, c->tau);
   printf("   t_out = %5.3fms,  f_in = %5.3f Hz/ch\n", c->t_out, c->f_in);
   printf("## CALCULATION:\n");
   printf("   bin_s = %5.1f, bin_e = %5.1f, bin_w = %5.1f\n", 
          c->bin_s, c->bin_e, c->bin_w);
   printf("   use_oor = %d\n", c->use_oor);
}
//------------------------------------------------------------------------------
void  FreeConfig(config_t* c)
{
   int i;
   for (i=0; i<=c->nbr_src-1; ++i) {
      free(c->src_name[i]);
      c->src_name[i] = NULL;
   }
   free(c->src_name); 
   c->src_name = NULL;

   for (i=0; i<=c->nbr_abs-1; ++i) {
      free(c->abs_mat[i]);
      c->abs_mat[i] = NULL;
   }
   free(c->abs_mat);
   c->abs_mat = NULL;

   for (i=0; i<=c->nbr_fluo-1; ++i) {
      free(c->fluo_mat[i]);
      c->fluo_mat[i] = NULL;
   }
   free(c->fluo_mat);
   c->fluo_mat = NULL;

   free(c->src_int); 
   c->src_int  = NULL;
   free(c->abs_thk);
   c->abs_thk  = NULL;
   free(c->fluo_int);
   c->fluo_int = NULL;
}
//------------------------------------------------------------------------------
void MakeSpecs(spec_t** LineSpec, spec_t** DetSpec, spec_t** ComSpec, 
               spec_t** PUSpec, int length)
{
   *LineSpec = (spec_t*)malloc(length * sizeof(spec_t));
   *DetSpec  = (spec_t*)malloc(length * sizeof(spec_t));
   *ComSpec  = (spec_t*)malloc(length * sizeof(spec_t));
   *PUSpec   = (spec_t*)malloc(length * sizeof(spec_t));
   if (*LineSpec == NULL) {printf("ERROR: cannot alloc LineSpec!\n"); exit(1);}
   if (*DetSpec  == NULL) {printf("ERROR: cannot alloc DetSpec!\n"); exit(1);}
   if (*ComSpec  == NULL) {printf("ERROR: cannot alloc ComSpec!\n"); exit(1);}
   if (*PUSpec   == NULL) {printf("ERROR: cannot alloc PUSpec!\n"); exit(1);}
}
//------------------------------------------------------------------------------
void InitSpecs(spec_t* LineSpec, spec_t* DetSpec, spec_t* ComSpec, 
               spec_t* PUSpec, int length)
{
   int i;
   float energy;
   for (i=0; i<=length-1; ++i) {
      LineSpec[i].in  = 0;
      DetSpec[i].in   = 0;
      ComSpec[i].in   = 0;
      PUSpec[i].in    = 0;
      energy = i2e(i);
      LineSpec[i].erg = energy;
      DetSpec[i].erg  = energy;
      ComSpec[i].erg  = energy;
      PUSpec[i].erg   = energy;
   }
}
//------------------------------------------------------------------------------
void FreeSpecs(spec_t** LineSpec, spec_t** DetSpec, spec_t** ComSpec, 
               spec_t** PUSpec)
{
   free(*LineSpec);
   *LineSpec=NULL;
   free(*DetSpec);
   *DetSpec=NULL;
   free(*ComSpec);
   *ComSpec=NULL;
   free(*PUSpec);
   *PUSpec=NULL;
}
//------------------------------------------------------------------------------
void ApplyAbsorber(spec_t* LineSpec, int length, config_t* c)
{
   int i, i_E; //index for LineSpec
   for (i=0; i<c->nbr_abs-1; ++i) {
      if (c->abs_thk[i] > TH) {
         for (i_E=0; i_E<= length-1; ++i_E) {
            if(LineSpec[i_E].in > TH) {
               //printf("i_E=%d  E = %f  -> %f\n", 
                //       i_E, LineSpec[i_E].erg, LineSpec[i_E].in);
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
   int i, i_E; //index for LineSpec
   int i_fluo; //index for fluorescence line 
   for (i_E=0; i_E<= length-1; ++i_E) {
      if(LineSpec[i_E].in > TH) {
         for (i=0; i<NBR_FLUO; ++i) {
            if (!strcmp(c.fluo_mat[i], "C_g")) {
               i_fluo = e2i(277);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += c.fluo_int[i];}
               continue;
            }
            if (!strcmp(c.fluo_mat[i], "C_d")) {
               i_fluo = e2i(277);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += c.fluo_int[i];}
               continue;
            }
            if (!strcmp(c.fluo_mat[i], "Al")) {
               i_fluo = e2i(1487);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += c.fluo_int[i];}
               continue;
            }
            if (!strcmp(c.fluo_mat[i], "Si")) {
               i_fluo = e2i(1740);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += c.fluo_int[i];}
               continue;
            }
            if (!strcmp(c.fluo_mat[i], "Ti")) {
               i_fluo = e2i(4511);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += c.fluo_int[i];}
               continue;
            }
            if (!strcmp(c.fluo_mat[i], "Fe")) {
               i_fluo = e2i(6406);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += c.fluo_int[i];}
               continue;
            }
            if (!strcmp(c.fluo_mat[i], "Cu")) {
               i_fluo = e2i(8048);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += c.fluo_int[i];}
               continue;
            }
            if (!strcmp(c.fluo_mat[i], "Ge")) {
               i_fluo = e2i(9886);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += c.fluo_int[i];}
               continue;
            }
            if (!strcmp(c.fluo_mat[i], "Ag")) {
               i_fluo = e2i(22163);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += c.fluo_int[i];}
               continue;
            }
            if (!strcmp(c.fluo_mat[i], "Cd")) {
               i_fluo = e2i(23174);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += c.fluo_int[i];}
               continue;
            }
            if (!strcmp(c.fluo_mat[i], "Sn")) {
               i_fluo = e2i(25271);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += c.fluo_int[i];}
               continue;
            }
            if (!strcmp(c.fluo_mat[i], "Te")) {
               i_fluo = e2i(27472);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += c.fluo_int[i];}
               continue;
            }
            if (!strcmp(c.fluo_mat[i], "W")) {
               i_fluo = e2i(59318);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += c.fluo_int[i];}
               continue;
            }
            if (!strcmp(c.fluo_mat[i], "Pt")) {
               i_fluo = e2i(66832);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += c.fluo_int[i];}
               continue;
            }
            if (!strcmp(c.fluo_mat[i], "Au")) {
               i_fluo = e2i(68804);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += c.fluo_int[i];}
               continue;
            }
            if (!strcmp(c.fluo_mat[i], "Pb")) {
               i_fluo = e2i(72804);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += 0.5*c.fluo_int[i];}
               i_fluo = e2i(74976);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += 1.0*c.fluo_int[i];}
               i_fluo = e2i(84936);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += 0.1*c.fluo_int[i];}
               i_fluo = e2i(10450);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += 0.1*c.fluo_int[i];}
               i_fluo = e2i(10552);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += 0.2*c.fluo_int[i];}
               i_fluo = e2i(12614);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += 0.1*c.fluo_int[i];}
               i_fluo = e2i(12623);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += 0.1*c.fluo_int[i];}
               i_fluo = e2i(14764);
               if (i_E > i_fluo) {LineSpec[i_fluo].in += 0.02*c.fluo_int[i];}
               continue;
            }
         }
      }
   }
}
////------------------------------------------------------------------------------
//void CalcCompt(spec_t* LineSpec, spec_t* ComSpec, int length, config_t c)
//{
//   int   i_E, ii_E;
//   float ComptEdge, LineIn, Erg;
//   for (i_E=0; i_E<=length-1; ++i_E) {
//      if (LineSpec[i_E].in > TH) {
//         LineIn    = LineSpec[i_E].in;
//         ComptEdge = i2e(i_E)*(1.0 - 1.0/(1.0 + (2.0*i2e(i_E)/511000.0)));
//         Erg = i2e(i_E);
//         for (ii_E=0; ii_E<=e2i(ComptEdge); ++ii_E) {
//            ComSpec[ii_E].in += LineIn * DsigIS_De(i2e(ii_E), Erg);
//         }
//      }
//   }
//}
////------------------------------------------------------------------------------
//void CalcCBS(spec_t* LineSpec, spec_t* ComSpec, int length, config_t c)
//{
//   int   i_E, ii_E, ii_E0;
//   float ComptEdge, LineIn, Erg;
//   for (i_E=0; i_E<=length-1; ++i_E) {
//      if (LineSpec[i_E].in > TH) {
//         LineIn    = LineSpec[i_E].in;
//         ComptEdge = i2e(i_E)*(1.0 - 1.0/(1.0 + (2.0*i2e(i_E)/511000.0)));
//         Erg = i2e(i_E);
//         ii_E0     = i_E - e2i(ComptEdge);
//         for (ii_E=ii_E0; ii_E<=i_E-1; ++ii_E) {
//            ComSpec[ii_E].in += LineIn * DsigIBS_De(i2e(ii_E), Erg);
//         }
//      }
//   }
//}
////------------------------------------------------------------------------------
//float DsigIS_De(float Erg, float k0)
//{
//   float k = k0 - Erg; //energy of sc. photon
//   return DsigCS_De(k, k0) * sin(theta(k, k0)) * theta_bar(k, k0) *
//          S("CdTe", theta(k, k0), k0);
//}
////------------------------------------------------------------------------------
//float DsigIBS_De(float Erg, float k0)
//{
//   float k = Erg; //energy of sc. photon
//   return DsigCS_De(k, k0) * sin(theta(k, k0)) * theta_bar(k, k0) *
//          S("CdTe", theta(k, k0), k0) * myabs(cos(theta(k, k0)));
//   //the additional cos(theta) accounts for the apparent detector area
//}
////------------------------------------------------------------------------------
//float DsigCS_De(float k, float k0)
//{
//   const double r_0_2 = 0.07740787; //sq of classical el. radius
//   const double kdk0  = k/k0;
//   return  M_PI * r_0_2 * (kdk0)*(kdk0) * 
//           (kdk0 + k0/k - (sin(theta(k, k0)))*(sin(theta(k, k0))));
//}
////------------------------------------------------------------------------------
//float S(char* element, float theta, float k0)
//{
//   const float k0_keV = k0 / 1000.0;
//   const float t_2    = theta / 2.0;
//   if (!strcmp(element, "Si"))
//      return 14.0 * (1-exp(-0.058942 * pow(k0_keV*sin(t_2), 1.1289))) / 
//                    (1-0.638*exp(-0.09432*k0_keV*sin(t_2)));
//   if (!strcmp(element, "Ge"))
//      return 32.0 * (1-exp(-0.039373 * pow(k0_keV * sin(t_2), 1.1002))) / 
//                    (1-0.570*exp(-0.05860*k0_keV*sin(t_2)));
//   if (!strcmp(element, "Sn") || !strcmp(element, "CdTe"))
//      return 50.0 * (1-exp(-0.025066 * pow(k0_keV * sin(t_2), 1.0714))) / 
//                    (1-0.702*exp(-0.03408*k0_keV*sin(t_2)));
//   return -1; //indicates error
//}
////------------------------------------------------------------------------------
//float theta(float k, float k0)
//{
//   return acos(1.0 - (511000.0*k0 / k-511000.0) / k0);
//}
////------------------------------------------------------------------------------
//float theta_bar(float k, float k0)
//{
//   return 1.0 / (sqrt(1.0 - (1.0 + 511000.0 * (k-k0)/(k*k0)) *
//                            (1.0 + 511000.0 * (k-k0)/(k*k0)))) * 511000.0/(k*k);
//}
////------------------------------------------------------------------------------
//void AddComptonization(spec_t* ComSpec, spec_t* DetSpec, int length)
//{
//   int i_D;
//   for (i_D=0; i_D<=length-1; ++i_D) {
//      DetSpec[i_D].in += ComSpec[i_D].in;
//   }
//}
////------------------------------------------------------------------------------
//float Mu_abs(int i_E, char* mat)
//{
//   char file_name[32];
//   sprintf(file_name, "SynSpec_material_%s.txt", mat);
//   FILE* fp_mat;
//   fp_mat = fopen(file_name, "r");
//   if (fp_mat==NULL) {
//      printf("ERROR: cannot open material definition file %s\n", file_name);
//   }
//
//   char  line[256];
//   float erg0=i2e(i_E)/1E6; //in MeV
//   float rho; //density in g/cm^3
//   float att; //in cm^2/g
//   float du;  //dump unused informations
//   float erg=0, last_erg=0;
//   float last_att=0;
//   while (fgets(line, sizeof(line), fp_mat)) {
//      if (line[0]!='#') {
//         sscanf(line, " density    = %f", &rho);
//         sscanf(line, " %f %f %f %f %f %f %f %f",
//                &erg, &du, &du, &att, &du, &du, &du, &du);
//      }
//      if ((last_erg < erg0) && (erg0 <= erg)) {//interpolate attenuation
//         att = last_att + (att-last_att)/(erg-last_erg)*(erg0-last_erg);
//         break;
//      }
//      last_att = att;
//      last_erg = erg;
//   }
// 
//   fclose(fp_mat);
//   return att*rho;
//}
//------------------------------------------------------------------------------
float Transmission(int i_E, float d, char* mat)
{
   char file_name[64];
   sprintf(file_name, "./materials/SynSpec_material_%s.txt", mat);
   FILE* fp_mat;
   fp_mat = fopen(file_name, "r");
   if (fp_mat==NULL) {
      printf("ERROR: cannot open material definition file %s\n", file_name);
   }

   char  line[256];
   float erg0=i2e(i_E)/1E6; //in MeV
   float rho; //density in g/cm^3
   float att; //in cm^2/g
   float du;  //dump unused informations
   float erg=0, last_erg=0;
   float last_att=0;
   while (fgets(line, sizeof(line), fp_mat)) {
      if (line[0]!='#') {
         sscanf(line, " density    = %f", &rho);
         sscanf(line, " %f %f %f %f %f %f %f %f",
                &erg, &du, &du, &att, &du, &du, &du, &du);
      }
      if ((last_erg < erg0) && (erg0 < erg)) {//interpolate attenuation
         att = last_att + (att-last_att)/(erg-last_erg)*(erg0-last_erg);
         break;
      }
      last_att = att;
      last_erg = erg;
   }
 
   float dcm = d/10.0; //thickness in cm
   fclose(fp_mat);
   //printf("i_E = %d, E0=%f, E=%f, att=%f, trans=%f\n", 
    //       i_E,      erg0,   erg, att,    exp(-att*rho*dcm));
   return exp(-att*rho*dcm);
}
//------------------------------------------------------------------------------
void ApplyEscape(spec_t* LineSpec, int length, config_t c)
{
   int   i_E, j, i_E0 = length;
   int   i_E_edg[8];           //edge energy for corresponding shell(index)
   int   i_E_esc[8];           //escaping energies arranged in order of size
   float P_esc[8];             //escape probabilities
   if (!strcmp(c.det_mat, "Si")) {
      i_E_edg[0] = e2i(1839);      
      i_E_esc[0] = e2i(1740);  //Si-Ka
      P_esc[0]   = 0.25*0.05;  //fluo yield 0.05. todo: precise the 0.25
      i_E_esc[1] = -1;         //indicates end of escape list
      i_E0 = i_E_edg[0];
      printf("apply Si escape...");
   }
   if (!strcmp(c.det_mat, "Ge")) {
      i_E_edg[0] = e2i(11103);
      i_E_esc[0] = e2i(9886);  //Ge-Ka
      P_esc[0]   = 0.25*0.52;  //fluo yield 0.52. todo: precise the 0.25
      i_E_esc[1] = -1;         //indicates end of escape list
      i_E0 = i_E_edg[0];
      printf("apply Ge escape...");
   }
   if (!strcmp(c.det_mat, "CdTe")) { 
      i_E_edg[0] = e2i(26711);
      i_E_esc[0] = e2i(23174); //Cd-Ka
      P_esc[0]   = 0.07*0.84;  //fluo yield 0.84. todo: precise the 0.07
                               //idea was: 0.25 (out of thesis) / 2 because of
                               //the chance to hit Cd or Te
      i_E_edg[1] = e2i(31814);  
      i_E_esc[1] = e2i(27472); //Te-Ka
      P_esc[1]   = 0.03*0.88;  //fluo yield 0.88. todo: precise the 0.03
                               //idea was: additional absorption due to Cd 
                               //absorption (Te_Ka > Cd_Kedge)
      i_E_esc[2] = -1;         //indicates end of escape list
      i_E0 = i_E_edg[0];
      printf("apply Cd and Te escape...");
   }
   if (i_E0 >= length) {
      printf("WARNING: escape peak out of range: %d>=%d\n", i_E0, length);
   }
   for (i_E=i_E0; i_E<=length-1; ++i_E) {
      if(LineSpec[i_E].in > TH) {
         for (j=0; j<=7; ++j) {
            if (i_E_esc[j]==-1) {break;}
            if (i_E >= i_E_edg[j]) {
               LineSpec[i_E-i_E_esc[j]].in += P_esc[j]*LineSpec[i_E].in;
               LineSpec[i_E].in            -= P_esc[j]*LineSpec[i_E].in;
            }
         }
      }
   }
   printf("done %s\n", c.det_mat);
}
//------------------------------------------------------------------------------
void ApplyEfficiency_PE(spec_t* LineSpec, int length, config_t c)
{
   int i_E;
   for (i_E=0; i_E<=length-1; ++i_E) {
      if(LineSpec[i_E].in > TH) {
         LineSpec[i_E].in *= Efficiency(i_E, c.det_d, c.det_mat); 
      }
   }
}
//------------------------------------------------------------------------------
void ApplyLT(spec_t* LineSpec, int length, const config_t* c)
{
   int i_E; 
   int i_end = e2i(c->det_LT);
   printf("deb: i_end = %d (%d)\n", i_end, length);
   if (i_end >= length-1) i_end = length-1;
   for (i_E=0; i_E<=i_end; ++i_E) {
      LineSpec[i_E].in = 0.0;
   }
}
//------------------------------------------------------------------------------
void ApplyDetRes(spec_t* LineSpec, spec_t* DetSpec, int length, config_t c)
{
   int i_E, i_D, i_D_start, i_D_stop; //index for LineSpec and DetSpec
   float var2, omega, F;
   float res0 = pow(c.det_res0/2.3548, 4);
   float root;
   float x, x0;
   GetDetParameter(&omega, &F, &c);
   float FOmega = c.det_FM*c.det_FM * F * omega;
   
   for (i_E=0; i_E<=length-1; ++i_E) {
      if(LineSpec[i_E].in > TH) {
         x0 = i2e(i_E);
         var2 = 2.0 * sqrt((x0 * FOmega) * (x0 * FOmega) + res0);
         root  = sqrt(M_PI*var2);
         i_D_start = i_E - e2i(3*root);
         if (i_D_start < 0) i_D_start = 0;
         i_D_stop  = i_E + e2i(3*root);
         if (i_D_stop >= length) i_D_stop = length-1;

         for (i_D=i_D_start; i_D<=i_D_stop; ++i_D) {
            x = i2e(i_D);
            DetSpec[i_D].in += LineSpec[i_E].in / root * 
                               exp(-(x-x0)*(x-x0) / var2);
         }
      }
   } 
}
//------------------------------------------------------------------------------
void GetDetParameter(float* omega, float* F, config_t* c)
{
   if (!strcmp(c->det_mat, "Si")) {
      *omega = 3.63;
      *F = 0.115;
      return;
   }
   if (!strcmp(c->det_mat, "Ge")) {
      *omega = 2.9;
      *F = 0.13;
      return;
   }
   if (!strcmp(c->det_mat, "CdTe")) {
      *omega = 4.43;
      *F = 0.15;
      return;
   }
}
//------------------------------------------------------------------------------
//void AddSpectra(spec_t* DetSpec, spec_t* ComSpec, int length)
//{
//   int i_E; //index for LineSpec
//   for (i_E=0; i_E<=length-1; ++i_E) DetSpec[i_E].in += ComSpec[i_E].in;
//}
//------------------------------------------------------------------------------
float Efficiency(int energy, float d, char* mat)
{
   return 1 - Transmission(energy, d, mat);
}
//------------------------------------------------------------------------------
//float myGauss(float x0, float x, float height, config_t c)
//{
//   float var, omega=1, F=1;
//   if (!strcmp(c.det_mat, "Si")) {
//      omega = 3.63;
//      F = 0.10;
//   } else {
//      if (!strcmp(c.det_mat, "CdTe")) {
//         omega = 4.43;
//         F = 0.15;
//      } else {
//         printf("ERROR: unknown detector material: %s\n", c.det_mat);
//      }
//   }
//   var = sqrt((c.det_FM * F * x0 * omega) * (c.det_FM * F * x0 * omega) + 
//              pow(c.det_res0/(2.3548*omega), 4));
//   return height * exp(-(x-x0)*(x-x0) / (2*var));
//}
////------------------------------------------------------------------------------
//void ApplyTail(spec_t* LineSpec, spec_t* ComSpec, int length, config_t* c)
//{
//   float L = c->mu * c->tau/1E6 * c->depl_volt / (c->det_d/10.0);
//   float h_min = 1-(c->det_d/10.0)/(2*L);
//   float h;         //relative pulse height
//   float s, sum=0;  //calc sum of all smeared signals for peak reduction 
//   float in_0, mu_abs;
//   int i_E, ii_E;
//   for (i_E=0; i_E<=length-1; ++i_E) {
//      if(LineSpec[i_E].in > TH) {
//         in_0 = LineSpec[i_E].in;
//         mu_abs = Mu_abs(i_E, c->det_mat);
//         printf("mu_abs = %f, L=%f, h_min = %f\n", mu_abs, L, h_min);
//         for (ii_E=(int)(h_min*i_E); ii_E<=i_E; ++ii_E) {
//            h = (float)ii_E/i_E;
//            s = exp(-mu_abs * sqrt((1.0-h)*2.0*c->det_d*L));
//            sum += s;
//         }
//         //to optimze :run to i_E-1 in prev loop and do here: sum +=1
//         for (ii_E=(int)(h_min*i_E); ii_E<=i_E-1; ++ii_E) {
//            h = (float)ii_E/i_E;
//            s = exp(-mu_abs * sqrt((1.0-h)*2.0*c->det_d*L));
//            LineSpec[ii_E].in += s/sum * in_0;
//         }
//         printf("sum = %f\n", sum);
//         LineSpec[i_E].in = 1.0/sum * in_0;
//      }
//   }
//}
//------------------------------------------------------------------------------
void Normalize(spec_t* DetSpec, int length, float norm)
{
   int i;
   if (norm < 0.0) {
      float max_val=0.0;
      for (i=0; i<=length-1; ++i) {
         if (DetSpec[i].in > max_val) {
            max_val = DetSpec[i].in;
         }
      } 
      if (max_val==0) {
         printf("ERROR: no counts in detector spectrum!\n");
         exit(1);
      }
      for (i=0; i<=length-1; ++i) {
         DetSpec[i].in = DetSpec[i].in/max_val;
      }
      return;
   }
   if (norm==0.0) {
      float sum_val=0.0;
      for (i=0; i<=length-1; ++i) {
         sum_val += DetSpec[i].in;
      } 
      if (sum_val==0.0) {
         printf("ERROR: no counts in detector spectrum!\n");
         exit(1);
      }
      for (i=0; i<=length-1; ++i) {
         DetSpec[i].in = DetSpec[i].in/sum_val;
      }
      return;
   }
   if (norm>0.0) {
      float max_val=0.0;
      for (i=0; i<=length-1; ++i) {
         if (DetSpec[i].in > max_val) {
            max_val = DetSpec[i].in;
         }
      } 
      if (max_val==0) {
         printf("ERROR: no counts in detector spectrum!\n");
         exit(1);
      }
      for (i=0; i<=length-1; ++i) {
         DetSpec[i].in = DetSpec[i].in/max_val * norm;
      }
      return;
   }
}
//------------------------------------------------------------------------------
void AddSource(char* src, float src_int, spec_t* LineSpec, int length, 
               int*udr_cnt, int* ovr_cnt, int use_oor)
{
   char file_name[64];
   sprintf(file_name, "./sources/SynSpec_source_%s.txt", src);
   FILE* fp_src;
   fp_src = fopen(file_name, "r");
   if (fp_src==NULL) {
      printf("ERROR: cannot open source definition file %s\n", file_name);
   }

   char  line[256];
   char  name[256];
   char  data[256];
   char  type[16];
   float erg, line_int; //energy, line_intensity
   int   index;
   while (fgets(line, sizeof(line), fp_src)) {
      sscanf(line, "%s ; %s", name, data);
      if (!strcmp(name, "Reference")) {
         rewind(fp_src);
         break;
      }
   }
   while (fgets(line, sizeof(line), fp_src)) {
      if (line[0]!='#') {
         if (sscanf(line, "%f", &erg)) { //if first part of line is a number
            if (!strcmp(data, "Lara")) {
               sscanf(line, "%f %*s %f %*s", &erg, &line_int);
            } else {
               sscanf(line, "%f ;%*[^;]; %f ;%*[^;]; %s ;",
                      &erg, &line_int, type);
               if (!strcmp(type, "a")) continue;
            }
            erg      *= 1000; //in eV
            line_int /= 100;  //0..1.0
            index = e2i(erg);
            if (index >=0 && index <=length-1) {
               LineSpec[index].in += line_int * src_int;
            } else {
               if (use_oor) {
                  if (index < 0) {
                     LineSpec[0].in += line_int * src_int;
                     ++udr_cnt;
                  } else {
                     LineSpec[length-1].in += line_int * src_int;
                     ++ovr_cnt;
                  }
               } //else: don't use events that are out of range
            }
         }
      }
   }
 
   fclose(fp_src);
}
//------------------------------------------------------------------------------
//void ApplyPileUp(spec_t* DetSpec, spec_t* PUSpec, int length, config_t* c)
//{
//   int i_E, ii_E;                              //index for DetSpec
//   float la, la_0 = c->t_out/1000 * c->f_in;   //lambda for intensity = 1.0
//   //printf("lambda_0 = %E\n", la_0);
//   for (i_E=0; i_E<=length-1; ++i_E) {
//      if (DetSpec[i_E].in > TH) {
//         for (ii_E=0; ii_E<=length-1; ++ii_E) {
//            if (DetSpec[ii_E].in > TH) {
//               if ((i_E + ii_E) <= length-1) {
//                  la = la_0 * DetSpec[ii_E].in;
//                  PUSpec[i_E + ii_E].in += DetSpec[i_E].in * la * exp(-la);
//               }
//            }      
//         }
//      }
//   }
//}
//------------------------------------------------------------------------------
void WriteOutput(spec_t* DetSpec, int length)
{
   int i;
   FILE *fpout;
   fpout = fopen("SynSpec.dat", "w");
   if (fpout==NULL) {
      printf("ERROR: cannot open SynSpec.dat");
   } else {
      fprintf(fpout, "#energy[eV]  rel.intensity\n");
      for (i=0; i<=length-1; ++i) {
         fprintf(fpout, "%f    %E\n", DetSpec[i].erg/1000.0, DetSpec[i].in);
      }
   fclose(fpout);
   }
}
//------------------------------------------------------------------------------
int e2i(float erg)
{
   return (int)(erg/binning);
}
//------------------------------------------------------------------------------
float i2e(int i)
{
   return (i+0.5) * binning;
}
//------------------------------------------------------------------------------
void AddAm241(spec_t* LineSpec, int length, double inten)
{
   if (length-1 < e2i(59541)) {
      printf("ERROR: length of LineSpec is not long enough!\n");
   } else {
      LineSpec[e2i( 13940)].in += inten * 0.07;
      LineSpec[e2i( 16910)].in += inten * 0.08;
      LineSpec[e2i( 17750)].in += inten * 0.15;
      LineSpec[e2i( 20780)].in += inten * 0.05;
      LineSpec[e2i( 21214)].in += inten * 0.05;
      LineSpec[e2i( 26345)].in += inten * 0.11;
      LineSpec[e2i( 32183)].in += inten * 0.0359;
      LineSpec[e2i( 33196)].in += inten * 0.035;
      LineSpec[e2i( 36368)].in += inten * 0.05;
      LineSpec[e2i( 59541)].in += inten * 1.0;
   }
}
//------------------------------------------------------------------------------
//float myabs(float a)
//{
//   if (a >= 0) {
//      return a;
//   } else {
//      return -a;
//   }
//}
