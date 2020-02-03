typedef struct {
   char   output[8]; 
   int    nbr_src;
   char** src_name;
   float* src_int;
   float  CdTe_absd, Si_absd;
   int    nbr_abs;
   char** abs_mat;
   float* abs_thk;
   int    nbr_fluo;
   char** fluo_mat;
   float* fluo_int;
   char   det_mat[8];
   float  det_d, det_FM, det_res0;
   float  depl_volt, mu, tau;
   float  t_out, f_in;
   float  bin_s, bin_e, bin_w;
   float  LT;
   float  norm;
   float  max;
} config_t;

typedef struct {
   float erg;
   float in;
} spec_t;


void  MakeConfig(char* file_name, config_t* c);
void  InitConfiguration(config_t* c); 
void  ReadConfiguration(char* file_name, config_t* c);
void  PrintConfiguration(config_t* c);
void  FreeConfig(config_t* c);
void  MakeSpecs(spec_t** LineSpec, spec_t** DetSpec, spec_t** ComSpec, 
                spec_t** PUSpec, int length);
void  InitSpecs(spec_t* LineSpec, spec_t* DetSpec, spec_t* ComSpec, 
                spec_t* PUSpec, int length);
void  FreeSpecs(spec_t** LineSpec, spec_t** DetSpec, spec_t** ComSpec, 
                spec_t** PUSpec);

int   e2i(float erg);
float i2e(int i);

void  AddSource(char* src, float src_int, spec_t* LineSpec, int length, 
               int*udr_cnt, int* ovr_cnt, float LT);
void  GetDetParameter(float* omega, float* F, config_t* c);
void  ApplyAbsorber(spec_t* LineSpec, int length, config_t* c);
void  AddFluorescence(spec_t* LineSpec, int length, config_t c);
//void  CalcCompt(spec_t* LineSpec, spec_t* ComSpec, int length, config_t c);
//void  CalcCBS(spec_t* LineSpec, spec_t* ComSpec, int length, config_t c);
//float DsigIS_De(float Erg, float k0);
//float DsigIBS_De(float Erg, float k0);
//float DsigCS_De(float k, float k0);
//float S(char* element, float theta, float k0);
//float theta(float k, float k0);
//float theta_bar(float k, float k0);
void  ApplyEfficiency_PE(spec_t* LineSpec, int length, config_t c);
float Transmission(int i_E, float d, char* material);
//float Mu_abs(int i_E, char* mat);
//void  ApplyTail(spec_t* LineSpec, spec_t* ComSpec, int length, config_t* c);
void  ApplyDetRes(spec_t* LineSpec, spec_t* DetSpec, int length, config_t c);
//void  AddSpectra(spec_t* DetSpec, spec_t* ComSpec, int length);
float Efficiency(int energy, float d, char* mat);
//float myGauss(float x0, float x, float height, config_t c);
void Normalize(spec_t* DetSpec, int length, float norm, float max);
void  ApplyEscape(spec_t* LineSpec, int length, config_t c);
//void  AddComptonization(spec_t* ComSpec, spec_t* DetSpec, int length);
//void  ApplyPileUp(spec_t* DetSpec, spec_t* PUSpec, int length, config_t* c);
void WriteOutput(spec_t* DetSpec, int length, float LT);
//float myabs(float a);
