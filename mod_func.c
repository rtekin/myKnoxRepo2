#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _HH2_reg();
extern void _IM_reg();
extern void _IT_reg();
extern void _IT2_reg();
extern void _Ih_reg();
extern void _ampa_reg();
extern void _cadecay_reg();
extern void _gabaa_reg();
extern void _gabab_reg();
extern void _kleak_reg();
extern void _nmda_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," HH2.mod");
fprintf(stderr," IM.mod");
fprintf(stderr," IT.mod");
fprintf(stderr," IT2.mod");
fprintf(stderr," Ih.mod");
fprintf(stderr," ampa.mod");
fprintf(stderr," cadecay.mod");
fprintf(stderr," gabaa.mod");
fprintf(stderr," gabab.mod");
fprintf(stderr," kleak.mod");
fprintf(stderr," nmda.mod");
fprintf(stderr, "\n");
    }
_HH2_reg();
_IM_reg();
_IT_reg();
_IT2_reg();
_Ih_reg();
_ampa_reg();
_cadecay_reg();
_gabaa_reg();
_gabab_reg();
_kleak_reg();
_nmda_reg();
}
