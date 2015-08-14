#ifdef __INTEL_COMPILER

#define DOSAMPLING initsampler_mp_dosamplingfromc_

#elif defined __GNUC__

#define DOSAMPLING __initsampler_MOD_dosamplingfromc

#else

#error Do not know how to link to Fortran libraries, check symbol table for your platform (nm libchord.so | grep dosamplingfromc)

#endif

namespace chord
{

        extern "C" {
                void DOSAMPLING(double (*Lfunc)(int &nDims, double *theta, int &nDerived, double *phi, void *context), int &Ndim, int &nDerived, int &nLive, int &Nchords,  double *PriorsArray, char *Froot, void *context, double *output, int &do_grades, int &maxgrade, int *grades, int *grade_repeats, int *hypercube_indices, int *physical_indices);
        }

        static void Sample(double (*Lfunc)(int &nDims, double *theta, int &nDerived, double *phi, void *context), int &Ndim, int &nDerived, int &nLive, int &Nchords,  double *PriorsArray, char *Froot, void *context, double *output, int &do_grades, int &maxgrade, int *grades, int *grade_repeats, int *hypercube_indices, int *physical_indices){

                int i;
                for (i = strlen(Froot); i < 100; i++) Froot[i] = ' ';

                DOSAMPLING(Lfunc, Ndim, nDerived, nLive, Nchords, PriorsArray, Froot, context, output, do_grades, maxgrade, grades, grade_repeats, hypercube_indices, physical_indices);
        }

}

