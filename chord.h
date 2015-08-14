namespace chord
{

        extern "C" {
                void initsampler_mp_dosamplingfromc_(double (*Lfunc)(int &nDims, double *theta, int &nDerived, double *phi, void *context), int &Ndim, int &nDerived, int &nLive, int &Nchords,  double *PriorsArray, char *Froot, void *context);
        }

        static void Sample(double (*Lfunc)(int &nDims, double *theta, int &nDerived, double *phi, void *context), int &Ndim, int &nDerived, int &nLive, int &Nchords,  double *PriorsArray, char *Froot, void *context){

                int i;
                for (i = strlen(Froot); i < 100; i++) Froot[i] = ' ';

                initsampler_mp_dosamplingfromc_(Lfunc, Ndim, nDerived, nLive, Nchords, PriorsArray, Froot, context);
        }

}

