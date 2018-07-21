#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <cuda.h>
#include <errno.h>

#define NUCBITSZ 3
bool* getbool(char base){
  bool *n = new bool[3] {false, false, false};
  bool *a = new bool[3] {false, false, true };
  bool *t = new bool[3] {false, true, false };
  bool *g = new bool[3] {false, true, true  };
  bool *c = new bool[3] {true, false, false };
  bool *x = new bool[3] {true, true, true   };

  switch(base){
    case 'A': return a; break;
    case 'T': return t; break;
    case 'G': return g; break;
    case 'C': return c; break;
    case 'N': return n; break;
    case 'a': return a; break;
    case 't': return t; break;
    case 'g': return g; break;
    case 'c': return c; break;
    case 'n': return n; break;
    default : return x; break;
  }
}

__global__
void myKernel(bool* d_genome, bool* d_guides, bool* d_output, int GUIDESBATCHSZ){
  int startgenome = blockIdx.x * blockDim.x;
  int startwindow = startgenome + threadIdx.x;
   __shared__ bool blockGuide[20*3];
  int guidecnt = 0;
  while (guidecnt < (GUIDESBATCHSZ/3)/20 && startwindow < blockDim.x * gridDim.x ){
    for (int i = 0; i < 20*3; i++)
      blockGuide[i] = d_guides[20*3*guidecnt + i];
    bool resArray[20*3];
    for (int i = 0; i < 20*3; i++)
      resArray[i] = !(blockGuide[i] ^ d_genome[startwindow*3 + i]);
    bool res = true;
    for (int i = 0; i < 20*3; i++)
      res &= resArray[i];
    if (res) d_output[(gridDim.x * blockDim.x) * guidecnt + startwindow] = true;
    guidecnt += 1;
  }
}



#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true){
  if (code != cudaSuccess){
    fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

int main(int argc, char **argv){

  FILE *fgenome, *fguides, *fexact, *fmis1, *fmis2;
  int GRIDSZ = 32, BLOCKSZ = 1024, GUIDESZ = 1024, GUIDELENGTH = 20, TOTGUIDES = 0, GENOMEBATCHSZ, GUIDESBATCHSZ, OUTPUTBATCHSZ;
  bool *h_genome, *h_guides, *h_output, *d_genome, *d_output, *d_guides;
  int inp, igenome, ibingen, ngenome = 0, nguides = 0, iguides, ibingui, absidxgui, absidxgen = 0;
  bool flaggenome = true, flagguides = true;
  int *exact, *mis1, *mis2;
 
  GENOMEBATCHSZ = (BLOCKSZ*GRIDSZ+GUIDELENGTH - 1)*NUCBITSZ;
  GUIDESBATCHSZ = GUIDESZ*GUIDELENGTH*NUCBITSZ;
  OUTPUTBATCHSZ = GUIDESZ*GRIDSZ*BLOCKSZ;

  cudaFree(0);
  gpuErrchk(cudaMalloc((void**)&d_genome, sizeof(bool)*GENOMEBATCHSZ));
  gpuErrchk(cudaMalloc((void**)&d_guides, sizeof(bool)*GUIDESBATCHSZ));
  gpuErrchk(cudaMalloc((void**)&d_output, sizeof(bool)*OUTPUTBATCHSZ));
  h_genome = (bool*)calloc(GENOMEBATCHSZ, sizeof(bool));
  h_guides = (bool*)calloc(GUIDESBATCHSZ, sizeof(bool));
  h_output = (bool*)calloc(OUTPUTBATCHSZ, sizeof(bool));
  
  TOTGUIDES = 136632259;
  exact = (int*) calloc(TOTGUIDES, sizeof(int));

  fexact  = fopen("../Results/exact.txt", "w+");
  fmis1   = fopen("../Results/mis1.txt", "w+");
  fmis2   = fopen("../Results/mis2.txt", "w+");

  fgenome = fopen("chr1.txt", "r"); // GENOME
  while (flaggenome){
    printf("Genome: %d...\n", absidxgen);
    igenome = -1;
    ibingen = 0;
    while ((inp = fgetc(fgenome))!= EOF 
            && igenome < (GENOMEBATCHSZ/NUCBITSZ)){
      if ((char)inp != '\n'){
        bool *base;
        base = getbool((char)inp);
        for (int i = 0; i < NUCBITSZ; i++) h_genome[ibingen + i] = base[i];
        igenome++;
        ibingen += 3;
      }
    }
    if (inp == EOF) flaggenome = false;
    gpuErrchk(cudaMemcpy(d_genome, &(h_genome[ngenome]), sizeof(bool)*GENOMEBATCHSZ, cudaMemcpyHostToDevice));
    nguides = 0;
    flagguides = true;
    absidxgui = 0; 
    fguides = fopen("guides.txt", "r"); // GUIDES
    while (flagguides){
      printf("Guides: %d...\n", absidxgui);
      iguides = -1;
      ibingui = 0;
      while ((inp = fgetc(fguides)) != EOF && iguides < (GUIDESBATCHSZ/NUCBITSZ)){
        if ((char)inp != '\n'){
          bool * base;
          base = getbool((char)inp);
          for (int  i = 0; i < NUCBITSZ; i++) h_guides[ibingui + i] = base[i];
          iguides++;
          ibingui += 3;
        }
      }
      if (inp == EOF) flagguides = false;
      gpuErrchk(cudaMemcpy(d_guides, &(h_guides[nguides]), sizeof(bool)*GUIDESBATCHSZ, cudaMemcpyHostToDevice));
      myKernel<<<8, BLOCKSZ>>>(d_genome, d_guides, d_output, GUIDESBATCHSZ);
      gpuErrchk(cudaMemcpy(h_output, d_output, sizeof(bool)*OUTPUTBATCHSZ, cudaMemcpyDeviceToHost));
      for (int i = 0; i < iguides/GUIDELENGTH; i++) 
        for (int j = 0; j < (igenome - GUIDELENGTH + 1); j++){
          //if (h_output[i*iguides+j]) fprintf(fexact, "%d,%d:", absidxgui + i, absidxgen + j);
          if (h_output[i*iguides+j]) exact[i+absidxgui] += 1;
        }
      nguides += ibingui;
      absidxgui += iguides;
    }
    /*for(int i = 0; i < GUIDESZ; i++)
      fprintf(fexact, "%d:%d,", absidxgui+i, exact[i]);*/
    fclose(fguides);
    ngenome += ibingen;
    absidxgen += igenome;
  }
  fclose(fgenome);
  printf("Number of GUIDES: %d", TOTGUIDES);
  for (int i = 0; i < TOTGUIDES; i++)
    if (exact[i] != 0) fprintf(fexact, "%d:%d,", i, exact[i]);
  fclose(fexact);
  fclose(fmis1);  
  fclose(fmis2);

  free(h_guides);
  free(h_genome);
  free(h_output);
  free(exact);
  cudaFree(d_genome);
  cudaFree(d_guides);
  cudaFree(d_output);
  return 0;
}
