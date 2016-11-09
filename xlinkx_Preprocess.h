/*
   Copyright 2012 University of Washington

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#ifndef _XLINKXPREPROCESS_H_
#define _XLINKXPREPROCESS_H_

class xlinkx_preprocess
{
public:
   xlinkx_preprocess();
   ~xlinkx_preprocess();

   static void Reset();
   static void LoadAndPreprocessSpectra(MSReader &mstReader,
                                        int iScanNumber,
                                        double dMZ1,
                                        double dMZ2);
   static bool DoneProcessingAllSpectra();
   static bool AllocateMemory(int maxNumThreads);
   static bool DeallocateMemory(int maxNumThreads);

private:

   // Private static methods
   static bool PreprocessSpectrum(Spectrum &spec,
                                  double dMZ1,
                                  double dMZ2,
                                  double *pdTmpRawData,
                                  double *pdTmpFastXcorrData,
                                  double *pdTmpCorrelationData);
   static bool CheckExistOutFile(int iCharge,
                                 int iScanNum);
   static bool AdjustMassTol(struct Query *pScoring);
   static bool CheckActivationMethodFilter(MSActivation act);
   static bool CheckExit(int iAnalysisType,
                         int iScanNum,
                         int iTotalScans,
                         int iLastScan,
                         int iReaderLastScan,
                         int iNumSpectraLoaded);
   static bool Preprocess(struct Query *pScoring,
                          Spectrum mstSpectrum,
                          double dMZ1,
                          double dMZ2,
                          double *pdTmpRawData,
                          double *pdTmpFastXcorrData,
                          double *pdTmpCorrelationData);
   static bool LoadIons(struct Query *pScoring,
                        double *pdTmpRawData,
                        Spectrum mstSpectrum,
                        struct PreprocessStruct *pPre,
                        double dMZ1,
                        double dMZ2);
   static void MakeCorrData(double *pdTmpRawData,
                            double *pdTmpCorrelationData,
                            struct Query *pScoring,
                            struct PreprocessStruct *pPre);
   static bool IsValidInputType(int inputType);

   // Private member variables
   static bool _bFirstScan;
   static bool _bDoneProcessingAllSpectra;

   //MH: Common memory to be shared by all threads during spectral processing
   static bool *pbMemoryPool;                 //MH: Regulator of memory use
   static double **ppdTmpRawDataArr;          //MH: Number of arrays equals threads
   static double **ppdTmpFastXcorrDataArr;    //MH: Ditto
   static double **ppdTmpCorrelationDataArr;  //MH: Ditto
};

#endif // _XLINKXPREPROCESS_H_
