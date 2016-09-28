/*
   Copyright 2016 University of Washington

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

#include "Common.h"
#include "xlinkx_Preprocess.h"
#include "xlinkx_DataInternal.h"

std::vector<Query*>           g_pvQuery;
std::vector<InputFileInfo *>  g_pvInputFiles;
StaticParams                  g_staticParams;
MassRange                     g_massRange;


bool xlinkx_preprocess::_bDoneProcessingAllSpectra;
bool xlinkx_preprocess::_bFirstScan;
bool *xlinkx_preprocess::pbMemoryPool;
double **xlinkx_preprocess::ppdTmpRawDataArr;
double **xlinkx_preprocess::ppdTmpFastXcorrDataArr;
double **xlinkx_preprocess::ppdTmpCorrelationDataArr;
double **xlinkx_preprocess::ppdTmpSmoothedSpectrumArr;
double **xlinkx_preprocess::ppdTmpPeakExtractedArr;

// Generate data for both sp scoring (pfSpScoreData) and xcorr analysis (FastXcorr).
xlinkx_preprocess::xlinkx_preprocess()
{
}


xlinkx_preprocess::~xlinkx_preprocess()
{
}

void xlinkx_preprocess::Reset()
{
    _bFirstScan = true;
    _bDoneProcessingAllSpectra = false;
}

void xlinkx_preprocess::LoadAndPreprocessSpectra(MSReader &mstReader,
                                                 int iScanNumber)
{
   Spectrum mstSpectrum;           // For holding spectrum.

   g_massRange.iMaxFragmentCharge = 0;

   // Loads in MSMS spectrum data.
   mstReader.readFile(NULL, mstSpectrum, iScanNumber);

   if (mstSpectrum.getScanNumber() != 0)   // should not be needed by quick sanity check to make sure scan is read
   {
      int iNumClearedPeaks = 0;

      // Clear out m/z range if clear_mz_range parameter is specified
      // Accomplish this by setting corresponding intensity to 0
      if (g_staticParams.options.clearMzRange.dEnd > 0.0
            && g_staticParams.options.clearMzRange.dStart <= g_staticParams.options.clearMzRange.dEnd)
      {
         int i=0;

         while (true)
         {
            if (i >= mstSpectrum.size() || mstSpectrum.at(i).mz > g_staticParams.options.clearMzRange.dEnd)
               break;

            if (mstSpectrum.at(i).mz >= g_staticParams.options.clearMzRange.dStart
                  && mstSpectrum.at(i).mz <= g_staticParams.options.clearMzRange.dEnd)
            {
               mstSpectrum.at(i).intensity = 0.0;
               iNumClearedPeaks++;
            }

            i++;
         }
      }

      if (mstSpectrum.size()-iNumClearedPeaks >= g_staticParams.options.iMinPeaks)
      {
         if (CheckActivationMethodFilter(mstSpectrum.getActivationMethod()))
         {
            int i=0;
            PreprocessSpectrum(mstSpectrum,
                  ppdTmpRawDataArr[i],
                  ppdTmpFastXcorrDataArr[i],
                  ppdTmpCorrelationDataArr[i],
                  ppdTmpSmoothedSpectrumArr[i],
                  ppdTmpPeakExtractedArr[i]);

         }
      }
   }
}


bool xlinkx_preprocess::DoneProcessingAllSpectra()
{
   return _bDoneProcessingAllSpectra;
}


bool xlinkx_preprocess::Preprocess(struct Query *pScoring,
                                 Spectrum mstSpectrum,
                                 double *pdTmpRawData,
                                 double *pdTmpFastXcorrData,
                                 double *pdTmpCorrelationData,
                                 double *pdTmpSmoothedSpectrum,
                                 double *pdTmpPeakExtracted)
{
   int i;
   int x;
   int y;
   struct msdata pTmpSpData[NUM_SP_IONS];
   struct PreprocessStruct pPre;

   pPre.iHighestIon = 0;
   pPre.dHighestIntensity = 0;

   //MH: Find appropriately sized array cushion based on user parameters. Fixes error found by Patrick Pedrioli for
   // very wide mass tolerance searches (i.e. 500 Da).
   double dCushion = 0.0;
   if (g_staticParams.tolerances.iMassToleranceUnits == 0) // amu
   {
      dCushion = g_staticParams.tolerances.dInputTolerance;

      if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
      {
        dCushion *= 8; //MH: hope +8 is large enough charge because g_staticParams.options.iEndCharge can be overridden.
      }
   }
   else if (g_staticParams.tolerances.iMassToleranceUnits == 1) // mmu
   {
      dCushion = g_staticParams.tolerances.dInputTolerance * 0.001;

      if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
      {
         dCushion *= 8; //MH: hope +8 is large enough charge because g_staticParams.options.iEndCharge can be overridden.
      }
   }
   else // ppm
   {
      dCushion = g_staticParams.tolerances.dInputTolerance * g_staticParams.options.dPeptideMassHigh / 1000000.0;
   }

   // initialize these temporary arrays before re-using
   size_t iTmp = (size_t)((g_staticParams.options.dPeptideMassHigh + dCushion + 2.0) * g_staticParams.dInverseBinWidth)*sizeof(double);

   memset(pdTmpRawData, 0, iTmp);
   memset(pdTmpFastXcorrData, 0, iTmp);
   memset(pdTmpCorrelationData, 0, iTmp);
   memset(pdTmpSmoothedSpectrum, 0, iTmp);
   memset(pdTmpPeakExtracted, 0, iTmp);

   // pdTmpRawData is a binned array holding raw data
   if (!LoadIons(pScoring, pdTmpRawData, mstSpectrum, &pPre))
   {
      return false;
   }

   try
   {
      pScoring->pfFastXcorrData = new float[pScoring->_spectrumInfoInternal.iArraySize]();
   }
   catch (std::bad_alloc& ba)
   {
      fprintf(stderr,  " Error - new(pfFastXcorrData[%d]). bad_alloc: %s.\n", pScoring->_spectrumInfoInternal.iArraySize, ba.what());
      fprintf(stderr, "Xlinkx ran out of memory. Look into \"spectrum_batch_size\"\n");
      fprintf(stderr, "parameters to address mitigate memory use.\n");
      return false;
   }

   if (g_staticParams.ionInformation.bUseNeutralLoss
         && (g_staticParams.ionInformation.iIonVal[ION_SERIES_A]
            || g_staticParams.ionInformation.iIonVal[ION_SERIES_B]
            || g_staticParams.ionInformation.iIonVal[ION_SERIES_Y]))
   {
      try
      {
         pScoring->pfFastXcorrDataNL = new float[pScoring->_spectrumInfoInternal.iArraySize]();
      }
      catch (std::bad_alloc& ba)
      {
         fprintf(stderr,  " Error - new(pfFastXcorrDataNL[%d]). bad_alloc: %s.\n", pScoring->_spectrumInfoInternal.iArraySize, ba.what());
         fprintf(stderr, "Xlinkx ran out of memory. Look into \"spectrum_batch_size\"\n");
         fprintf(stderr, "parameters to address mitigate memory use.\n");
         return false;
      }
   }

   // Create data for correlation analysis.
   // pdTmpRawData intensities are normalized to 100; pdTmpCorrelationData is windowed
   MakeCorrData(pdTmpRawData, pdTmpCorrelationData, pScoring, &pPre);

   // Make fast xcorr spectrum.
   double dSum=0.0;
   int iTmpRange = 2*g_staticParams.iXcorrProcessingOffset + 1;
   double dTmp = 1.0 / (double)(iTmpRange - 1);

   dSum=0.0;
   for (i=0; i<g_staticParams.iXcorrProcessingOffset; i++)
      dSum += pdTmpCorrelationData[i];
   for (i=g_staticParams.iXcorrProcessingOffset; i < pScoring->_spectrumInfoInternal.iArraySize + g_staticParams.iXcorrProcessingOffset; i++)
   {
      if (i<pScoring->_spectrumInfoInternal.iArraySize)
         dSum += pdTmpCorrelationData[i];
      if (i>=iTmpRange)
         dSum -= pdTmpCorrelationData[i-iTmpRange];
      pdTmpFastXcorrData[i-g_staticParams.iXcorrProcessingOffset] = (dSum - pdTmpCorrelationData[i-g_staticParams.iXcorrProcessingOffset])* dTmp;
   }

   pScoring->pfFastXcorrData[0] = 0.0;
   for (i=1; i<pScoring->_spectrumInfoInternal.iArraySize; i++)
   {
      double dTmp = pdTmpCorrelationData[i] - pdTmpFastXcorrData[i];

      pScoring->pfFastXcorrData[i] = (float)dTmp;

      // Add flanking peaks if used
      if (g_staticParams.ionInformation.iTheoreticalFragmentIons == 0)
      {
         int iTmp;

         iTmp = i-1;
         pScoring->pfFastXcorrData[i] += (float) ((pdTmpCorrelationData[iTmp] - pdTmpFastXcorrData[iTmp])*0.5);

         iTmp = i+1;
         if (iTmp < pScoring->_spectrumInfoInternal.iArraySize)
            pScoring->pfFastXcorrData[i] += (float) ((pdTmpCorrelationData[iTmp] - pdTmpFastXcorrData[iTmp])*0.5);
      }

      // If A, B or Y ions and their neutral loss selected, roll in -17/-18 contributions to pfFastXcorrDataNL
      if (g_staticParams.ionInformation.bUseNeutralLoss
            && (g_staticParams.ionInformation.iIonVal[ION_SERIES_A]
               || g_staticParams.ionInformation.iIonVal[ION_SERIES_B]
               || g_staticParams.ionInformation.iIonVal[ION_SERIES_Y]))
      {
         int iTmp;

         pScoring->pfFastXcorrDataNL[i] = pScoring->pfFastXcorrData[i];

         iTmp = i-g_staticParams.precalcMasses.iMinus17;
         if (iTmp>= 0)
         {
            pScoring->pfFastXcorrDataNL[i] += (float)((pdTmpCorrelationData[iTmp] - pdTmpFastXcorrData[iTmp]) * 0.2);
         }

         iTmp = i-g_staticParams.precalcMasses.iMinus18;
         if (iTmp>= 0)
         {
            pScoring->pfFastXcorrDataNL[i] += (float)((pdTmpCorrelationData[iTmp] - pdTmpFastXcorrData[iTmp]) * 0.2);
         }

      }
   }

   // Using sparse matrix which means we free pScoring->pfFastXcorrData, ->pfFastXcorrDataNL here
   // If A, B or Y ions and their neutral loss selected, roll in -17/-18 contributions to pfFastXcorrDataNL.
   if (g_staticParams.ionInformation.bUseNeutralLoss
         && (g_staticParams.ionInformation.iIonVal[ION_SERIES_A]
            || g_staticParams.ionInformation.iIonVal[ION_SERIES_B]
            || g_staticParams.ionInformation.iIonVal[ION_SERIES_Y]))
   {
      pScoring->iFastXcorrDataNL=pScoring->_spectrumInfoInternal.iArraySize/SPARSE_MATRIX_SIZE+1;

      try
      {
         pScoring->ppfSparseFastXcorrDataNL = new float*[pScoring->iFastXcorrDataNL]();
      }
      catch (std::bad_alloc& ba)
      {
         fprintf(stderr,  " Error - new(pScoring->ppfSparseFastXcorrDataNL[%d]). bad_alloc: %s.", pScoring->iFastXcorrDataNL, ba.what());
         fprintf(stderr, "Xlinkx ran out of memory. Look into \"spectrum_batch_size\"\n");
         fprintf(stderr, "parameters to address mitigate memory use.\n");
         return false;
      }

      for (i=1; i<pScoring->_spectrumInfoInternal.iArraySize; i++)
      {
         if (pScoring->pfFastXcorrDataNL[i]>FLOAT_ZERO || pScoring->pfFastXcorrDataNL[i]<-FLOAT_ZERO)
         {
            x=i/SPARSE_MATRIX_SIZE;
            if (pScoring->ppfSparseFastXcorrDataNL[x]==NULL)
            {
               try
               {
                  pScoring->ppfSparseFastXcorrDataNL[x] = new float[SPARSE_MATRIX_SIZE]();
               }
               catch (std::bad_alloc& ba)
               {
                  fprintf(stderr,  " Error - new(pScoring->ppfSparseFastXcorrDataNL[%d][%d]). bad_alloc: %s.\n", x, SPARSE_MATRIX_SIZE, ba.what());
                  fprintf(stderr, "Xlinkx ran out of memory. Look into \"spectrum_batch_size\"\n");
                  fprintf(stderr, "parameters to address mitigate memory use.\n");
                  return false;
               }
               for (y=0; y<SPARSE_MATRIX_SIZE; y++)
                  pScoring->ppfSparseFastXcorrDataNL[x][y]=0;
            }
            y=i-(x*SPARSE_MATRIX_SIZE);
            pScoring->ppfSparseFastXcorrDataNL[x][y] = pScoring->pfFastXcorrDataNL[i];
         }
      }

      delete[] pScoring->pfFastXcorrDataNL;
      pScoring->pfFastXcorrDataNL = NULL;

   }

   pScoring->iFastXcorrData=pScoring->_spectrumInfoInternal.iArraySize/SPARSE_MATRIX_SIZE+1;

   //MH: Fill sparse matrix
   try
   {
      pScoring->ppfSparseFastXcorrData = new float*[pScoring->iFastXcorrData]();
   }
   catch (std::bad_alloc& ba)
   {
      fprintf(stderr,  " Error - new(pScoring->ppfSparseFastXcorrData[%d]). bad_alloc: %s.\n", pScoring->iFastXcorrData, ba.what());
      fprintf(stderr, "Xlinkx ran out of memory. Look into \"spectrum_batch_size\"\n");
      fprintf(stderr, "parameters to address mitigate memory use.\n");
      return false;
   }

   for (i=1; i<pScoring->_spectrumInfoInternal.iArraySize; i++)
   {
      if (pScoring->pfFastXcorrData[i]>FLOAT_ZERO || pScoring->pfFastXcorrData[i]<-FLOAT_ZERO)
      {
         x=i/SPARSE_MATRIX_SIZE;
         if (pScoring->ppfSparseFastXcorrData[x]==NULL)
         {
            try
            {
               pScoring->ppfSparseFastXcorrData[x] = new float[SPARSE_MATRIX_SIZE]();
            }
            catch (std::bad_alloc& ba)
            {
               fprintf(stderr,  " Error - new(pScoring->ppfSparseFastXcorrData[%d][%d]). bad_alloc: %s.\n", x, SPARSE_MATRIX_SIZE, ba.what());
               fprintf(stderr, "Xlinkx ran out of memory. Look into \"spectrum_batch_size\"\n");
               fprintf(stderr, "parameters to address mitigate memory use.\n");
               return false;
            }
            for (y=0; y<SPARSE_MATRIX_SIZE; y++)
               pScoring->ppfSparseFastXcorrData[x][y]=0;
         }
         y=i-(x*SPARSE_MATRIX_SIZE);
         pScoring->ppfSparseFastXcorrData[x][y] = pScoring->pfFastXcorrData[i];
      }
   }

   delete[] pScoring->pfFastXcorrData;
   pScoring->pfFastXcorrData = NULL;

   // Create data for sp scoring.
   // Arbitrary bin size cutoff to do smoothing, peak extraction.
   if (g_staticParams.tolerances.dFragmentBinSize >= 0.10)
   {
      if (!Smooth(pdTmpRawData, pScoring->_spectrumInfoInternal.iArraySize, pdTmpSmoothedSpectrum))
      {
         return false;
      }

      if (!PeakExtract(pdTmpRawData, pScoring->_spectrumInfoInternal.iArraySize, pdTmpPeakExtracted))
      {
         return false;
      }
   }

   for (i=0; i<NUM_SP_IONS; i++)
   {
      pTmpSpData[i].dIon = 0.0;
      pTmpSpData[i].dIntensity = 0.0;
   }

   GetTopIons(pdTmpRawData, &(pTmpSpData[0]), pScoring->_spectrumInfoInternal.iArraySize);

   qsort(pTmpSpData, NUM_SP_IONS, sizeof(struct msdata), QsortByIon);

   // Modify for Sp data.
   StairStep(pTmpSpData);

   try
   {
      pScoring->pfSpScoreData = new float[pScoring->_spectrumInfoInternal.iArraySize]();
   }
   catch (std::bad_alloc& ba)
   {
      fprintf(stderr,  " Error - new(pfSpScoreData[%d]). bad_alloc: %s.\n", pScoring->_spectrumInfoInternal.iArraySize, ba.what());
      fprintf(stderr, "Xlinkx ran out of memory. Look into \"spectrum_batch_size\"\n");
      fprintf(stderr, "parameters to address mitigate memory use.\n");
      return false;
   }

   // note that pTmpSpData[].dIon values are already BIN'd
   for (i=0; i<NUM_SP_IONS; i++)
      pScoring->pfSpScoreData[(int)(pTmpSpData[i].dIon)] = (float) pTmpSpData[i].dIntensity;

   // MH: Fill sparse matrix for SpScore
   pScoring->iSpScoreData = pScoring->_spectrumInfoInternal.iArraySize / SPARSE_MATRIX_SIZE + 1;

   try
   {
      pScoring->ppfSparseSpScoreData = new float*[pScoring->iSpScoreData]();
   }
   catch (std::bad_alloc& ba)
   {
      fprintf(stderr,  " Error - new(pScoring->ppfSparseSpScoreData[%d]). bad_alloc: %s.\n", pScoring->iSpScoreData, ba.what());
      fprintf(stderr, "Xlinkx ran out of memory. Look into \"spectrum_batch_size\"\n");
      fprintf(stderr, "parameters to address mitigate memory use.\n");
      return false;
   }

   for (i=0; i<pScoring->_spectrumInfoInternal.iArraySize; i++)
   {
      if (pScoring->pfSpScoreData[i] > FLOAT_ZERO)
      {
         x=i/SPARSE_MATRIX_SIZE;
         if (pScoring->ppfSparseSpScoreData[x]==NULL)
         {
            try
            {
               pScoring->ppfSparseSpScoreData[x] = new float[SPARSE_MATRIX_SIZE]();
            }
            catch (std::bad_alloc& ba)
            {
               fprintf(stderr,  " Error - new(pScoring->ppfSparseSpScoreData[%d][%d]). bad_alloc: %s.\n", x, SPARSE_MATRIX_SIZE, ba.what());
               fprintf(stderr, "Xlinkx ran out of memory. Look into \"spectrum_batch_size\"\n");
               fprintf(stderr, "parameters to address mitigate memory use.\n");
               return false;
            }
            for (y=0; y<SPARSE_MATRIX_SIZE; y++)
               pScoring->ppfSparseSpScoreData[x][y]=0;
         }
         y=i-(x*SPARSE_MATRIX_SIZE);
         pScoring->ppfSparseSpScoreData[x][y] = pScoring->pfSpScoreData[i];
      }
   }

   delete[] pScoring->pfSpScoreData;
   pScoring->pfSpScoreData = NULL;

   return true;
}


bool xlinkx_preprocess::CheckActivationMethodFilter(MSActivation act)
{
   bool bSearchSpectrum = true;

   // Check possible activation method filter.
   if (strcmp(g_staticParams.options.szActivationMethod, "ALL")!=0 && (act != mstNA))
   {
      if (!strcmp(g_staticParams.options.szActivationMethod, "CID") && (act != mstCID))
      {
         bSearchSpectrum = 0;
      }
      else if (!strcmp(g_staticParams.options.szActivationMethod, "HCD") && (act != mstHCD))
      {
         bSearchSpectrum = 0;
      }
      else if (!strcmp(g_staticParams.options.szActivationMethod, "ETD") && (act != mstETD))
      {
         bSearchSpectrum = 0;
      }
      else if (!strcmp(g_staticParams.options.szActivationMethod, "ECD") && (act != mstECD))
      {
         bSearchSpectrum = 0;
      }
      else if (!strcmp(g_staticParams.options.szActivationMethod, "PQD") && (act != mstPQD))
      {
         bSearchSpectrum = 0;
      }
      else if (!strcmp(g_staticParams.options.szActivationMethod, "IRMPD") && (act != mstIRMPD))
      {
         bSearchSpectrum = 0;
      }
   }

   return bSearchSpectrum;
}

bool xlinkx_preprocess::CheckExit(int iAnalysisType,
                                int iScanNum,
                                int iTotalScans,
                                int iLastScan,
                                int iReaderLastScan,
                                int iNumSpectraLoaded)
{
   if (iAnalysisType == AnalysisType_SpecificScan)
   {
      _bDoneProcessingAllSpectra = true;
      return true;
   }

   if (iAnalysisType == AnalysisType_SpecificScanRange)
   {
      if (iLastScan > 0)
      {
         if (iScanNum >= iLastScan)
         {
            _bDoneProcessingAllSpectra = true;
            return true;
         }
      }
   }

   if (iAnalysisType == AnalysisType_EntireFile
         && IsValidInputType(g_staticParams.inputFile.iInputType)
         && iScanNum == 0)
   {
      _bDoneProcessingAllSpectra = true;
      return true;
   }

   // Horrible way to exit as this typically requires a quick cycle through
   // while loop but not sure what else to do when getScanNumber() returns 0
   // for non MS/MS scans.
   if (IsValidInputType(g_staticParams.inputFile.iInputType)
         && iTotalScans > iReaderLastScan)
   {
      _bDoneProcessingAllSpectra = true;
      return true;
   }

   if ((g_staticParams.options.iSpectrumBatchSize != 0)
         && (iNumSpectraLoaded >= g_staticParams.options.iSpectrumBatchSize))
   {
      return true;
   }

   return false;
}


bool xlinkx_preprocess::PreprocessSpectrum(Spectrum &spec,
                                         double *pdTmpRawData,
                                         double *pdTmpFastXcorrData,
                                         double *pdTmpCorrelationData,
                                         double *pdTmpSmoothedSpectrum,
                                         double *pdTmpPeakExtracted)
{
   int z;
   int zStop;

   int iScanNumber = spec.getScanNumber();

   int iAddCharge = -1;  // specifies how charge states are going to be applied
                         //-1 = should never apply
                         // 0 = use charge in file, else use range
                         // 1 = use precursor_charge range
                         // 2 = search only charge state in file within precursor_charge
                         // 3 = use charge in file else use 1+ or precursor_charge range
   int bFileHasCharge = 0;
   if (spec.sizeZ() > 0)
      bFileHasCharge = 1;

   if (g_staticParams.options.iStartCharge > 0)
   {
      if (!bFileHasCharge)    // override_charge specified, no charge in file
      {
         if (g_staticParams.options.bOverrideCharge == 0
               || g_staticParams.options.bOverrideCharge == 1
               || g_staticParams.options.bOverrideCharge == 2)
         {
            iAddCharge = 2;
         }
         else if (g_staticParams.options.bOverrideCharge == 3)
         {
            iAddCharge = 3;
         }
      }
      else                    // have a charge from file //
      {
         // bOverrideCharge == 0, 2, 3 are not relevant here
         if (g_staticParams.options.bOverrideCharge == 1)
         {
            iAddCharge = 2;
         }
         else
         {
            iAddCharge = 0;   // do nothing
         }
      }
   }
   else  // precursor_charge range not specified
   {
      if (!bFileHasCharge)    // no charge in file
      {
         iAddCharge = 1;
      }
      else
      {
         iAddCharge = 0;      // have a charge from file; nothing to do
      }
   }

   if (iAddCharge == 2)  // use specific charge range
   {
      // if charge is already specified, don't re-add it here when overriding charge range
      int iChargeFromFile = 0;
      if (bFileHasCharge)                  // FIX: no reason that file only has one charge so iChargeFromFile should be an array
         iChargeFromFile= spec.atZ(0).z;   // should read all charge states up to spec.sizeZ();

      for (z=g_staticParams.options.iStartCharge; z<=g_staticParams.options.iEndCharge; z++)
      {
         if (z != iChargeFromFile)
            spec.addZState(z, spec.getMZ() * z - (z-1) * PROTON_MASS);
      }
   }
   else if (iAddCharge == 1 || iAddCharge == 3)  // do 1+ or charge range rule
   {
      int i=0;
      double dSumBelow = 0.0;
      double dSumTotal = 0.0;

      while (true)
      {
         if (i >= spec.size())
            break;

         dSumTotal += spec.at(i).intensity;

         if (spec.at(i).mz < spec.getMZ())
            dSumBelow += spec.at(i).intensity;

         i++;
      }

      if (isEqual(dSumTotal, 0.0) || ((dSumBelow/dSumTotal) > 0.95))
      {
         z = 1;
         spec.addZState(z, spec.getMZ() * z - (z-1) * PROTON_MASS);
      }
      else
      {
         if (iAddCharge == 1)  // 2+ and 3+
         {
            z=2;
            spec.addZState(z, spec.getMZ() * z - (z-1) * PROTON_MASS);
            z=3;
            spec.addZState(z, spec.getMZ() * z - (z-1) * PROTON_MASS);
         }
         else // iAddCharge == 3
         {
            // This option will either use charge from file or
            // charge_range with 1+ rule.  So no redundant addition
            // of charges possible

            for (z=g_staticParams.options.iStartCharge; z<=g_staticParams.options.iEndCharge; z++)
            {
               spec.addZState(z, spec.getMZ() * z - (z-1) * PROTON_MASS);
            }
         }
      }
   }
   else if (iAddCharge == -1)  // should never get here
   {
      fprintf(stderr,  " Error - iAddCharge=%d\n",  iAddCharge);
      return false;
   }

   // Set our boundaries for multiple z lines.
   zStop = spec.sizeZ();

   for (z=0; z<zStop; z++)
   {
      if (g_staticParams.options.bOverrideCharge == 2 && g_staticParams.options.iStartCharge > 0)
      {
         // ignore spectra that aren't 2+ or 3+.
         if (spec.atZ(z).z < g_staticParams.options.iStartCharge || z > g_staticParams.options.iEndCharge)
         {
            continue;
         }
      }

      int iPrecursorCharge = spec.atZ(z).z;  // I need this before iChargeState gets assigned.
      double dMass = spec.atZ(z).mh;

      if (!g_staticParams.options.bOverrideCharge
            || g_staticParams.options.iStartCharge == 0
            || (g_staticParams.options.bOverrideCharge == 3 && bFileHasCharge)
            || (g_staticParams.options.bOverrideCharge
               && (g_staticParams.options.iStartCharge > 0
                  && ((iPrecursorCharge>=g_staticParams.options.iStartCharge
                        && iPrecursorCharge<=g_staticParams.options.iEndCharge )
                     || (iPrecursorCharge == 1 && g_staticParams.options.bOverrideCharge == 3)))))
      {
         if (CheckExistOutFile(iPrecursorCharge, iScanNumber)
               && (isEqual(g_staticParams.options.dPeptideMassLow, 0.0)
                  || ((dMass >= g_staticParams.options.dPeptideMassLow)
                     && (dMass <= g_staticParams.options.dPeptideMassHigh)))
               && iPrecursorCharge <= g_staticParams.options.iMaxPrecursorCharge)
         {
            Query *pScoring = new Query();

            pScoring->_pepMassInfo.dExpPepMass = dMass;
            pScoring->_spectrumInfoInternal.iChargeState = iPrecursorCharge;
            pScoring->_spectrumInfoInternal.dTotalIntensity = 0.0;
            pScoring->_spectrumInfoInternal.dRTime = 60.0*spec.getRTime();;
            pScoring->_spectrumInfoInternal.iScanNumber = iScanNumber;

            if (iPrecursorCharge == 1)
               pScoring->_spectrumInfoInternal.iMaxFragCharge = 1;
            else
            {
               pScoring->_spectrumInfoInternal.iMaxFragCharge = iPrecursorCharge - 1;

               if (pScoring->_spectrumInfoInternal.iMaxFragCharge > g_staticParams.options.iMaxFragmentCharge)
                  pScoring->_spectrumInfoInternal.iMaxFragCharge = g_staticParams.options.iMaxFragmentCharge;
            }

            //MH: Find appropriately sized array cushion based on user parameters. Fixes error found by Patrick Pedrioli for
            // very wide mass tolerance searches (i.e. 500 Da).
            double dCushion = 0.0;
            if (g_staticParams.tolerances.iMassToleranceUnits == 0) // amu
            {
               dCushion = g_staticParams.tolerances.dInputTolerance;

               if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
               {
                  dCushion *= pScoring->_spectrumInfoInternal.iChargeState;
               }
            }
            else if (g_staticParams.tolerances.iMassToleranceUnits == 1) // mmu
            {
               dCushion = g_staticParams.tolerances.dInputTolerance * 0.001;

               if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
               {
                  dCushion *= pScoring->_spectrumInfoInternal.iChargeState;
               }
            }
            else // ppm
            {
               dCushion = g_staticParams.tolerances.dInputTolerance * dMass / 1000000.0;
            }
            pScoring->_spectrumInfoInternal.iArraySize = (int)((dMass + dCushion + 2.0) * g_staticParams.dInverseBinWidth);

            // g_massRange.iMaxFragmentCharge is global maximum fragment ion charge across all spectra.
            if (pScoring->_spectrumInfoInternal.iMaxFragCharge > g_massRange.iMaxFragmentCharge)
            {
               g_massRange.iMaxFragmentCharge = pScoring->_spectrumInfoInternal.iMaxFragCharge;
            }

            if (!AdjustMassTol(pScoring))
            {
               return false;
            }

            // Populate pdCorrelation data.
            // NOTE: there must be a good way of doing this just once per spectrum instead
            //       of repeating for each charge state.
            if (!Preprocess(pScoring, spec, pdTmpRawData, pdTmpFastXcorrData,
                     pdTmpCorrelationData, pdTmpSmoothedSpectrum, pdTmpPeakExtracted))
            {
               return false;
            }

            g_pvQuery.push_back(pScoring);
         }
      }
   }

   return true;
}


// Skip repeating a search if output exists only works for .out files
bool xlinkx_preprocess::CheckExistOutFile(int iCharge,
                                        int iScanNum)
{
   bool bSearchSpectrum = 1;

   if (g_staticParams.options.bOutputOutFiles
         && g_staticParams.options.bSkipAlreadyDone
         && !g_staticParams.options.bOutputPepXMLFile
         && !g_staticParams.options.bOutputPercolatorFile)
   {
      char szOutputFileName[SIZE_FILE];
      char *pStr;
      FILE *fpcheck;

      if ( (pStr = strrchr(g_staticParams.inputFile.szBaseName, '\\')) == NULL
            && (pStr = strrchr(g_staticParams.inputFile.szBaseName, '/')) == NULL)
      {
         pStr = g_staticParams.inputFile.szBaseName;
      }
      else
         (*pStr)++;

      sprintf(szOutputFileName, "%s/%s.%.5d.%.5d.%d.out",
            g_staticParams.inputFile.szBaseName,
            pStr,
            iScanNum,
            iScanNum,
            iCharge);

      // Check existence of .out file.
      if ((fpcheck = fopen(szOutputFileName, "r")) != NULL)
      {
         bSearchSpectrum = 0;
         fclose(fpcheck);
      }
   }

   return bSearchSpectrum;
}


bool xlinkx_preprocess::AdjustMassTol(struct Query *pScoring)
{
   if (g_staticParams.tolerances.iMassToleranceUnits == 0) // amu
   {
      pScoring->_pepMassInfo.dPeptideMassTolerance = g_staticParams.tolerances.dInputTolerance;

      if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
      {
         pScoring->_pepMassInfo.dPeptideMassTolerance *= pScoring->_spectrumInfoInternal.iChargeState;
      }
   }
   else if (g_staticParams.tolerances.iMassToleranceUnits == 1) // mmu
   {
      pScoring->_pepMassInfo.dPeptideMassTolerance = g_staticParams.tolerances.dInputTolerance * 0.001;

      if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
      {
         pScoring->_pepMassInfo.dPeptideMassTolerance *= pScoring->_spectrumInfoInternal.iChargeState;
      }
   }
   else // ppm
   {
      pScoring->_pepMassInfo.dPeptideMassTolerance = g_staticParams.tolerances.dInputTolerance
         * pScoring->_pepMassInfo.dExpPepMass / 1000000.0;
   }

   if (g_staticParams.tolerances.iIsotopeError == 0)
   {
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus = pScoring->_pepMassInfo.dExpPepMass
         - pScoring->_pepMassInfo.dPeptideMassTolerance;

      pScoring->_pepMassInfo.dPeptideMassTolerancePlus = pScoring->_pepMassInfo.dExpPepMass
         + pScoring->_pepMassInfo.dPeptideMassTolerance;
   }
   else if (g_staticParams.tolerances.iIsotopeError == 1) // search 0, +1 isotope windows
   {
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus = pScoring->_pepMassInfo.dExpPepMass
         - pScoring->_pepMassInfo.dPeptideMassTolerance - C13_DIFF * PROTON_MASS;

      pScoring->_pepMassInfo.dPeptideMassTolerancePlus = pScoring->_pepMassInfo.dExpPepMass
         + pScoring->_pepMassInfo.dPeptideMassTolerance;
   }
   else if (g_staticParams.tolerances.iIsotopeError == 2) // search 0, +1, +2 isotope windows
   {
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus = pScoring->_pepMassInfo.dExpPepMass
         - pScoring->_pepMassInfo.dPeptideMassTolerance - 2.0 * C13_DIFF * PROTON_MASS;

      pScoring->_pepMassInfo.dPeptideMassTolerancePlus = pScoring->_pepMassInfo.dExpPepMass
         + pScoring->_pepMassInfo.dPeptideMassTolerance;
   }
   else if (g_staticParams.tolerances.iIsotopeError == 3) // search 0, +1, +2, +3 isotope windows
   {
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus = pScoring->_pepMassInfo.dExpPepMass
         - pScoring->_pepMassInfo.dPeptideMassTolerance - 3.0 * C13_DIFF * PROTON_MASS;

      pScoring->_pepMassInfo.dPeptideMassTolerancePlus = pScoring->_pepMassInfo.dExpPepMass
         + pScoring->_pepMassInfo.dPeptideMassTolerance;
   }
   else if (g_staticParams.tolerances.iIsotopeError == 4) // search -1 to +1, +2, +3 isotope windows
   {
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus = pScoring->_pepMassInfo.dExpPepMass
         - pScoring->_pepMassInfo.dPeptideMassTolerance - 3.0 * C13_DIFF * PROTON_MASS;

      pScoring->_pepMassInfo.dPeptideMassTolerancePlus = pScoring->_pepMassInfo.dExpPepMass
         + pScoring->_pepMassInfo.dPeptideMassTolerance + 1.0 * C13_DIFF * PROTON_MASS;
   }
   else if (g_staticParams.tolerances.iIsotopeError == 5) // search -8, -4, 0, 4, 8 windows
   {
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus = pScoring->_pepMassInfo.dExpPepMass
         - pScoring->_pepMassInfo.dPeptideMassTolerance - 8.1;

      pScoring->_pepMassInfo.dPeptideMassTolerancePlus = pScoring->_pepMassInfo.dExpPepMass
         + pScoring->_pepMassInfo.dPeptideMassTolerance + 8.1;
   }
   else  // Should not get here.
   {
      fprintf(stderr,  " Error - iIsotopeError=%d\n",  g_staticParams.tolerances.iIsotopeError);
      return false;
   }

   if (g_staticParams.vectorMassOffsets.size() > 0)
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus -= g_staticParams.vectorMassOffsets[g_staticParams.vectorMassOffsets.size()-1];

   if (pScoring->_pepMassInfo.dPeptideMassTolerancePlus > g_staticParams.options.dPeptideMassHigh)
      pScoring->_pepMassInfo.dPeptideMassTolerancePlus = g_staticParams.options.dPeptideMassHigh;

   if (pScoring->_pepMassInfo.dPeptideMassToleranceMinus < g_staticParams.options.dPeptideMassLow)
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus = g_staticParams.options.dPeptideMassLow;

   if (pScoring->_pepMassInfo.dPeptideMassToleranceMinus < 100.0)
      pScoring->_pepMassInfo.dPeptideMassToleranceMinus = 100.0;

   return true;
}


//  Reads MSMS data file as ASCII mass/intensity pairs.
bool xlinkx_preprocess::LoadIons(struct Query *pScoring,
                               double *pdTmpRawData,
                               Spectrum mstSpectrum,
                               struct PreprocessStruct *pPre)
{
   int  i;
   double dIon,
          dIntensity;

   i = 0;
   while(true)
   {
      if (i >= mstSpectrum.size())
         break;

      dIon = mstSpectrum.at(i).mz;
      dIntensity = mstSpectrum.at(i).intensity;
      i++;

      pScoring->_spectrumInfoInternal.dTotalIntensity += dIntensity;

      if ((dIntensity >= g_staticParams.options.dMinIntensity) && (dIntensity > 0.0))
      {
         if (dIon < (pScoring->_pepMassInfo.dExpPepMass + 50.0))
         {
            int iBinIon = BIN(dIon);

            dIntensity = sqrt(dIntensity);

            if (iBinIon > pPre->iHighestIon)
               pPre->iHighestIon = iBinIon;

            if ((iBinIon < pScoring->_spectrumInfoInternal.iArraySize)
                  && (dIntensity > pdTmpRawData[iBinIon]))
            {
               if (g_staticParams.options.iRemovePrecursor == 1)
               {
                  double dMZ = (pScoring->_pepMassInfo.dExpPepMass
                        + (pScoring->_spectrumInfoInternal.iChargeState - 1) * PROTON_MASS)
                     / (double)(pScoring->_spectrumInfoInternal.iChargeState);

                  if (fabs(dIon - dMZ) > g_staticParams.options.dRemovePrecursorTol)
                  {
                     if (dIntensity > pdTmpRawData[iBinIon])
                        pdTmpRawData[iBinIon] = dIntensity;

                     if (pdTmpRawData[iBinIon] > pPre->dHighestIntensity)
                        pPre->dHighestIntensity = pdTmpRawData[iBinIon];
                  }
               }
               else if (g_staticParams.options.iRemovePrecursor == 2)
               {
                  int j;
                  int bNotPrec = 1;

                  for (j=1; j <= pScoring->_spectrumInfoInternal.iChargeState; j++)
                  {
                     double dMZ;

                     dMZ = (pScoring->_pepMassInfo.dExpPepMass + (j - 1)*PROTON_MASS) / (double)(j);
                     if (fabs(dIon - dMZ) < g_staticParams.options.dRemovePrecursorTol)
                     {
                        bNotPrec = 0;
                        break;
                     }
                  }
                  if (bNotPrec)
                  {
                     if (dIntensity > pdTmpRawData[iBinIon])
                        pdTmpRawData[iBinIon] = dIntensity;

                     if (pdTmpRawData[iBinIon] > pPre->dHighestIntensity)
                        pPre->dHighestIntensity = pdTmpRawData[iBinIon];
                  }
               }
               else // iRemovePrecursor==0
               {
                  if (dIntensity > pdTmpRawData[iBinIon])
                     pdTmpRawData[iBinIon] = dIntensity;

                  if (pdTmpRawData[iBinIon] > pPre->dHighestIntensity)
                     pPre->dHighestIntensity = pdTmpRawData[iBinIon];
               }
            }
         }
      }
   }

   return true;
}


// pdTmpRawData now holds raw data, pdTmpCorrelationData is windowed data after this function
void xlinkx_preprocess::MakeCorrData(double *pdTmpRawData,
                                   double *pdTmpCorrelationData,
                                   struct Query *pScoring,
                                   struct PreprocessStruct *pPre)
{
   int  i,
        ii,
        iBin,
        iWindowSize,
        iNumWindows=10;
   double dMaxWindowInten,
          dTmp1,
          dTmp2;

   iWindowSize = (int)((pPre->iHighestIon)/iNumWindows) + 1;

   for (i=0; i<iNumWindows; i++)
   {
      dMaxWindowInten = 0.0;

      for (ii=0; ii<iWindowSize; ii++)    // Find max inten. in window.
      {
         iBin = i*iWindowSize+ii;
         if (iBin < pScoring->_spectrumInfoInternal.iArraySize)
         {
            if (pdTmpRawData[iBin] > dMaxWindowInten)
               dMaxWindowInten = pdTmpRawData[iBin];
         }
      }

      if (dMaxWindowInten > 0.0)
      {
         dTmp1 = 50.0 / dMaxWindowInten;
         dTmp2 = 0.05 * pPre->dHighestIntensity;

         for (ii=0; ii<iWindowSize; ii++)    // Normalize to max inten. in window.
         {
            iBin = i*iWindowSize+ii;
            if (iBin < pScoring->_spectrumInfoInternal.iArraySize)
            {
               if (pdTmpRawData[iBin] > dTmp2)
                  pdTmpCorrelationData[iBin] = pdTmpRawData[iBin]*dTmp1;
            }
         }
      }
   }
}


// Smooth input data over 5 points.
bool xlinkx_preprocess::Smooth(double *data,
                             int iArraySize,
                             double *pdTmpSmoothedSpectrum)
{
   int  i;

   data[0] = 0.0;
   data[1] = 0.0;
   data[iArraySize-1] = 0.0;
   data[iArraySize-2] = 0.0;

   for (i=2; i<iArraySize-2; i++)
   {
      // *0.0625 is same as divide by 16.
      pdTmpSmoothedSpectrum[i] = (data[i-2]+4.0*data[i-1]+6.0*data[i]+4.0*data[i+1]+data[i+2]) * 0.0625;
   }

   memcpy(data, pdTmpSmoothedSpectrum, iArraySize*sizeof(double));

   return true;
}


// Run 2 passes through to pull out peaks.
bool xlinkx_preprocess::PeakExtract(double *data,
                                  int iArraySize,
                                  double *pdTmpPeakExtracted)
{
   int  i,
        ii,
        iStartIndex,
        iEndIndex;
   double dStdDev,
          dAvgInten;

   // 1st pass, choose only peak greater than avg + dStdDev.
   for (i=0; i<iArraySize; i++)
   {
      pdTmpPeakExtracted[i] = 0.0;
      dAvgInten = 0.0;

      iStartIndex = i-50;
      if (i-50 < 0)
         iStartIndex = 0;

      iEndIndex = i+50;
      if (i+50 > iArraySize-1)
         iEndIndex = iArraySize-1;

      for (ii=iStartIndex; ii<=iEndIndex; ii++)
         dAvgInten += (double)data[ii];
      dAvgInten /= iEndIndex-iStartIndex;

      dStdDev = 0.0;
      for (ii=iStartIndex; ii<=iEndIndex; ii++)
         dStdDev += (data[ii]-dAvgInten)*(data[ii]-dAvgInten);
      dStdDev = sqrt(dStdDev/(iEndIndex-iStartIndex+1));

      if ((i > 0) && (i < iArraySize-1))
      {
         if (data[i] > (dAvgInten+dStdDev))
         {
            pdTmpPeakExtracted[i] = data[i] - dAvgInten + dStdDev;
            data[i] = 0;     // Remove the peak before 2nd pass.
         }
      }
   }

   // 2nd pass, choose only peak greater than avg + 2*dStdDev.
   for (i=0; i<iArraySize; i++)
   {
      dAvgInten = 0.0;

      iStartIndex = i-50;
      if (i-50 < 0)
         iStartIndex = 0;

      iEndIndex = i+50;
      if (i+50 > iArraySize-1)
         iEndIndex = iArraySize-1;

      for (ii=iStartIndex; ii<=iEndIndex; ii++)
         dAvgInten += (double)data[ii];
      dAvgInten /= iEndIndex-iStartIndex;

      dStdDev = 0.0;
      for (ii=iStartIndex; ii<=iEndIndex; ii++)
         dStdDev += (data[ii]-dAvgInten)*(data[ii]-dAvgInten);
      dStdDev = sqrt(dStdDev/(iEndIndex-iStartIndex+1));

      if ((i > 0) && (i < iArraySize-1))
      {
         if (data[i] > (dAvgInten + 2*dStdDev))
            pdTmpPeakExtracted[i] = data[i] - dAvgInten + dStdDev;
      }
   }

   memcpy(data, pdTmpPeakExtracted, (size_t)iArraySize*sizeof(double));

   return true;
}


// Pull out top # ions for intensity matching in search.
void xlinkx_preprocess::GetTopIons(double *pdTmpRawData,
                                 struct msdata *pTmpSpData,
                                 int iArraySize)
{
   int  i,
        ii,
        iLowestIntenIndex=0;
   double dLowestInten=0.0,
          dMaxInten=0.0;

   for (i=0; i<iArraySize; i++)
   {
      if (pdTmpRawData[i] > dLowestInten)
      {
         (pTmpSpData+iLowestIntenIndex)->dIntensity = (double)pdTmpRawData[i];
         (pTmpSpData+iLowestIntenIndex)->dIon = (double)i;

         if ((pTmpSpData+iLowestIntenIndex)->dIntensity > dMaxInten)
            dMaxInten = (pTmpSpData+iLowestIntenIndex)->dIntensity;

         dLowestInten = (pTmpSpData+0)->dIntensity;
         iLowestIntenIndex = 0;

         for (ii=1; ii<NUM_SP_IONS; ii++)
         {
            if ((pTmpSpData+ii)->dIntensity < dLowestInten)
            {
               dLowestInten = (pTmpSpData+ii)->dIntensity;
               iLowestIntenIndex=ii;
            }
         }
      }
   }

   if (dMaxInten > FLOAT_ZERO)
   {
      for (i=0; i<NUM_SP_IONS; i++)
         (pTmpSpData+i)->dIntensity = (((pTmpSpData+i)->dIntensity)/dMaxInten)*100.0;
   }
}


int xlinkx_preprocess::QsortByIon(const void *p0, const void *p1)
{
   if ( ((struct msdata *) p1)->dIon < ((struct msdata *) p0)->dIon )
      return (1);
   else if ( ((struct msdata *) p1)->dIon > ((struct msdata *) p0)->dIon )
      return (-1);
   else
      return (0);
}


// Works on Sp data.
void xlinkx_preprocess::StairStep(struct msdata *pTmpSpData)
{
   int  i,
        ii,
        iii;
   double dMaxInten,
          dGap;

   i=0;
   while (i < NUM_SP_IONS-1)
   {
      ii = i;
      dMaxInten = (pTmpSpData+i)->dIntensity;
      dGap = 0.0;

      while (dGap<=g_staticParams.tolerances.dFragmentBinSize && ii<NUM_SP_IONS-1)
      {
         ii++;
         dGap = (pTmpSpData+ii)->dIon - (pTmpSpData+ii-1)->dIon;

         // Finds the max intensity for adjacent points.
         if (dGap<=g_staticParams.tolerances.dFragmentBinSize)
         {
            if ((pTmpSpData+ii)->dIntensity > dMaxInten)
               dMaxInten = (pTmpSpData+ii)->dIntensity;
         }
      }

      // Sets the adjacent points to the dMaxInten.
      for (iii=i; iii<ii; iii++)
         (pTmpSpData+iii)->dIntensity = dMaxInten;

      i = ii;
   }
}


//MH: This function allocates memory to be shared by threads for spectral processing
bool xlinkx_preprocess::AllocateMemory(int maxNumThreads)
{
   int i;

   //MH: Find appropriately sized array cushion based on user parameters. Fixes error found by Patrick Pedrioli for
   // very wide mass tolerance searches (i.e. 500 Da).
   double dCushion = 0.0;
   if (g_staticParams.tolerances.iMassToleranceUnits == 0) // amu
   {
      dCushion = g_staticParams.tolerances.dInputTolerance;

      if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
      {
        dCushion *= 8; //MH: hope +8 is large enough charge because g_staticParams.options.iEndCharge can be overridden.
      }
   }
   else if (g_staticParams.tolerances.iMassToleranceUnits == 1) // mmu
   {
      dCushion = g_staticParams.tolerances.dInputTolerance * 0.001;

      if (g_staticParams.tolerances.iMassToleranceType == 1)  // precursor m/z tolerance
      {
         dCushion *= 8; //MH: hope +8 is large enough charge because g_staticParams.options.iEndCharge can be overridden.
      }
   }
   else // ppm
   {
      dCushion = g_staticParams.tolerances.dInputTolerance * g_staticParams.options.dPeptideMassHigh / 1000000.0;
   }

   //MH: Must be equal to largest possible array
   int iArraySize = (int)((g_staticParams.options.dPeptideMassHigh + dCushion + 2.0) * g_staticParams.dInverseBinWidth);

   //MH: Initally mark all arrays as available (i.e. false=not inuse).
   pbMemoryPool = new bool[maxNumThreads];
   for (i=0; i<maxNumThreads; i++)
   {
      pbMemoryPool[i] = false;
   }

   //MH: Allocate arrays
   ppdTmpRawDataArr = new double*[maxNumThreads]();
   for (i=0; i<maxNumThreads; i++)
   {
      try
      {
         ppdTmpRawDataArr[i] = new double[iArraySize]();
      }
      catch (std::bad_alloc& ba)
      {
         fprintf(stderr,  " Error - new(pdTmpRawData[%d]). bad_alloc: %s.\n", iArraySize, ba.what());
         fprintf(stderr, "Xlinkx ran out of memory. Look into \"spectrum_batch_size\"\n");
         fprintf(stderr, "parameters to address mitigate memory use.\n");
         return false;
      }
   }

   //MH: Allocate arrays
   ppdTmpFastXcorrDataArr = new double*[maxNumThreads]();
   for (i=0; i<maxNumThreads; i++)
   {
      try
      {
         ppdTmpFastXcorrDataArr[i] = new double[iArraySize]();
      }
      catch (std::bad_alloc& ba)
      {
         fprintf(stderr,  " Error - new(pdTmpFastXcorrData[%d]). bad_alloc: %s.\n", iArraySize, ba.what());
         fprintf(stderr, "Xlinkx ran out of memory. Look into \"spectrum_batch_size\"\n");
         fprintf(stderr, "parameters to address mitigate memory use.\n");
         return false;
      }
   }

   //MH: Allocate arrays
   ppdTmpCorrelationDataArr = new double*[maxNumThreads]();
   for (i=0; i<maxNumThreads; i++)
   {
      try
      {
         ppdTmpCorrelationDataArr[i] = new double[iArraySize]();
      }
      catch (std::bad_alloc& ba)
      {
         fprintf(stderr,  " Error - new(pdTmpCorrelationData[%d]). bad_alloc: %s.\n", iArraySize, ba.what());
         fprintf(stderr, "Xlinkx ran out of memory. Look into \"spectrum_batch_size\"\n");
         fprintf(stderr, "parameters to address mitigate memory use.\n");
         return false;
      }
   }

   //MH: Allocate arrays
   ppdTmpSmoothedSpectrumArr = new double*[maxNumThreads]();
   for (i=0; i<maxNumThreads; i++)
   {
      try
      {
         ppdTmpSmoothedSpectrumArr[i] = new double[iArraySize]();
      }
      catch (std::bad_alloc& ba)
      {
         fprintf(stderr,  " Error - new(pdTmpSmoothedSpectrum[%d]). bad_alloc: %s.\n", iArraySize, ba.what());
         fprintf(stderr, "Xlinkx ran out of memory. Look into \"spectrum_batch_size\"\n");
         fprintf(stderr, "parameters to address mitigate memory use.\n");
         return false;
      }
   }

   //MH: Allocate arrays
   ppdTmpPeakExtractedArr = new double*[maxNumThreads]();
   for (i=0; i<maxNumThreads; i++)
   {
      try
      {
         ppdTmpPeakExtractedArr[i] = new double[iArraySize]();
      }
      catch (std::bad_alloc& ba)
      {
         fprintf(stderr,  " Error - new(pdTmpSmoothedSpectrum[%d]). bad_alloc: %s.\n", iArraySize, ba.what());
         fprintf(stderr, "Xlinkx ran out of memory. Look into \"spectrum_batch_size\"\n");
         fprintf(stderr, "parameters to address mitigate memory use.\n");
         return false;
      }
   }

   return true;
}


//MH: Deallocates memory shared by threads during spectral processing.
bool xlinkx_preprocess::DeallocateMemory(int maxNumThreads)
{
   int i;

   delete [] pbMemoryPool;

   for (i=0; i<maxNumThreads; i++)
   {
      delete [] ppdTmpRawDataArr[i];
      delete [] ppdTmpFastXcorrDataArr[i];
      delete [] ppdTmpCorrelationDataArr[i];
      delete [] ppdTmpSmoothedSpectrumArr[i];
      delete [] ppdTmpPeakExtractedArr[i];
   }

   delete [] ppdTmpRawDataArr;
   delete [] ppdTmpFastXcorrDataArr;
   delete [] ppdTmpCorrelationDataArr;
   delete [] ppdTmpSmoothedSpectrumArr;
   delete [] ppdTmpPeakExtractedArr;

   return true;
}

bool xlinkx_preprocess::IsValidInputType(int inputType)
{
   return (inputType == InputType_MZXML || inputType == InputType_RAW);
}
