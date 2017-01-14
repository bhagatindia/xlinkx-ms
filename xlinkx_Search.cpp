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
#include "xlinkx.h"
#include "xlinkx_Search.h"
#include "xlinkx_DataInternal.h"
#include "xlinkx_Preprocess.h"
#include "CometDecoys.h"

// Generate data for both sp scoring (pfSpScoreData) and xcorr analysis (FastXcorr).
xlinkx_Search::xlinkx_Search()
{
}
#define NUMPEPTIDES 10

xlinkx_Search::~xlinkx_Search()
{
}

void insert_pep_pq(char *pepArray[], float xcorrArray[], char *ins_pep, float ins_xcorr)
{
   int i; 

   // check for duplicates
   for (i = 0; i< NUMPEPTIDES - 1; i++) {
      if (pepArray[i] != NULL) {
         if (!strcmp(pepArray[i], ins_pep)) {
             //cout << "" << pepArray[i] << " " << ins_pep << " Found duplicate" << endl;
             return;
         }
      }
   }

   // Insert first into the array
   if (xcorrArray[NUMPEPTIDES - 1] < ins_xcorr) {
      xcorrArray[NUMPEPTIDES - 1] = ins_xcorr;
      pepArray[NUMPEPTIDES - 1] = ins_pep;
   } else return;

   // shiffle
   for (i = NUMPEPTIDES - 1; i > 0; i--) {
      if (xcorrArray[i] > xcorrArray[i-1]) {
         // Swap
         float temp = xcorrArray[i];
         char *temp_pep = pepArray[i];
         xcorrArray[i] = xcorrArray[i-1];
         pepArray[i] = pepArray[i-1];
         xcorrArray[i-1] = temp;
         pepArray[i-1] = temp_pep;
      } else break;
   }
}

#define BIN_SIZE 0.1
#define MAX_XCORR_VALUE 20
#define NUM_BINS (int)(MAX_XCORR_VALUE/BIN_SIZE + 1)

inline int xlinkx_get_histogram_bin_num(float value)
{
    if (value > MAX_XCORR_VALUE)
       value = MAX_XCORR_VALUE;
    else if (value < 0)
       value = 0;
    return value/BIN_SIZE;
}

void xlinkx_print_histogram(int hist_pep[])
{
   for (int i = 0; i <NUM_BINS; i++)
      printf("%d ", hist_pep[i]);
   printf("\n");
}

#define DECOY_SIZE 3000
#define MAX_DECOY_PEP_LEN 40

bool xlinkx_Search::CalculateEValue(int *hist_pep,
                                    int iMatchPepCount,
                                    double *dSlope,
                                    double *dIntercept,
                                    double dNeutralPepMass,
                                    int iScanNumber)
{
   int iMaxCorr;
   int iStartCorr;
   int iNextCorr;

   if (iMatchPepCount < DECOY_SIZE)
   {
      if (!GenerateXcorrDecoys(dNeutralPepMass, iMatchPepCount, hist_pep, iScanNumber))
      {
         return false;
      }
   }

   LinearRegression(hist_pep, dSlope, dIntercept, &iMaxCorr, &iStartCorr, &iNextCorr);
   *dSlope *= 10.0; // Used in pow() function so do multiply outside of for loop.

   return true;
}


void xlinkx_Search::LinearRegression(int *piHistogram,
                                     double *slope,
                                     double *intercept,
                                     int *iMaxXcorr,
                                     int *iStartXcorr,
                                     int *iNextXcorr)
{
   double Sx, Sxy;      // Sum of square distances.
   double Mx, My;       // means
   double b, a;
   double SumX, SumY;   // Sum of X and Y values to calculate mean.

   double dCummulative[HISTO_SIZE];  // Cummulative frequency at each xcorr value.

   int i;
   int iNextCorr;    // 2nd best xcorr index
   int iMaxCorr=0;   // max xcorr index
   int iStartCorr;
   int iNumPoints;

   // Find maximum correlation score index.
   for (i=HISTO_SIZE-2; i>=0; i--)
   {
      if (piHistogram[i] > 0)
         break;
   }
   iMaxCorr = i;

   iNextCorr = 0;
   for (i=0; i<iMaxCorr; i++)
   {
      if (piHistogram[i]==0)
      {
         // register iNextCorr if there's a histo value of 0 consecutively
         if (piHistogram[i+1]==0 || i+1 == iMaxCorr)
         {
            if (i>0)
               iNextCorr = i-1;
            break;
         }
      }
   }

   if (i==iMaxCorr)
   {
      iNextCorr = iMaxCorr;
      if (iMaxCorr>12)
         iNextCorr = iMaxCorr-2;
   }


   // Create cummulative distribution function from iNextCorr down, skipping the outliers.
   dCummulative[iNextCorr] = piHistogram[iNextCorr];
   for (i=iNextCorr-1; i>=0; i--)
   {
      dCummulative[i] = dCummulative[i+1] + piHistogram[i];
      if (piHistogram[i+1] == 0)
         dCummulative[i+1] = 0.0;
   }

   // log10
   for (i=iNextCorr; i>=0; i--)
   {
      piHistogram[i] = (int)dCummulative[i];  // First store cummulative in histogram.
      dCummulative[i] = log10(dCummulative[i]);
   }

   iStartCorr = 1;
   if (iNextCorr >= 30)
      iStartCorr = (int)(iNextCorr - iNextCorr*0.25);
   else if (iNextCorr >= 15)
      iStartCorr = (int)(iNextCorr - iNextCorr*0.5);

   Mx=My=a=b=0.0;

   while (iStartCorr >= 0)
   {
      Sx=Sxy=SumX=SumY=0.0;
      iNumPoints=0;

      // Calculate means.
      for (i=iStartCorr; i<=iNextCorr; i++)
      {
         if (piHistogram[i] > 0)
         {
            SumY += (float)dCummulative[i];
            SumX += i;
            iNumPoints++;
         }
      }

      if (iNumPoints > 0)
      {
         Mx = SumX / iNumPoints;
         My = SumY / iNumPoints;
      }
      else
         Mx = My = 0.0;

      // Calculate sum of squares.
      for (i=iStartCorr; i<=iNextCorr; i++)
      {
         if (dCummulative[i] > 0)
         {
            double dX;
            double dY;

            dX = i - Mx;
            dY = dCummulative[i] - My;

            Sx  += dX*dX;
            Sxy += dX*dY;
         }
      }

      if (Sx > 0)
         b = Sxy / Sx;   // slope
      else
         b = 0;

      if (b < 0.0)
         break;
      else
         iStartCorr--;
   }

   a = My - b*Mx;  // y-intercept

   *slope = b;
   *intercept = a;
   *iMaxXcorr = iMaxCorr;
   *iStartXcorr = iStartCorr;
   *iNextXcorr = iNextCorr;
}


// Make synthetic decoy spectra to fill out correlation histogram by going
// through each candidate peptide and rotating spectra in m/z space.
bool xlinkx_Search::GenerateXcorrDecoys(double dNeutralPepMass,
                                        int iMatchPepCount,
                                        int *hist_pep,
                                        int iScanNumber)
{
   int i;
   int ii;
   int j;
   int bin_num;
   int iMaxFragCharge;
   int ctCharge;
   double dBion;
   double dYion;
   double dFastXcorr;
   double dFragmentIonMass;

   int *piHistogram;

   int iFragmentIonMass;
   int iWhichQuery;

   for (iWhichQuery=0; iWhichQuery<(int)g_pvQuery.size(); iWhichQuery++)
   {
      if (g_pvQuery.at(iWhichQuery)->_spectrumInfoInternal.iScanNumber == iScanNumber)
         break;
   }
   if (iWhichQuery < (int)g_pvQuery.size() && g_pvQuery.at(iWhichQuery)->_spectrumInfoInternal.iScanNumber == iScanNumber)
   {
      Query* pQuery = g_pvQuery.at(iWhichQuery);

      piHistogram = hist_pep;

      //iMaxFragCharge = pQuery->_spectrumInfoInternal.iMaxFragCharge;
      iMaxFragCharge = 1;  //FIX only considering 1+ charges now

      // DECOY_SIZE is the minimum # of decoys required or else this function is
      // called.  So need generate iLoopMax more xcorr scores for the histogram.
      int iLoopMax = DECOY_SIZE - iMatchPepCount;
      int iLastEntry;

      iLastEntry = iMatchPepCount;

      if (iLastEntry > g_staticParams.options.iNumStored)
         iLastEntry = g_staticParams.options.iNumStored;

      j=0;
      for (i=0; i<iLoopMax; i++)  // iterate through required # decoys
      {
         dFastXcorr = 0.0;

         for (j=0; j<MAX_DECOY_PEP_LEN; j++)  // iterate through decoy fragment ions
         {
            dBion = decoyIons[i].pdIonsN[j];
            dYion = decoyIons[i].pdIonsC[j];

            for (ii=0; ii<2; ii++)
            {
               dFragmentIonMass =  0.0;
               switch (ii)
               {
                  case 0:
                     dFragmentIonMass = dBion;
                     break;
                  case 1:
                     dFragmentIonMass = dYion;
                     break;
               }

               for (ctCharge=1; ctCharge<=iMaxFragCharge; ctCharge++)
               {
                  dFragmentIonMass = (dFragmentIonMass + (ctCharge-1)*PROTON_MASS)/ctCharge;

                  if (dFragmentIonMass < dNeutralPepMass)
                  {
                     iFragmentIonMass = BIN(dFragmentIonMass);

                     if (iFragmentIonMass < pQuery->_spectrumInfoInternal.iArraySize && iFragmentIonMass >= 0)
                     {
                        int x = iFragmentIonMass / SPARSE_MATRIX_SIZE;
                        if (pQuery->ppfSparseFastXcorrData[x]!=NULL)
                        {
                           int y = iFragmentIonMass - (x*SPARSE_MATRIX_SIZE);
                           dFastXcorr += pQuery->ppfSparseFastXcorrData[x][y];
                        }
                     }
                     else
                     {
                        char szErrorMsg[256];
                        sprintf(szErrorMsg,  " Error - XCORR DECOY: dFragMass %f, iFragMass %d, ArraySize %d, InputMass %f, scan %d, z %d",
                              dFragmentIonMass,
                              iFragmentIonMass,
                              pQuery->_spectrumInfoInternal.iArraySize,
                              pQuery->_pepMassInfo.dExpPepMass,
                              pQuery->_spectrumInfoInternal.iScanNumber,
                              ctCharge);

                        string strErrorMsg(szErrorMsg);
                        logerr(szErrorMsg);
                        return false;
                     }
                  }

               }
            }
         }

         dFastXcorr *= 0.005;
         bin_num = xlinkx_get_histogram_bin_num(dFastXcorr);
         piHistogram[bin_num] += 1;
      }

      return true;
   }

   return false;

}


void xlinkx_Search::SearchForPeptides(char *szMZXML,
                                      const char *protein_file,
                                      enzyme_cut_params params,
                                      const char *pep_hash_file,
                                      const char *pep_output_file)
{
   int i;
   int ii;
   double dTolerance;
   double dPPM = 20.0;  // use 20ppm tolerance for now
   int hist_pep1[NUM_BINS],
       hist_pep2[NUM_BINS],
       hist_combined[NUM_BINS],
       num_pep1,
       num_pep2,
       num_pep_combined;

   char szProt1[128];
   char szProt2[128];

#define LYSINE_MOD 197.032422

   char *toppep1[NUMPEPTIDES], *toppep2[NUMPEPTIDES], *toppepcombined[NUMPEPTIDES];
   float xcorrPep1[NUMPEPTIDES], xcorrPep2[NUMPEPTIDES], xcorrCombined[NUMPEPTIDES];

   FILE *fptxt;
   if ((fptxt=fopen(pep_output_file, "w")) == NULL)
   {
      printf(" Error - cannot write txt output %s\n", pep_output_file);
      return;
   }

   MSReader mstReader;
   Spectrum mstSpectrum;
   // We want to read only MS2 scans.
   vector<MSSpectrumType> msLevel;
   msLevel.push_back(MS2);

   mstReader.setFilter(msLevel);
   mstReader.readFile(szMZXML, mstSpectrum, 1);
 
   xlinkx_preprocess::AllocateMemory(1);

   // If PeptideHash not present, generate it now; otherwise open the hash file.
   protein_hash_db_t phdp = phd_retrieve_hash_db(protein_file, params, pep_hash_file);

   fprintf(fptxt, "scan\texp_mass1\texp_mass2\tpeptide1\txcorr1\tevalue1\tcalcmass1\tpeptide2\txcorr2\tevalue2\tcalcmass2\tcombinedxcorr\tcombinedevalue\n"); 

   FILE *fpxml;
   char szOutput[1024];
   char szBaseName[1024];

   strcpy(szBaseName, szMZXML);
   if (!strcmp(szBaseName+strlen(szBaseName)-6, ".mzXML"))
      szBaseName[strlen(szBaseName)-6]='\0';
   else if (!strcmp(szBaseName+strlen(szBaseName)-6, ".mzML"))
      szBaseName[strlen(szBaseName)-5]='\0';
   sprintf(szOutput, "%s.pep.xml", szBaseName);

   if ((fpxml=fopen(szOutput, "w")) == NULL)
   {
      printf(" Error - cannot write pepXML output %s\n", szOutput);
      return;
   }
   int iIndex=0;

   WritePepXMLHeader(fpxml, szBaseName);

   if (!g_staticParams.options.bVerboseOutput)
   {
      printf("percent complete: ");
      fflush(stdout);
   }

   for (i=0; i<(int)pvSpectrumList.size(); i++)
   {
      for (ii=0; ii<(int)pvSpectrumList.at(i).pvdPrecursors.size(); ii++)
      {
         for (int j = 0; j < NUM_BINS; j++)
            hist_pep1[j] = hist_pep2[j] = hist_combined[j] = 0;

         num_pep1 = num_pep2 = num_pep_combined = 0;

         double dMZ1 =  (pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1
               + pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge1 * PROTON_MASS)/ pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge1;
         double dMZ2 =  (pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2
               + pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge2 * PROTON_MASS)/ pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge2;

         xlinkx_preprocess::LoadAndPreprocessSpectra(mstReader, pvSpectrumList.at(i).iScanNumber, dMZ1, dMZ2);
 
         for (int li = 0; li < NUMPEPTIDES; li++) {
            xcorrPep1[li] = xcorrPep2[li] = xcorrCombined[li] = -99999;
            toppep1[li] = toppep2[li] = toppepcombined[li] = NULL;
         }

         if (g_staticParams.options.bVerboseOutput)
         {
            printf("Scan %d, retrieving peptides of mass %0.4f (%d+ %0.4f) and %0.4f (%d+ %0.4f)\n",
                  pvSpectrumList.at(i).iScanNumber,
                  pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1,
                  pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge1,
                  dMZ1,
                  pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2,
                  pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge2,
                  dMZ2);
         }

         double pep_mass1 = pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1 - LYSINE_MOD - g_staticParams.precalcMasses.dOH2;
         double pep_mass2 = pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2 - LYSINE_MOD - g_staticParams.precalcMasses.dOH2;

         if (pep_mass1 <= 0)
         {
            cout << "Peptide mass1 is coming out to be zero after removing Lysine residue" << endl;
            exit(1);
         }

         if (pep_mass2 <= 0)
         {
            cout << "Peptide mass2 is coming out to be zero after removing Lysine residue" << endl;
            exit(1);
         }

         fprintf(fptxt, "%d\t%f\t%f", pvSpectrumList.at(i).iScanNumber, pep_mass1, pep_mass2);

         double dXcorr = 0.0;
         vector<double> vdXcorr_pep1;
         vector<double> vdXcorr_pep2;

         if (g_staticParams.options.bVerboseOutput)
         {
            cout << "After Lysine residue reduction the peptide of mass " << pep_mass1 << " are being extracted";
            cout << " (" << pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1 << ")" << endl;
         }
         dTolerance = (dPPM * pep_mass1) / 1e6;
         //vector<string*> *peptides = phdp->phd_get_peptides_ofmass(pep_mass1);
         vector<string*> *peptides1 = phdp->phd_get_peptides_ofmass_tolerance(pep_mass1, 1.0);

         for (string *peptide : *peptides1)
         {
            char *szPeptide = new char[(*peptide).length() + 1];
            strcpy(szPeptide, (*peptide).c_str() );

            // sanity check to ignore peptides w/unknown AA residues
            // should not be needed now that this is addressed in the hash building
            if (strchr(szPeptide, 'B') || strchr(szPeptide, 'X') || strchr(szPeptide, 'J') || strchr(szPeptide, 'Z'))
               dXcorr = 0.0;
            else
               dXcorr = XcorrScore(szPeptide, pvSpectrumList.at(i).iScanNumber);

            vdXcorr_pep1.push_back(dXcorr);

            hist_pep1[xlinkx_get_histogram_bin_num(dXcorr)]++;
            insert_pep_pq(toppep1, xcorrPep1, szPeptide, dXcorr);
            num_pep1++;
            if (g_staticParams.options.bVerboseOutput)
               cout << "pep1: " << *peptide << "  xcorr " << dXcorr << endl;
         }

         if (g_staticParams.options.bVerboseOutput)
         {
            cout << "After Lysine residue reduction the peptide of mass " << pep_mass2 << " are being extracted";
            cout << " (" << pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2 << ")" << endl;
         }
         dTolerance = (dPPM * pep_mass2) / 1e6;
         //peptides = phdp->phd_get_peptides_ofmass(pep_mass2);
         vector<string*> *peptides2 = phdp->phd_get_peptides_ofmass_tolerance(pep_mass2, 1.0);

         for (string *peptide : *peptides2)
         {
            char *szPeptide = new char[(*peptide).length() + 1];
            strcpy(szPeptide, (*peptide).c_str() );

            // sanity check to ignore peptides w/unknown AA residues
            // should not be needed now that this is addressed in the hash building
            if (strchr(szPeptide, 'B') || strchr(szPeptide, 'X') || strchr(szPeptide, 'J') || strchr(szPeptide, 'Z'))
               dXcorr = 0.0;
            else
               dXcorr = XcorrScore(szPeptide, pvSpectrumList.at(i).iScanNumber);

            vdXcorr_pep2.push_back(dXcorr);

            hist_pep2[xlinkx_get_histogram_bin_num(dXcorr)]++;
            insert_pep_pq(toppep2, xcorrPep2, szPeptide, dXcorr);
            num_pep2++;
            if (g_staticParams.options.bVerboseOutput)
               cout << "pep2: " << *peptide << "  xcorr " << dXcorr << endl;
         }

         if (toppep1[0] == NULL || toppep2[0] == NULL)
            continue;

         double dSlope;
         double dIntercept;
         double dExpect;

         // return dSlope and dIntercept for histogram
         CalculateEValue(hist_pep1, num_pep1, &dSlope, &dIntercept,
               pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1, pvSpectrumList.at(i).iScanNumber);

         if (g_staticParams.options.bVerboseOutput)
         {
            xlinkx_print_histogram(hist_pep1);
            cout << "Top "<< NUMPEPTIDES << " pep1 peptides for this scan are " << endl;
         }

         double dExpect1 = 999;;
         for (int li = 0 ; li < NUMPEPTIDES; li++)
         {
            if (toppep1[li] != NULL)
            {
               if (dSlope > 0)
                  dExpect = 999;
               else
                  dExpect = pow(10.0, dSlope * xcorrPep1[li] + dIntercept);

               if (li == 0)
                  dExpect1 = dExpect;

               if (g_staticParams.options.bVerboseOutput)
                  cout << "pep1_top: " << toppep1[li] << " xcorr " << xcorrPep1[li] << " expect " << dExpect << endl;
            }
         }

         if (toppep1[0] != NULL)
            fprintf(fptxt, "\t%s\t%f\t%0.3E\t%f", toppep1[0], xcorrPep1[0], dExpect1, phdp->phd_calculate_mass_peptide(string(toppep1[0])));
         else
            fprintf(fptxt, "\t-\t0\t999\t0");


         CalculateEValue(hist_pep2, num_pep2, &dSlope, &dIntercept,
               pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2, pvSpectrumList.at(i).iScanNumber);

         if (g_staticParams.options.bVerboseOutput)
         {
            xlinkx_print_histogram(hist_pep2);
            cout << "Top "<< NUMPEPTIDES << " pep2 peptides for this scan are " << endl;
         }
         double dExpect2 = 999;
         for (int li = 0; li < NUMPEPTIDES; li++)
         {
            if (toppep2[li] != NULL)
            {
               if (dSlope > 0)
                  dExpect = 999;
               else
                  dExpect = pow(10.0, dSlope * xcorrPep2[li] + dIntercept);

               if (li == 0)
                  dExpect2 = dExpect;

               if (g_staticParams.options.bVerboseOutput)
                  cout << "pep2_top: " << toppep2[li] << " xcorr " << xcorrPep2[li] << " expect " << dExpect << endl;
            }
         }

         if (toppep2[0] != NULL)
            fprintf(fptxt, "\t%s\t%f\t%0.3E\t%f", toppep2[0], xcorrPep2[0], dExpect2, phdp->phd_calculate_mass_peptide(string(toppep2[0])));
         else
            fprintf(fptxt, "\t-\t0\t999\t0");

         if (g_staticParams.options.bVerboseOutput)
            cout << "Size of peptide1 list is " << num_pep1 << " and the size of peptide2 list is " << num_pep2 << endl;

         // Compute histogram of combined scores;
         for (int x=0; x<NUM_BINS; x++)
            hist_combined[x] = 0;
 
         for (vector<double>::iterator x = vdXcorr_pep1.begin(); x != vdXcorr_pep1.end(); ++x)
         {
            for (vector<double>::iterator y = vdXcorr_pep2.begin(); y != vdXcorr_pep2.end(); ++y)
            {
               hist_combined[xlinkx_get_histogram_bin_num(*x + *y)]++;
            }
         }

         CalculateEValue(hist_combined, num_pep_combined, &dSlope, &dIntercept,
               pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2, pvSpectrumList.at(i).iScanNumber);

         if (g_staticParams.options.bVerboseOutput)
            xlinkx_print_histogram(hist_combined);

         // take all combinations of top pep1 and pep2 and store best
         for (int x = 0; x< NUMPEPTIDES - 1; x++)
         {
            if (toppep1[x] != NULL)
            {
               for (int y = 0; y< NUMPEPTIDES - 1; y++)
               {
                  if (toppep2[y] != NULL)
                  {
                     char *combinedPep = new char[512];
                     sprintf(combinedPep, "%s + %s", toppep1[x], toppep2[y]);

                     double dCombinedXcorr = xcorrPep1[x] + xcorrPep2[y];

                     insert_pep_pq(toppepcombined, xcorrCombined, combinedPep, dCombinedXcorr);
                  }
               }
            }
         }

         double dExpectCombined = 999;
         for (int li = 0; li < NUMPEPTIDES; li++)
         {
            if (toppepcombined[li] != NULL)
            {
               if (dSlope > 0)
                  dExpect = 999;
               else
                  dExpect = pow(10.0, dSlope * xcorrCombined[li] + dIntercept);

               if (g_staticParams.options.bVerboseOutput)
                  cout << "combined: " << toppepcombined[li] << " xcorr " << xcorrCombined[li] << " expect " << dExpect << endl;

               if (li == 0)
               {
                  fprintf(fptxt, "\t%f\t%0.3E\n",  xcorrCombined[li], dExpect);
                  dExpectCombined = dExpect;
               }
            }
         }


         strcpy(szProt1, "prot1");
         strcpy(szProt2, "prot2");
      
         WriteSpectrumQuery(fpxml, szBaseName,
               pep_mass1, pep_mass2,
               xcorrPep1[0], xcorrPep2[0],
               dExpect1, dExpect2,
               phdp->phd_calculate_mass_peptide(string(toppep1[0])), phdp->phd_calculate_mass_peptide(string(toppep2[0])),
               xcorrCombined[0], dExpectCombined,
               toppep1[0], toppep2[0],
               szProt1, szProt2,
               pvSpectrumList.at(i).iPrecursorCharge, iIndex, pvSpectrumList.at(i).iScanNumber);
      }

      if (!g_staticParams.options.bVerboseOutput)
      {
         printf("%3d%%", (int)(100.0*i/pvSpectrumList.size()));
         fflush(stdout);
         printf("\b\b\b\b");
      }
   }

   fprintf(fpxml, "  </msms_run_summary>\n");
   fprintf(fpxml, "</msms_pipeline_analysis>\n");


   if (!g_staticParams.options.bVerboseOutput)
     printf("\n\nDone.\n");

   fclose(fptxt);
   fclose(fpxml);
}


double xlinkx_Search::XcorrScore(const char *szPeptide,
                                 int iScanNumber)
{

   int iWhichQuery;
   int iLenPeptide = strlen(szPeptide);
   double dXcorr = 0.0;

   // find which
   for (iWhichQuery=0; iWhichQuery<(int)g_pvQuery.size(); iWhichQuery++)
   {
      if (g_pvQuery.at(iWhichQuery)->_spectrumInfoInternal.iScanNumber == iScanNumber)
         break;
   }

   if (iWhichQuery < (int)g_pvQuery.size() && g_pvQuery.at(iWhichQuery)->_spectrumInfoInternal.iScanNumber == iScanNumber)
   {
      int bin, x, y;
      int iMax = g_pvQuery.at(iWhichQuery)->_spectrumInfoInternal.iArraySize/SPARSE_MATRIX_SIZE + 1;

      double dBion = g_staticParams.precalcMasses.dNtermProton;
      double dYion = g_staticParams.precalcMasses.dCtermOH2Proton;

      bool bBionLysine = false; // set to true after first b-ion lysine is modified
      bool bYionLysine = false; // set to true after first y-ion lysine is modified

      for (int i=0; i<iLenPeptide-1; i++) // will ignore multiple fragment ion charge states for now
      {
         dBion += g_staticParams.massUtility.pdAAMassFragment[(int)szPeptide[i]];
         if (szPeptide[i] == 'K' && !bBionLysine)
         {
            dBion += LYSINE_MOD;
            bBionLysine = true;
         }

         bin = BIN(dBion);
         x =  bin / SPARSE_MATRIX_SIZE;
         if (!(g_pvQuery.at(iWhichQuery)->ppfSparseFastXcorrData[x]==NULL || x>iMax)) // x should never be > iMax so this is just a safety check
         {
            y = bin - (x*SPARSE_MATRIX_SIZE);
            dXcorr += g_pvQuery.at(iWhichQuery)->ppfSparseFastXcorrData[x][y];
         }

         dYion += g_staticParams.massUtility.pdAAMassFragment[(int)szPeptide[iLenPeptide -1 - i]];
         if (szPeptide[iLenPeptide -1 - i] == 'K' && !bYionLysine)
         {
            dYion += LYSINE_MOD;
            bYionLysine = true;
         }

         bin = BIN(dYion);
         x =  bin / SPARSE_MATRIX_SIZE;
         if (!(g_pvQuery.at(iWhichQuery)->ppfSparseFastXcorrData[x]==NULL || x>iMax)) // x should never be > iMax so this is just a safety check
         {
            y = bin - (x*SPARSE_MATRIX_SIZE);
            dXcorr += g_pvQuery.at(iWhichQuery)->ppfSparseFastXcorrData[x][y];
         }
      }

      if (!bBionLysine && !bYionLysine) // sanity check
      {
         cout << " Error, no internal lysine: " << szPeptide << endl;
         exit(1);
      }

      dXcorr *= 0.005;
   }

   if (dXcorr < 0.0)
      dXcorr = 0.0;

   return dXcorr;
}


bool xlinkx_Search::WithinTolerance(double dMass1,
                                    double dMass2)
{
   if (dMass1>0.0 &&  (1E6 * fabs(dMass1 - dMass2)/dMass1) <= g_staticParams.tolerances.dInputTolerance)
      return true;
   else
      return false;
}


void xlinkx_Search::WritePepXMLHeader(FILE *fpxml,
                                      char *szBaseName)
{

   fprintf(fpxml, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
   fprintf(fpxml, "<msms_pipeline_analysis date= \"2016-11-16T15:59:37\" summary_xml=\"/net/pr/vol1/ProteomicsResource/search/engj/20161216-high-lo-xlinkx/%s.pep.xml\" xmlns=\"http://regis-web.systemsbiology.net/pepXML\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://regis-web.systemsbiology.net/pepXML http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v120.xsd\">\n", szBaseName);
   fprintf(fpxml, " <msms_run_summary base_name=\"/net/pr/vol1/ProteomicsResource/search/engj/20161216-high-lo-xlinkx/%s\" raw_data_type=\"raw\" raw_data=\".mzXML\">\n", szBaseName);
   fprintf(fpxml, "  <sample_enzyme name=\"trypsin\">\n");
   fprintf(fpxml, "   <specificity cut=\"KR\" no_cut=\"P\" sense=\"C\"/>\n");
   fprintf(fpxml, "  </sample_enzyme>\n");
   fprintf(fpxml, "  <search_summary base_name=\"/net/pr/vol1/ProteomicsResource/search/engj/20161216-high-lo-xlinkx/%s\" search_engine=\"Kojak\" search_engine_version=\"1.0.0\" precursor_mass_type=\"monoisotopic\" fragment_mass_type=\"monoisotopic\" search_id=\"1\">\n", szBaseName);
   fprintf(fpxml, "   <search_database local_path=\"/net/pr/vol1/ProteomicsResource/search/engj/20161216-high-lo-xlinkx/HUMAN.fasta.20160308\" type=\"AA\"/>\n");
   fprintf(fpxml, "   <enzymatic_search_contstraint enzyme=\"Trypsin\" max_num_internal_cleavages=\"1\" min_number_termini=\"2\"/>\n");
   fprintf(fpxml, "   <parameter name=\"ion_series_A\" value=\"0\"/>\n");
   fprintf(fpxml, "   <parameter name=\"ion_series_B\" value=\"1\"/>\n");
   fprintf(fpxml, "   <parameter name=\"ion_series_C\" value=\"0\"/>\n");
   fprintf(fpxml, "   <parameter name=\"ion_series_X\" value=\"0\"/>\n");
   fprintf(fpxml, "   <parameter name=\"ion_series_Y\" value=\"1\"/>\n");
   fprintf(fpxml, "   <parameter name=\"ion_series_Z\" value=\"0\"/>\n");
   fprintf(fpxml, "  </search_summary>\n");
}


void xlinkx_Search::WriteSpectrumQuery(FILE *fpxml,
                                       char *szBaseName,
                                       double dExpMass1,
                                       double dExpMass2,
                                       double dXcorr1,
                                       double dXcorr2,
                                       double dExpect1,
                                       double dExpect2,
                                       double dCalcMass1,
                                       double dCalcMass2,
                                       double dXcorrCombined,
                                       double dExpectCombined,
                                       char *szPep1,
                                       char *szPep2,
                                       char *szProt1,
                                       char *szProt2,
                                       int iCharge,
                                       int iIndex,
                                       int iScan) 
{                         
   int i;
   double xl = 300.0;
   
   fprintf(fpxml, "  <spectrum_query spectrum=\"%s.%d.%d.%d\" start_scan=\"%d\" end_scan=\"%d\" precursor_neutral_mass=\"%0.6f\" assumed_charge=\"%d\" index=\"%d\">\n",
         szBaseName, iScan, iScan, iCharge, iScan, iScan, dExpMass1+dExpMass2+xl, iCharge, ++iIndex);
   fprintf(fpxml, "   <search_result>\n");
   fprintf(fpxml, "    <search_hit hit_rank=\"1\" peptide=\"-\" peptide_prev_aa=\"-\" peptide_next_aa=\"-\" protein=\"-\" num_tot_proteins=\"1\" calc_neutral_pep_mass=\"%0.6f\" massdiff=\"%0.6f\" xlink_type=\"xl\">\n",
         dCalcMass1+dCalcMass2+xl, (dCalcMass1+dCalcMass2)-(dExpMass1+dExpMass2));
   fprintf(fpxml, "     <xlink identifier=\"BDP-NHP\" mass=\"200.00\">\n");
   fprintf(fpxml, "      <linked_peptide peptide=\"%s\" peptide_prev_aa=\"-\" peptide_next_aa=\"-\" protein=\"%s\" num_tot_proteins=\"1\" calc_neutral_pep_mass=\"%0.6f\" complement_mass=\"%0.6f\" designation=\"alpha\">\n",
         szPep1, szProt1, dCalcMass1, dCalcMass2+xl);
   fprintf(fpxml, "       <modification_info>\n");
   for (i=0; i<strlen(szPep1); i++)
      if (szPep1[i]=='K')
         break;
   fprintf(fpxml, "        <mod_aminoacid_mass position=\"%d\" mass=\"325.127385\"/>\n", i+1);
   fprintf(fpxml, "       </modification_info>\n");
   fprintf(fpxml, "       <xlink_score name=\"score\" value=\"%0.3E\"/>\n", dExpect1);
// fprintf(fpxml, "       <xlink_score name=\"xcorr\" value=\"%0.3f\"/>\n", dXcorr1);
   fprintf(fpxml, "      </linked_peptide>\n");
   fprintf(fpxml, "      <linked_peptide peptide=\"%s\" peptide_prev_aa=\"-\" peptide_next_aa=\"-\" protein=\"%s\" num_tot_proteins=\"1\" calc_neutral_pep_mass=\"%0.6f\" complement_mass=\"%0.6f\" designation=\"beta\">\n",
         szPep2, szProt2, dCalcMass2, dCalcMass1+xl);
   fprintf(fpxml, "       <modification_info>\n");
   for (i=0; i<strlen(szPep2); i++)
      if (szPep2[i]=='K')
         break;
   fprintf(fpxml, "        <mod_aminoacid_mass position=\"%d\" mass=\"325.127385\"/>\n", i+1);
   fprintf(fpxml, "       </modification_info>\n");
   fprintf(fpxml, "       <xlink_score name=\"score\" value=\"%0.3E\"/>\n", dExpect2);
// fprintf(fpxml, "       <xlink_score name=\"delta_score\" value=\"%0.3f\"/>\n", dXcorr2);
   fprintf(fpxml, "      </linked_peptide>\n");
   fprintf(fpxml, "     </xlink>\n");
   fprintf(fpxml, "     <search_score name=\"kojak_score\" value=\"%0.3E\"/>\n", dExpectCombined);
   fprintf(fpxml, "     <search_score name=\"delta_score\" value=\"%0.3f\"/>\n", dXcorrCombined);
   fprintf(fpxml, "     <search_score name=\"ppm_error\" value=\"0.0\"/>\n");
   fprintf(fpxml, "    </search_hit>\n");
   fprintf(fpxml, "   </search_result>\n");
   fprintf(fpxml, "  </spectrum_query>\n");
   
}

