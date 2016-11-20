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
         if (!strcmp(pepArray[i], ins_pep)) return;
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
         xcorrArray[i] = xcorrArray[i-1]; pepArray[i] = pepArray[i-1];
         xcorrArray[i-1] = temp; pepArray[i-1] = temp_pep;
      } else break;
   }
}

#define BIN_SIZE 0.1
#define MAX_XCORR_VALUE 20
#define NUM_BINS (int)(MAX_XCORR_VALUE/BIN_SIZE + 1)

inline int xlinkx_get_histogram_bin_num(float value)
{
    if (value > MAX_XCORR_VALUE) value = MAX_XCORR_VALUE;
    return value/BIN_SIZE;
}

void xlinkx_print_histogram(int hist_pep[])
{
    for (int i = 0; i < NUM_BINS; i++) printf("%d ", hist_pep[i]);
    printf("\n");
}

void xlinkx_Search::SearchForPeptides(char *szMZXML,
                                      const char *protein_file,
                                      enzyme_cut_params params,
                                      const char *pep_hash_file)
{
   int i;
   int ii;
   double dTolerance;
   double dPPM = 20.0;  // use 20ppm tolerance for now
   int hist_pep1[NUM_BINS], hist_pep2[NUM_BINS], hist_combined[NUM_BINS], num_pep1, num_pep2;

#define LYSINE_MOD 197.032422

   char *toppep1[NUMPEPTIDES], *toppep2[NUMPEPTIDES];
   float xcorrPep1[NUMPEPTIDES], xcorrPep2[NUMPEPTIDES];

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

   for (i=0; i<(int)pvSpectrumList.size(); i++)
   {
//if (pvSpectrumList.at(i).iScanNumber == 24686)
      for (ii=0; ii<(int)pvSpectrumList.at(i).pvdPrecursors.size(); ii++)
      {

         for (int i = 0; i < NUM_BINS; i++) hist_pep1[i] = hist_pep2[i] = hist_combined[i] = 0;
         num_pep1 = num_pep2 = 0;
         // JKE: need to move load and preprocess spectra here
         // Need charge states for mass1 & mass2 to exclude peaks from spectra
         // before processing.
         double dMZ1 =  (pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1
               + pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge1 * PROTON_MASS)/ pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge1;
         double dMZ2 =  (pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2
               + pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge2 * PROTON_MASS)/ pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge2;

         xlinkx_preprocess::LoadAndPreprocessSpectra(mstReader, pvSpectrumList.at(i).iScanNumber, dMZ1, dMZ2);
 
         for (int li = 0; li < NUMPEPTIDES; li++) {
            xcorrPep1[li] = xcorrPep2[li] = -99999;
            toppep1[li] = toppep2[li] = NULL;
         }

         printf("Scan %d, retrieving peptides of mass %0.4f (%d+ %0.4f) and %0.4f (%d+ %0.4f)\n",
               pvSpectrumList.at(i).iScanNumber,
               pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1,
               pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge1,
               dMZ1,
               pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2,
               pvSpectrumList.at(i).pvdPrecursors.at(ii).iCharge2,
               dMZ2);

         double pep_mass1 = pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1 - LYSINE_MOD - g_staticParams.precalcMasses.dOH2;
         cout << "After Lysine residue reduction the peptide of mass " << pep_mass1 << " are being extracted";
         cout << " (" << pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1 << ")" << endl;
         if (pep_mass1 <= 0)
         {
            cout << "Peptide mass is coming out to be zero after removing Lysine resideu" << endl;
            exit(1);
         }

         double dXcorr = 0.0;
         double dXcorr1 = 0.0;
         double dXcorr2 = 0.0;

         dTolerance = (dPPM * pep_mass1) / 1e6;
         //vector<string*> *peptides = phdp->phd_get_peptides_ofmass(pep_mass1);
         vector<string*> *peptides1 = phdp->phd_get_peptides_ofmass_tolerance(pep_mass1, 1);
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

            int bin_num = xlinkx_get_histogram_bin_num(dXcorr);
            hist_pep1[bin_num]++;
            insert_pep_pq(toppep1, xcorrPep1, szPeptide, dXcorr);
            num_pep1++;
         }

         xlinkx_print_histogram(hist_pep1);
         cout << "Top "<< NUMPEPTIDES << " pep1 peptides for this scan are " << endl;

         for (int li = 0 ; li < NUMPEPTIDES; li++) cout << "pep1_top: " << ((toppep1[li] != NULL)? toppep1[li]: "") << " xcorr " << xcorrPep1[li] << endl;

         double pep_mass2 = pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2 - LYSINE_MOD - g_staticParams.precalcMasses.dOH2;
         cout << "After Lysine residue reduction the peptide of mass " << pep_mass2 << " are being extracted";
         cout << " (" << pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2 << ")" << endl;
         if (pep_mass2 <= 0)
         {
            cout << "Peptide mass is coming out to be less than or equal to zero after removing lysine residue" << endl;
            exit(1);
         }

//pep_mass2 -= 1;  // account for weird 58 mod on Cysteine for this data

         dTolerance = (dPPM * pep_mass2) / 1e6;
         //peptides = phdp->phd_get_peptides_ofmass(pep_mass2);
         vector<string*> *peptides2 = phdp->phd_get_peptides_ofmass_tolerance(pep_mass2, 1);
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

            hist_pep2[xlinkx_get_histogram_bin_num(dXcorr)]++;
//          cout << "pep2: " << *peptide << "  xcorr " << dXcorr << endl;
            insert_pep_pq(toppep2, xcorrPep2, szPeptide, dXcorr);
            num_pep2++;
         }

         xlinkx_print_histogram(hist_pep2);
         cout << "Top "<< NUMPEPTIDES << " pep2 peptides for this scan are " << endl;

         for (int li = 0; li < NUMPEPTIDES; li++)
            cout << "pep2_top: " << (toppep2[li] != NULL? toppep2[li] : "") << " xcorr " << xcorrPep2[li] << endl;

         cout << "Size of peptide1 list is " << num_pep1 << " and the size of peptide2 list is " << num_pep2 << endl;
         // Computing the combined histogram of xcorr   
         for (string *peptide1 : *peptides1)
         {

            // sanity check to ignore peptides w/unknown AA residues
            // should not be needed now that this is addressed in the hash building
            
               dXcorr1 = XcorrScore(peptide1->c_str(), pvSpectrumList.at(i).iScanNumber);

               for (string *peptide2 : *peptides2)
               {
                   dXcorr2 = XcorrScore(peptide2->c_str(), pvSpectrumList.at(i).iScanNumber);

                   dXcorr = dXcorr1 + dXcorr2;

                   int bin_num = xlinkx_get_histogram_bin_num(dXcorr);
                   hist_combined[bin_num]++;
               }
         }

         xlinkx_print_histogram(hist_combined);
      }
   }
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

   if (iWhichQuery < (int)g_pvQuery.size())
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
