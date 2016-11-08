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

void xlinkx_Search::SearchForPeptides(const char *protein_file, enzyme_cut_params params, const char *pep_hash_file)
{
   int i;
   int ii;
   double dTolerance;
   double dPPM = 20.0;  // use 20ppm tolerance for now

#define LYSINE_MOD 197.0324

   char *toppep1[NUMPEPTIDES], *toppep2[NUMPEPTIDES];
   float xcorrPep1[NUMPEPTIDES], xcorrPep2[NUMPEPTIDES];

   // If PeptideHash not present, generate it now; otherwise open the hash file.
   protein_hash_db_t phdp = phd_retrieve_hash_db(protein_file, params, pep_hash_file);

   for (i=0; i<(int)pvSpectrumList.size(); i++)
   {
      for (ii=0; ii<(int)pvSpectrumList.at(i).pvdPrecursors.size(); ii++)
      {

         for (int li = 0; li < NUMPEPTIDES; li++) {
            xcorrPep1[li] = xcorrPep2[li] = -99999;
            toppep1[li] = toppep2[li] = NULL;
         }

         printf("scan %d, %f, %f\n", pvSpectrumList.at(i).iScanNumber,
               pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1,
               pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2);

         cout << "Retrieving peptides of mass " << pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1 << 
            " and " << pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2 << endl;

         double pep_mass1 = pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1 - LYSINE_MOD - g_staticParams.precalcMasses.dOH2;
         cout << "After Lysine residue reduction the peptide of mass " << pep_mass1 << " are being extracted" << endl;
         if (pep_mass1 <= 0) {
            cout << "Peptide mass is coming out to be zero after removing Lysine resideu" << endl;
            exit(1);
         }

         dTolerance = (dPPM * pep_mass1) / 1e6;
         vector<string*> *peptides = phdp->phd_get_peptides_ofmass(pep_mass1);
         for (string *peptide : *peptides)
         {
            char *szPeptide = new char[(*peptide).length() + 1];
            strcpy(szPeptide, (*peptide).c_str() );

            double dXcorr = XcorrScore(szPeptide, pvSpectrumList.at(i).iScanNumber);

            cout << "pep1: " << *peptide << "  xcorr " << dXcorr << endl;
            insert_pep_pq(toppep1, xcorrPep1, szPeptide, dXcorr);
         }

         cout << "Top "<< NUMPEPTIDES << " pep1 peptides for this scan are " << endl;

         for (int li = 0 ; li < NUMPEPTIDES; li++) cout << "pep1_top: " << ((toppep1[li] != NULL)? toppep1[li]: "") << " xcorr " << xcorrPep1[li] << endl;

         double pep_mass2 = pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2 - LYSINE_MOD - g_staticParams.precalcMasses.dOH2;
         cout << "After Lysine residue reduction the peptide of mass " << pep_mass2 << " are being extracted" << endl;
         if (pep_mass2 <= 0) {
            cout << "Peptide mass is coming out to be less than or equal to zero after removing lysine residue" << endl;
            exit(1);
         }
         peptides = phdp->phd_get_peptides_ofmass(pep_mass2);
         for (string *peptide : *peptides)
         {

            char *szPeptide = new char[(*peptide).length() + 1];
            strcpy(szPeptide, (*peptide).c_str() );

            double dXcorr = XcorrScore(szPeptide, pvSpectrumList.at(i).iScanNumber);

            cout << "pep2: " << *peptide << "  xcorr " << dXcorr << endl;
            insert_pep_pq(toppep2, xcorrPep2, szPeptide, dXcorr);
         }

         cout << "Top "<< NUMPEPTIDES << " pep2 peptides for this scan are " << endl;

         for (int li = 0; li < NUMPEPTIDES; li++) cout << "pep2_top: " << (toppep2[li] != NULL? toppep2[li] : "") << " xcorr " << xcorrPep2[li] << endl;
      }
   }
}


double xlinkx_Search::XcorrScore(char *szPeptide,
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
         if (g_pvQuery.at(iWhichQuery)->ppfSparseFastXcorrData[x]==NULL || x>iMax) // x should never be > iMax so this is just a safety check
            continue;
         y = bin - (x*SPARSE_MATRIX_SIZE);
         dXcorr += g_pvQuery.at(iWhichQuery)->ppfSparseFastXcorrData[x][y];


         dYion += g_staticParams.massUtility.pdAAMassFragment[(int)szPeptide[iLenPeptide -1 - i]];
         if (szPeptide[iLenPeptide -1 - i] == 'K' && !bYionLysine)
         {
            dYion += LYSINE_MOD;
            bYionLysine = true;
         }

         bin = BIN(dYion);
         x =  bin / SPARSE_MATRIX_SIZE;
         if (g_pvQuery.at(iWhichQuery)->ppfSparseFastXcorrData[x]==NULL || x>iMax) // x should never be > iMax so this is just a safety check
            continue;
         y = bin - (x*SPARSE_MATRIX_SIZE);
         dXcorr += g_pvQuery.at(iWhichQuery)->ppfSparseFastXcorrData[x][y];
      }

      if (!bBionLysine && !bYionLysine) // sanity check
      {
         cout << " Error, no internal lysine: " << szPeptide << endl;
         exit(1);
      }

      dXcorr *= 0.005;
   }

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
