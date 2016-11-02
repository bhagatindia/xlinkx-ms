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


xlinkx_Search::~xlinkx_Search()
{
}


void xlinkx_Search::SearchForPeptides(const char *protein_file, enzyme_cut_params params, const char *pep_hash_file)
{
   int i;
   int ii;

   // If PeptideHash not present, generate it now; otherwise open the hash file.
   protein_hash_db_t phdp = phd_retrieve_hash_db(protein_file, params, pep_hash_file);

   for (i=0; i<(int)pvSpectrumList.size(); i++)
   {
      for (ii=0; ii<(int)pvSpectrumList.at(i).pvdPrecursors.size(); ii++)
      {

         printf("scan %d, %f, %f\n", pvSpectrumList.at(i).iScanNumber,
               pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1,
               pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2);


         cout << "Retrieving peptides of mass " << pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1 << 
            " and " << pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2 << endl;

         vector<string*> *peptides = phdp->phd_get_peptides_ofmass(pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1);
         for (string *peptide : *peptides)
         {

            char *szPeptide = new char[(*peptide).length() + 1];
            strcpy(szPeptide, (*peptide).c_str() );

            double dXcorr = XcorrScore(szPeptide, pvSpectrumList.at(i).iScanNumber);

            cout << "pep1: " << *peptide << "  xcorr " << dXcorr << endl;
         }

         peptides = phdp->phd_get_peptides_ofmass(pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2);
         for (string *peptide : *peptides)
         {

            char *szPeptide = new char[(*peptide).length() + 1];
            strcpy(szPeptide, (*peptide).c_str() );

            double dXcorr = XcorrScore(szPeptide, pvSpectrumList.at(i).iScanNumber);

            cout << "pep2: " << *peptide << "  xcorr " << dXcorr << endl;
         }
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

      for (int i=0; i<iLenPeptide; i++) // will ignore multiple fragment ion charge states for now
      {
         dBion += g_staticParams.massUtility.pdAAMassFragment[(int)szPeptide[i]];
         if (szPeptide[i] == 'K' && dBionLysine==false)
         {
            dBion += 325.12918305;
            dBionLysine = true;
         }

         bin = BIN(dBion);
         x =  bin / SPARSE_MATRIX_SIZE;
         if (g_pvQuery.at(iWhichQuery)->ppfSparseFastXcorrData[x]==NULL || x>iMax) // x should never be > iMax so this is just a safety check
            continue;
         y = bin - (x*SPARSE_MATRIX_SIZE);
         dXcorr += g_pvQuery.at(iWhichQuery)->ppfSparseFastXcorrData[x][y];


         dYion += g_staticParams.massUtility.pdAAMassFragment[(int)szPeptide[iLenPeptide -1 - i]];
         if (szPeptide[iLenPeptide -1 - i] == 'K' && dYionLysine==false)
         {
            dYion += 325.12918305;
            dYionLysine = true;
         }

         bin = BIN(dYion);
         x =  bin / SPARSE_MATRIX_SIZE;
         if (g_pvQuery.at(iWhichQuery)->ppfSparseFastXcorrData[x]==NULL || x>iMax) // x should never be > iMax so this is just a safety check
            continue;
         y = bin - (x*SPARSE_MATRIX_SIZE);
         dXcorr += g_pvQuery.at(iWhichQuery)->ppfSparseFastXcorrData[x][y];
      }
   }

   if (bBionLysine == false || bYionLysine == false) // sanity check
   {
      cout << " Error, no internal lysine: " << szPeptide << endl;
      exit(1);
   }

   dXcorr *= 0.005;

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
