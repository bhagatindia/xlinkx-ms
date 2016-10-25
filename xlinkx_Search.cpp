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
	 for (string *peptide : *peptides) {
		 cout << *peptide << endl;
	 }

	 peptides = phdp->phd_get_peptides_ofmass(pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2);
	 for (string *peptide : *peptides) {
		 cout << *peptide << endl;
	 }

         // Grab all peptides of mass pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass1

         // Grab all peptides of mass pvSpectrumList.at(i).pvdPrecursors.at(ii).dNeutralMass2

         // Now loop through both sets of peptides, analyzing combinations of each
         // - account for the modification on internal lysing residue
         // - calculate fragments
         // - generate xcorr score
         

/*
         strcpy(szPep1, "DLRSKTDAAQISK");
         strcpy(szPep2, "IGSEKTADLK");

         XcorrScore(szPep1, szPep2, 
*/
      }
   }
}


bool xlinkx_Search::WithinTolerance(double dMass1,
                                    double dMass2)
{
   if (dMass1>0.0 &&  (1E6 * fabs(dMass1 - dMass2)/dMass1) <= g_staticParams.tolerances.dInputTolerance)
      return true;
   else
      return false;
}
