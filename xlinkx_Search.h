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

#ifndef _XLINKXREADFASTA_H_
#define _XLINKXREADFASTA_H_

#include "Common.h"
#include "xlinkx_DataInternal.h"

#include "hash/xlinkx-hash.h"

class xlinkx_Search
{
public:
   xlinkx_Search();
   ~xlinkx_Search();

   static void SearchForPeptides(char *szMZXML,
                                 const char *,
                                 enzyme_cut_params,
                                 const char *,
                                 const char *);
   static bool WithinTolerance(double dMass1,
                               double dMass2);

private:

   static double XcorrScore(const char *szPeptide,
                     int iScanNumber);

   static bool CalculateEValue(int *hist_pep,
                               int iMatchPepCount,
                               double *dSlope,
                               double *dIntercept,
                               double dNeutralPepMass,
                               int iScanNumber);

   static void LinearRegression(int *piHistogram,
                                double *slope,
                                double *intercept,
                                int *iMaxXcorr,
                                int *iStartXcorr,
                                int *iNextXcorr);

   static bool GenerateXcorrDecoys(double dNeutralPepMass,
                                   int iMatchPepCount,
                                   int *hist_pep,
                                   int iScanNumber);

   static void WritePepXMLHeader(FILE *fpxml,
                                 char *szBaseName);

   static void WriteSpectrumQuery(FILE *fpxml,
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
                                  int iScan);

   // Private static methods

   // Private member variables
};

#endif // _XLINKXREADFASTA_H_
