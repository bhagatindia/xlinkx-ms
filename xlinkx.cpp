#include <stdio.h>
#include <stdlib.h>

#include "Common.h"
#include "xlinkx.h"
#include "xlinkx_DataInternal.h"
#include "xlinkx_Preprocess.h"
#include "xlinkx_Search.h"
#include "xlinkx_MassSpecUtils.h"


// This program takes in an mzXML as input.
// If corresponding .hk file does not exist, it runs Hardklor to create it.
// Then reading each MS/MS scan in mzXML, find 2 peaks in .hk file
// that add up to ReACT relationship; need mzXML for precursor mass.

vector<ScanDataStruct> pvSpectrumList;

int main(int argc, char **argv)
{
   char szHK1[SIZE_FILE];
   char szHK2[SIZE_FILE];
   char szMZXML[SIZE_FILE];

   if (argc != 4)
   {
      printf("  USAGE:  xlinkx input.mzXML protein_database.txt peptide_hash_database.txt .\n\n");
      printf("  Need to enter mzXML file on command line.\n\n");
      exit(1);
   }

   // First read through mzXML file, grab all MS/MS scans and
   // their precursor m/z.  Require precursor charge be 4+
   // or larger?
   strcpy(szMZXML, argv[1]);

   strcpy(szHK1, szMZXML);
   szHK1[strlen(szHK1)-5]='\0';
   strcpy(szHK2, szHK1);
   strcat(szHK1, "hk1");        // ms1 hardklor run
   strcat(szHK2, "hk2");        // ms2 hardklor run


   // This first pass read simply gets all ms/ms scans and their measured precursor m/z
   READ_MZXMLSCANS(szMZXML);

   // Next, go to Hardklor .hk1 file to get accurate precursor m/z
   READ_HK1(szHK1);

   // Now, read through .hk2 file to find accurate peptide masses that add up to precursor
   READ_HK2(szHK2);

   // Load and preprocess all MS/MS scans that have a pair of peptide masses that add up to precursor
   g_staticParams.tolerances.dFragmentBinSize = 1.0005;
   g_staticParams.dInverseBinWidth = 1.0 /g_staticParams.tolerances.dFragmentBinSize;
   LOAD_SPECTRA(szMZXML);

   int iCount=0;
   for (int ii=0; ii<(int)pvSpectrumList.size(); ii++)
   {
      if (pvSpectrumList.at(ii).pvdPrecursors.size() > 0)
         iCount++;
   }
   printf(" #spectra with relationship: %d    #total spectra:  %d\n", iCount, (int)pvSpectrumList.size());


   // Apply some settings that might be better applied somewhere else
   xlinkx_MassSpecUtils::AssignMass(g_staticParams.massUtility.pdAAMassParent,
                                    g_staticParams.massUtility.bMonoMassesParent,
                                    &g_staticParams.massUtility.dOH2parent);

   xlinkx_MassSpecUtils::AssignMass(g_staticParams.massUtility.pdAAMassFragment,
                                    g_staticParams.massUtility.bMonoMassesFragment,
                                    &g_staticParams.massUtility.dOH2fragment);

   g_staticParams.precalcMasses.iMinus17 = BIN(g_staticParams.massUtility.dH2O);
   g_staticParams.precalcMasses.iMinus18 = BIN(g_staticParams.massUtility.dNH3);

   g_staticParams.precalcMasses.dNtermProton = g_staticParams.staticModifications.dAddNterminusPeptide
      + PROTON_MASS;

   g_staticParams.precalcMasses.dCtermOH2Proton = g_staticParams.staticModifications.dAddCterminusPeptide
      + g_staticParams.massUtility.dOH2fragment
      + PROTON_MASS;

   g_staticParams.precalcMasses.dOH2ProtonCtermNterm = g_staticParams.massUtility.dOH2parent
      + PROTON_MASS
      + g_staticParams.staticModifications.dAddCterminusPeptide
      + g_staticParams.staticModifications.dAddNterminusPeptide;

   strcpy(g_staticParams.enzymeInformation.szSearchEnzymeName, "trypsin");
   strcpy(g_staticParams.enzymeInformation.szSearchEnzymeBreakAA, "KR");
   strcpy(g_staticParams.enzymeInformation.szSearchEnzymeNoBreakAA, "P");
   g_staticParams.enzymeInformation.iSearchEnzymeOffSet = 1;

   g_staticParams.options.iEnzymeTermini = ENZYME_DOUBLE_TERMINI;
   g_staticParams.options.bNoEnzymeSelected = false;

   enzyme_cut_params params;
   params.semi_tryptic = 0;
   params.precut_amino = "-";
   params.prenocut_amino = "-";
   params.missed_cleavage = 0;
   params.postcut_amino = "KR";
   params.postnocut_amino = "P";

   // Now open fasta file and get a list of all peptides with masses close to 
   xlinkx_Search::SearchForPeptides(argv[2], params, argv[3]);

   return 0;

} //main


void READ_MZXMLSCANS(char *szMZXML)
{
   MSReader mstReader;
   Spectrum mstSpectrum;
   int iFileLastScan;
   int iScanNumber=1;
   int iMS1ScanNumber=0;

   printf(" reading %s ... ", szMZXML); fflush(stdout);

   // We want to read only MS1/MS2 scans.
   vector<MSSpectrumType> msLevel;
   msLevel.push_back(MS1);
   msLevel.push_back(MS2);

   mstReader.setFilter(msLevel);
   mstReader.readFile(szMZXML, mstSpectrum, 1);
   iFileLastScan = mstReader.getLastScan();

   for (iScanNumber = 1 ; iScanNumber <= iFileLastScan; iScanNumber++)
   {
      // Loads in MSMS spectrum data.
      mstReader.readFile(NULL, mstSpectrum, iScanNumber);

      if (mstSpectrum.getMsLevel() == 1)
      {
         iMS1ScanNumber = iScanNumber;
      }
      else if (mstSpectrum.getMsLevel() == 2)
      {
         struct ScanDataStruct pData;
        
         if (mstSpectrum.sizeZ() > 0)
            pData.iPrecursorCharge = mstSpectrum.atZ(0).z;
         else
            pData.iPrecursorCharge = 0;

         pData.dPrecursorMZ = mstSpectrum.getMZ();
         pData.iScanNumber = iScanNumber;
         pData.iPrecursorScanNumber = iMS1ScanNumber;
         pData.dHardklorPrecursorNeutralMass = 0.0;
         pData.dNeutralMass1 = 0.0;
         pData.dNeutralMass2 = 0.0;

         pvSpectrumList.push_back(pData);
      }
      else
      {
         printf(" Error ... should never get here (msLevel).\n");
         exit(1);
      }

      if (!(iScanNumber%200))
      {
         printf("%3d%%", (int)(100.0*iScanNumber/iFileLastScan));
         fflush(stdout);
         printf("\b\b\b\b");
      }
   }

   printf("100%%\n");
}


void READ_HK1(char *szHK)
{
   FILE *fp;
   char szBuf[SIZE_BUF];
   int iListCt;       // this will keep an index of pvSpectrumList
   long lFP = 0;
   long lEndFP;

   printf(" reading %s ... ", szHK); fflush(stdout);

   if ( (fp=fopen(szHK, "r"))== NULL)
   {
      printf("Please generate .hk1 file with same base name as .mzXML input.\n");
      exit(1);

      // Instead of reporting error above, generate HK1 file here.
      // This requires creating a Hardklor params file and executing hardklor.

      // run Hardklor

      // now try to re-open .hk file
      if ( (fp=fopen(szHK, "r"))== NULL)
      {
         printf(" Error ... cannot create/read %s file.\n", szHK);
         exit(1);
      }
   }

   fseek(fp, 0, SEEK_END);
   lEndFP = ftell(fp);  // store current
   rewind(fp);

   iListCt = 0;
   while (fgets(szBuf, SIZE_BUF, fp))
   {
      if (szBuf[0]=='S')
      {
         int iScanNumber;

         sscanf(szBuf, "S\t%d\t", &iScanNumber);

         while (pvSpectrumList.at(iListCt).iPrecursorScanNumber < iScanNumber)
            iListCt++;

         if (pvSpectrumList.at(iListCt).iPrecursorScanNumber == iScanNumber)
         {
            int iMaxIntensity = 0;

            // Get accurate precursor m/z from deconvoluted MS1 peaks
            lFP = ftell(fp);  // store current

            while (fgets(szBuf, SIZE_BUF, fp))
            {
               if (szBuf[0] == 'S')
               {
                  fseek(fp, lFP, SEEK_SET);
                  break;
               }

               lFP = ftell(fp);

               if (szBuf[0] == 'P')
               {
                  double dMass;
                  double dMZ;
                  int iCharge;
                  int iIntensity;

                  sscanf(szBuf, "P\t%lf\t%d\t%d\t", &dMass, &iCharge, &iIntensity);

                  dMZ = (dMass + (iCharge*PROTON))/iCharge;

                  if (fabs(dMZ -  pvSpectrumList.at(iListCt).dPrecursorMZ) < 0.5 && iIntensity > iMaxIntensity)
                  {
                     iMaxIntensity = iIntensity;
                     pvSpectrumList.at(iListCt).dHardklorPrecursorNeutralMass = dMass;
                  }
               }
            }
         }
      
         if (iScanNumber%500)
         {
            printf("%3d%%", (int)(100.0*lFP/lEndFP));
            fflush(stdout);
            printf("\b\b\b\b");
         }

      }
   }

   printf("100%%\n");
   fclose(fp);
}


void READ_HK2(char *szHK)
{
   FILE *fp;
   char szBuf[SIZE_BUF];
   int iListCt;       // this will keep an index of pvSpectrumList
   long lFP = 0;
   long lEndFP;

   printf(" reading %s ... ", szHK); fflush(stdout);

   if ( (fp=fopen(szHK, "r"))== NULL)
   {
      printf("Please generate .hk2 file with same base name as .mzXML input.\n");
      exit(1);

      // generate HK file here

      // run Hardklor

      // now try to re-open .hk file
      if ( (fp=fopen(szHK, "r"))== NULL)
      {
         printf(" Error ... cannot create/read %s file.\n", szHK);
         exit(1);
      }
   }

   fseek(fp, 0, SEEK_END);
   lEndFP = ftell(fp);  // store current
   rewind(fp);

   iListCt = 0;
   while (fgets(szBuf, SIZE_BUF, fp))
   {
      if (szBuf[0]=='S')
      {
         int iScanNumber;

         sscanf(szBuf, "S\t%d\t", &iScanNumber);

         while (pvSpectrumList.at(iListCt).iScanNumber < iScanNumber)
            iListCt++;

         if (pvSpectrumList.at(iListCt).iScanNumber == iScanNumber)
         {
            int iNumPeaks=0;

            struct
            {
               double dNeutralMass;
               int iCharge;
               int iIntensity;
            } pPeaks[MAX_PEAKS];

            lFP = ftell(fp);  // store current

            // Read in all deconvoluted peaks in this MS/MS scan
            while (fgets(szBuf, SIZE_BUF, fp))
            {
               if (szBuf[0] == 'S')
               {
                  fseek(fp, lFP, SEEK_SET);
                  break;
               }

               lFP = ftell(fp);

               if (szBuf[0] == 'P')
               {
                  double dMass;
                  int iCharge;
                  int iIntensity;

                  sscanf(szBuf, "P\t%lf\t%d\t%d\t", &dMass, &iCharge, &iIntensity);

                  if (iNumPeaks == MAX_PEAKS)
                  {
                     printf(" Error, scan %d has MAX_PEAKS (%d) number of deconvoluted peaks.\n", iScanNumber, MAX_PEAKS);
                     printf(" Need to increase definition of MAX_PEAKS\n");
                     exit(1);
                  }

                  pPeaks[iNumPeaks].dNeutralMass = dMass;
                  pPeaks[iNumPeaks].iCharge = iCharge;

                  iNumPeaks++;
               }
            }

            // At this point, we know precursor m/z and precursor charge and have read all ms/ms peaks.
            // Must find 2 peptides that add up to intact cross-link (ms1+ms2+reporter)
            // where the charge states of the peptides can't be larger than precursor charge.
            int i;
            int ii;
            double dReporter = 751.406065;   // need to make this a parameter entry

            for (i=0; i<iNumPeaks; i++)
            {
               for (ii=i; ii<iNumPeaks; ii++)
               {
                  // Placing check here that the peptide masses must be greater than some minimum
                  if (pPeaks[i].dNeutralMass >= 500.0 && pPeaks[ii].dNeutralMass >= 500.0)
                  {
                     double dCombinedMass = pPeaks[i].dNeutralMass + pPeaks[ii].dNeutralMass + dReporter;

                     if (WITHIN_TOLERANCE(dCombinedMass, pvSpectrumList.at(iListCt).dHardklorPrecursorNeutralMass))
                     {
                        if (pPeaks[i].iCharge + pPeaks[ii].iCharge <= pvSpectrumList.at(iListCt).iPrecursorCharge
                              && pPeaks[i].iCharge < pvSpectrumList.at(iListCt).iPrecursorCharge
                              && pPeaks[ii].iCharge < pvSpectrumList.at(iListCt).iPrecursorCharge)
                        {
                           struct PrecursorsStruct pTmp;

                           if (pPeaks[i].dNeutralMass < pPeaks[ii].dNeutralMass)
                           {
                              pTmp.dNeutralMass1 = pPeaks[i].dNeutralMass;
                              pTmp.dNeutralMass2 = pPeaks[ii].dNeutralMass;
                              pTmp.iCharge1 = pPeaks[i].iCharge;
                              pTmp.iCharge2 = pPeaks[i].iCharge;
                           }
                           else
                           {
                              pTmp.dNeutralMass2 = pPeaks[i].dNeutralMass;
                              pTmp.dNeutralMass1 = pPeaks[ii].dNeutralMass;
                              pTmp.iCharge2 = pPeaks[i].iCharge;
                              pTmp.iCharge1 = pPeaks[i].iCharge;
                           }

                           // Do a quick check here and push_back only if the two
                           // masses are not very similar to existing masses
                           
                           bool bMassesAlreadyPresent = false;

                           for (int iii=0; iii<(int)pvSpectrumList.at(iListCt).pvdPrecursors.size(); iii++)
                           {
                              if (WITHIN_TOLERANCE(pTmp.dNeutralMass1, pvSpectrumList.at(iListCt).pvdPrecursors.at(iii).dNeutralMass1)
                                    && WITHIN_TOLERANCE(pTmp.dNeutralMass2, pvSpectrumList.at(iListCt).pvdPrecursors.at(iii).dNeutralMass2))
                              {
                                 bMassesAlreadyPresent = true;
                                 break;
                              }
                           }

                           if (!bMassesAlreadyPresent)
                              pvSpectrumList.at(iListCt).pvdPrecursors.push_back(pTmp);
                        }
                     }
                  }
               }
            }
         }
      
         printf("%3d%%", (int)(100.0*lFP/lEndFP));
         fflush(stdout);
         printf("\b\b\b\b");

      }
   }

   printf("100%%\n");
   fclose(fp);
}


int WITHIN_TOLERANCE(double dMass1, double dMass2)
{
   double dPPM = 20.0;

   if (1E6 * fabs(dMass1 - dMass2)/dMass2 <= dPPM)
      return 1;
   else
      return 0;
}


void LOAD_SPECTRA(char *szMZXML)
{
   MSReader mstReader;
   Spectrum mstSpectrum;

   printf(" reading %s ... ", szMZXML); fflush(stdout);

   // We want to read only MS2 scans.
   vector<MSSpectrumType> msLevel;
   msLevel.push_back(MS2);

   mstReader.setFilter(msLevel);
   mstReader.readFile(szMZXML, mstSpectrum, 1);

   xlinkx_preprocess::AllocateMemory(1);

   int iReadScans = 0;
   int iSkippedScans = 0;
   int iLast = (int)pvSpectrumList.size();
   for (int i=0; i<iLast; i++)
   {
      if ((int)pvSpectrumList.at(i).pvdPrecursors.size() > 0) // this means pair of peptide masses found
      {
         xlinkx_preprocess::LoadAndPreprocessSpectra(mstReader, pvSpectrumList.at(i).iScanNumber);
         iReadScans++;
      }
      else
         iSkippedScans++;

      if (i%200)
      {  
         printf("%3d%%", (int)(100.0*i/iLast));
         fflush(stdout);
         printf("\b\b\b\b");
      }
   }

   printf("100%%\n");

   //printf("iRead %d, iSkipped %d\n\n", iReadScans, iSkippedScans);
}
