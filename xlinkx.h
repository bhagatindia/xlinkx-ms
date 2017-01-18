#ifndef _XLINKX_H_
#define _XLINKX_H_

#define SIZE_BUF    8192
#define SIZE_FILE   512
#define MAX_PEAKS   1000

using namespace std;

#include "MSReader.h"
#include "Spectrum.h"
#include "MSObject.h"
#include "math.h"
#include <vector>
#include <cfloat>

using namespace MSToolkit;

struct PrecursorsStruct
{
   int    iCharge1;
   int    iCharge2;
   int    iIntensity1;
   int    iIntensity2;
   double dNeutralMass1;
   double dNeutralMass2;
};

// MS/MS scan information
struct ScanDataStruct
{
   int    iScanNumber;

// int    iCharge1;           // charge state of released peptide1
// int    iCharge2;           // charge state of released peptide2
   double dNeutralMass1;      // precursor m/z of released peptide1
   double dNeutralMass2;      // precursor m/z of released peptide2

   int    iPrecursorScanNumber;
   int    iPrecursorCharge;   // charge state of intact precursor
   double dPrecursorMZ;       // intact precursor m/z
   double dPrecursorMH;       // intact precursor mh+
   double dHardklorPrecursorNeutralMass;  // precursor m/z found from Hardklor

   vector<PrecursorsStruct> pvdPrecursors;
};

extern vector<ScanDataStruct> pvSpectrumList;


void READ_MZXMLSCANS(char *szMZXML);
void READ_HK1(char *szHK);
void READ_HK2(char *szHK);
int WITHIN_TOLERANCE(double dMass1, double dMass2);


#endif // _XLINKX_H_
