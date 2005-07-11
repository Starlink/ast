#define astCLASS
#include "axis.h"
#include "box.h"
#include "channel.h"
#include "circle.h"
#include "cmpframe.h"
#include "cmpmap.h"
#include "cmpregion.h"
#include "dsbspecframe.h"
#include "dssmap.h"
#include "ellipse.h"
#include "fitschan.h"
#include "fluxframe.h"
#include "timeframe.h"
#include "timemap.h"
#include "frame.h"
#include "frameset.h"
#include "grismmap.h"
#include "interval.h"
#include "intramap.h"
#include "keymap.h"
#include "loader.h"
#include "lutmap.h"
#include "mapping.h"
#include "mathmap.h"
#include "matrixmap.h"
#include "nullregion.h"
#include "object.h"
#include "pcdmap.h"
#include "permmap.h"
#include "plot.h"
#include "pointlist.h"
#include "pointset.h"
#include "polygon.h"
#include "polymap.h"
#include "prism.h"
#include "normmap.h"
#include "ratemap.h"
#include "region.h"
#include "shiftmap.h"
#include "skyaxis.h"
#include "skyframe.h"
#include "slamap.h"
#include "specfluxframe.h"
#include "specframe.h"
#include "specmap.h"
#include "sphmap.h"
#include "tranmap.h"
#include "unitmap.h"
#include "wcsmap.h"
#include "winmap.h"
#include "xmlchan.h"
#include "zoommap.h"
#include "stc.h"
#include "stcresourceprofile.h"
#include "stcsearchlocation.h"
#include "stccatalogentrylocation.h"
#include "stcobsdatalocation.h"

#include "error.h"
#include "ast_err.h"
#include <stddef.h>
#include <string.h>

/*
*+
*  Copyright:
*     <COPYRIGHT_STATEMENT>

*  Authors:
*     RFWS: R.F. Warren-Smith (Starlink)
*     DSB: David S. Berry (Starlink)

*  History:
*     18-NOV-1997 (RFWS):
*        Original version.
*     18-MAR-1998 (RFWS):
*        Added the IntraMap class.
*     3-JUN-1999 (RFWS):
*        Added the PcdMap class.
*     17-AUG-1999 (RFWS):
*        Added the MathMap class.
*     8-JAN-2003 (DSB):
*        Added the SpecMap and SpecFrame classes.
*     15-JUL-2003 (DSB):
*        Added the GrsimMap class.
*-
*/

AstLoaderType *astGetLoader( const char *class ) {
   if ( !astOK ) return NULL;

#define LOAD(name) \
if ( !strcmp( class, #name ) ) return (AstLoaderType *) astLoad##name##_

   LOAD(Axis);
   LOAD(Box);
   LOAD(Channel);
   LOAD(Circle);
   LOAD(CmpFrame);
   LOAD(CmpMap);
   LOAD(CmpRegion);
   LOAD(DSBSpecFrame);
   LOAD(DssMap);
   LOAD(Ellipse);
   LOAD(FitsChan);
   LOAD(FluxFrame);
   LOAD(Frame);
   LOAD(FrameSet);
   LOAD(GrismMap);
   LOAD(Interval);
   LOAD(IntraMap);
   LOAD(KeyMap);
   LOAD(LutMap);
   LOAD(Mapping);
   LOAD(MathMap);
   LOAD(MatrixMap);
   LOAD(NullRegion);
   LOAD(Object);
   LOAD(PcdMap);
   LOAD(PermMap);
   LOAD(Plot);
   LOAD(PointList);
   LOAD(PointSet);
   LOAD(PolyMap);
   LOAD(Polygon);
   LOAD(Prism);
   LOAD(NormMap); 
   LOAD(RateMap); 
   LOAD(Region);
   LOAD(ShiftMap);
   LOAD(SkyAxis);
   LOAD(SkyFrame);
   LOAD(SlaMap);
   LOAD(SpecFluxFrame);
   LOAD(SpecFrame);
   LOAD(SpecMap);
   LOAD(SphMap);
   LOAD(TimeFrame);
   LOAD(TimeMap);
   LOAD(TranMap);
   LOAD(UnitMap);
   LOAD(WcsMap);
   LOAD(WinMap);
   LOAD(XmlChan);
   LOAD(ZoomMap);

   LOAD(Stc);
   LOAD(StcResourceProfile);
   LOAD(StcSearchLocation);
   LOAD(StcCatalogEntryLocation);
   LOAD(StcObsDataLocation);

   astError( AST__OCLUK, "astGetLoader: Object of unknown class \"%s\" cannot "
                         "be loaded.", class );
   return NULL;
#undef LOAD
}
