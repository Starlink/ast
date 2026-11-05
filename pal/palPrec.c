/*
*+
*  Name:
*     palPrec

*  Purpose:
*     Form the matrix of precession between two epochs (IAU 2006)

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     palPrec( double ep0, double ep1, double rmatp[3][3] )

*  Arguments:
*     ep0 = double (Given)
*        Beginning epoch
*     ep1 = double (Given)
*        Ending epoch
*     rmatp = double[3][3] (Returned)
*        Precession matrix

*  Description:
*     The IAU 2006 precession matrix from ep0 to ep1 is found and
*     returned. The matrix is in the sense  V(EP1)  =  RMATP * V(EP0).
*     The epochs are TDB (loosely TT) Julian epochs.
*
*     The matrix is formed from the IAU 2006 bias-precession matrices
*     for the two epochs, from GCRS to the mean equator and equinox of
*     each, as  RMATP = PB(EP1) * transpose( PB(EP0) ).  The frame bias
*     cancels, so the result is precession alone, and the identity
*     matrix when the two epochs are equal.
*
*     Though the matrix method itself is rigorous, the precession
*     angles are expressed through canonical polynomials which are
*     valid only for a limited time span of a few hundred years around
*     the current epoch.

*  Authors:
*     PTW: Pat Wallace (STFC)
*     DSB: David Berry (JAC, Hawaii)
*     TIMJ: Tim Jenness (Rubin Observatory)
*     {enter_new_authors_here}

*  History:
*     2012-02-10 (DSB):
*        Initial version with documentation taken from Fortran SLA
*        Adapted with permission from the Fortran SLALIB library.
*     2026-11-05 (TIMJ):
*        Always form PB(EP1) * transpose( PB(EP0) ).  When either epoch
*        was exactly 2000.0 the bias-precession matrix of the other was
*        returned on its own, so the frame bias (up to 23 mas) was
*        applied, and the matrix for 2000.0 to 2000.0 was not the
*        identity.  For other epochs the two matrices were multiplied
*        in the reverse order, giving errors of up to 6 mas between
*        1950 and 2050, and 47 mas between 1900 and 2100.
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 1996 Rutherford Appleton Laboratory
*     Copyright (C) 2012 Science and Technology Facilities Council.
*     All Rights Reserved.

*  Licence:
*     This program is free software: you can redistribute it and/or
*     modify it under the terms of the GNU Lesser General Public
*     License as published by the Free Software Foundation, either
*     version 3 of the License, or (at your option) any later
*     version.
*
*     This program is distributed in the hope that it will be useful,
*     but WITHOUT ANY WARRANTY; without even the implied warranty of
*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*     GNU Lesser General Public License for more details.
*
*     You should have received a copy of the GNU Lesser General
*     License along with this program.  If not, see
*     <http://www.gnu.org/licenses/>.

*  Bugs:
*     {note_any_bugs_here}
*-
*/

#include "pal.h"
#include "pal1sofa.h"

void palPrec( double ep0, double ep1, double rmatp[3][3] ){

/* Local Variables: */
   double rmatq[3][3];
   double ep0_days;
   double ep1_days;

/* Convert supplied dates to days since J2000 */
   ep0_days = ( ep0 - 2000.0 )*ERFA_DJY;
   ep1_days = ( ep1 - 2000.0 )*ERFA_DJY;

/* Get the bias-precession matrices from GCRS to the mean equator and
   equinox of each epoch. Each includes the frame bias, which is all of
   the matrix at J2000. */
   eraPmat06( ERFA_DJ00, ep0_days, rmatq );
   eraPmat06( ERFA_DJ00, ep1_days, rmatp );

/* Rotate from the mean of EP0 to GCRS (the transpose of the first
   matrix) and then from GCRS to the mean of EP1. The frame bias
   cancels. */
   eraTr( rmatq, rmatq );
   eraRxr( rmatp, rmatq, rmatp );

}
