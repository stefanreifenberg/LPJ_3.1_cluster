///////////////////////////////////////////////////////////////////////////////////////
/// \file landcover.h
/// \brief Functions handling landcover aspects, such as creating or resizing Stands
///
/// $Date: 2013-07-17 09:22:52 +0200 (Mi, 17 Jul 2013) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_LANDCOVER_H
#define LPJ_GUESS_LANDCOVER_H

#include "guess.h"
#include "inputmodule.h"

///	Creates stands for landcovers present in the gridcell
void landcover_init(Gridcell& gridcell, InputModule* input_module);

/// Handles changes in the landcover fractions from year to year
/** This function will for instance kill or create new stands
 *  if needed.
 */
void landcover_dynamics(Gridcell& gridcell, InputModule* input_module);

#endif // LPJ_GUESS_LANDCOVER_H
