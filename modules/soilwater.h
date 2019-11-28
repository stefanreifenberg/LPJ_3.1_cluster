///////////////////////////////////////////////////////////////////////////////////////
/// \file soilwater.h
/// \brief Soil hydrology and snow
///
/// Version including evaporation from soil surface, based on work by Dieter Gerten,
/// Sibyll Schaphoff and Wolfgang Lucht, Potsdam
///
/// \author Ben Smith
/// $Date: 2011-12-08 12:57:43 +0100 (Do, 08 Dez 2011) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_SOILWATER_H
#define LPJ_GUESS_SOILWATER_H

#include "guess.h"
void initial_infiltration(Patch& patch, Climate& climate);
void soilwater(Patch& patch, Climate& climate);

#endif // LPJ_GUESS_SOILWATER_H
