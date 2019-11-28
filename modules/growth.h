///////////////////////////////////////////////////////////////////////////////////////
/// \file growth.h
/// \brief The growth module header file
///
/// Vegetation C allocation, litter production, tissue turnover
/// leaf phenology, allometry and growth.
///
/// \author Ben Smith
/// $Date: 2011-05-16 09:39:42 +0200 (Mo, 16 Mai 2011) $
///
///////////////////////////////////////////////////////////////////////////////////////

// WHAT SHOULD THIS FILE CONTAIN?
// Module header files need normally contain only declarations of functions defined in
// the module that are to be accessible to the calling framework or to other modules.

#ifndef LPJ_GUESS_GROWTH_H
#define LPJ_GUESS_GROWTH_H

#include "guess.h"

double fracmass_lpj(double fpc_low,double fpc_high,Individual& indiv);
void leaf_phenology(Patch& patch,Climate& climate);
bool allometry(Individual& indiv); // guess2008 - now returns bool instead of void
void allocation_init(double bminit,double ltor,Individual& indiv);
void growth(Stand& stand,Patch& patch);

#endif // LPJ_GUESS_GROWTH_H
