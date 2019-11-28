///////////////////////////////////////////////////////////////////////////////////////
/// \file parameters.h
/// \brief The parameters module is responsible for reading in the instruction file
///
/// This module defines and makes available a lot of the instruction file parameters
/// used by the model, but also lets other modules define their own parameters or
/// access "custom" parameters without defining them.
///
/// A new parameter can be added by creating a new global variable here (or a new
/// Pft member variable if it's a PFT parameter), and then declaring it in
/// plib_declarations in parameters.cpp. See the many existing examples, and
/// documentation in the PLIB library for further documentation about this.
///
/// Sometimes, adding a new parameter shouldn't (or can't) be done here however.
/// A parameter specific for a certain input module, should only be declared if
/// that input module is used. In this case the input module should declare its
/// own parameters when it is created. This can also be a good idea simply to 
/// make modules more independent. For parameters like this, we can either use
/// the "custom" parameters (\see Paramlist) which don't need to be declared at
/// all, or the parameters can be declared with the declare_parameter family of
/// functions.
///
/// \author Joe Siltberg
/// $Date: 2015-04-08 13:56:22 +0200 (Mi, 08 Apr 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_PARAMETERS_H
#define LPJ_GUESS_PARAMETERS_H

#include "gutil.h"
#include <string>


///////////////////////////////////////////////////////////////////////////////////////
// Enums needed by some of the global instruction file parameters defined below


/// Vegetation 'mode', i.e. what each Individual object represents
/** Can be one of: 
 *  1. The average characteristics of all individuals comprising a PFT
 *     population over the modelled area (standard LPJ mode)
 *  2. A cohort of individuals of a PFT that are roughly the same age
 *  3. An individual plant
 */
typedef enum {NOVEGMODE, INDIVIDUAL, COHORT, POPULATION} vegmodetype;

/// Land cover type of a stand. NLANDCOVERTYPES keeps count of number of items.
typedef enum {URBAN, CROPLAND, PASTURE, FOREST, NATURAL, PEATLAND, NLANDCOVERTYPES} landcovertype;

/// Water uptake parameterisations
/** \see water_uptake in canexch.cpp
  */
typedef enum {WR_WCONT, WR_ROOTDIST, WR_SMART, WR_SPECIESSPECIFIC} wateruptaketype;


///////////////////////////////////////////////////////////////////////////////////////
// Global instruction file parameters

/// Title for this run
extern xtring title;

/// Vegetation mode (population, cohort or individual)
extern vegmodetype vegmode;

/// Number of patches in each stand
/** Should always be 1 in population mode,
 *  cropland stands always have 1 patch.
 */
extern int npatch;

/// Patch area (m2) (individual and cohort mode only)
extern double patcharea;

// SR_20161123
extern int plant_year;

// f_js_20170126
extern int ndiamclass;
extern double diamclasswidth;

/// Whether background establishment enabled (individual, cohort mode)
extern bool ifbgestab;

/// Whether spatial mass effect enabled for establishment (individual, cohort mode)
extern bool ifsme;

/// Whether establishment stochastic (individual, cohort mode)
extern bool ifstochestab;

/// Whether mortality stochastic (individual, cohort mode)
extern bool ifstochmort;

/// Whether fire enabled
extern bool iffire;

/// Whether "generic" patch-destroying disturbance enabled (individual, cohort mode)
extern bool ifdisturb;

/// Generic patch-destroying disturbance interval (individual, cohort mode)
extern double distinterval;

/// Whether SLA calculated from leaf longevity (alt: prescribed)
extern bool ifcalcsla;

/// Whether leaf C:N ratio minimum calculated from leaf longevity (alt: prescribed)
extern bool ifcalccton;

/// Establishment interval in cohort mode (years)
extern int estinterval;

/// Whether C debt (storage between years) permitted
extern bool ifcdebt;

/// Water uptake parameterisation
extern wateruptaketype wateruptake;

/// whether CENTURY SOM dynamics (otherwise uses standard LPJ formalism)
extern bool ifcentury;
/// whether plant growth limited by available N	
extern bool ifnlim;
/// number of years to allow spinup without nitrogen limitation	
extern int freenyears;
/// fraction of nitrogen relocated by plants from roots and leaves	
extern double nrelocfrac;
/// first term in nitrogen fixation eqn (Cleveland et al 1999)	
extern double nfix_a;
/// second term in nitrogen fixation eqn (Cleveland et al 1999)	
extern double nfix_b;

/// Whether other landcovers than natural vegetation are simulated.
extern bool run_landcover;

/// Whether a specific landcover type is simulated (URBAN, CROPLAND, PASTURE, FOREST, NATURAL, PEATLAND).
extern bool run[NLANDCOVERTYPES];

/// Whether landcover fractions are read from ins-file.
extern bool lcfrac_fixed;

/// Set to false by initio( ) if fraction input files have yearly data.
extern bool all_fracs_const;

/// If a slow harvested product pool is included in patchpft.
extern bool ifslowharvestpool;

/// number of spinup years
extern int nyear_spinup;

/// silt/sand fractions per soiltype
extern bool textured_soil;



///////////////////////////////////////////////////////////////////////////////////////
// Settings controlling the saving and loading from state files

/// Location of state files
extern xtring state_path;

/// Whether to restart from state files
extern bool restart;

/// Whether to save state files
extern bool save_state;

/// Save/restart year
extern int state_year;

/// whether to vary mort_greff smoothly with growth efficiency (1) or to use the standard step-function (0)
extern bool ifsmoothgreffmort;

/// whether establishment is limited by growing season drought
extern bool ifdroughtlimitedestab;

/// rain on wet days only (1, true), or a little every day (0, false);
extern bool ifrainonwetdaysonly;

/// whether BVOC calculations are included
extern bool ifbvoc;




///////////////////////////////////////////////////////////////////////////////////////
// The Paramlist class (and Paramtype) 
//

/// Represents one custom "param" item
/** \see Paramlist */
struct Paramtype {
	xtring name;
	xtring str;
	double num;
};

/// List for the "custom" parameters
/** Functionality for storing and retrieving custom "param" items from the instruction
 *  script. "Custom" parameters can be accessed by other modules without the need to 
 *  define them beforehand. This of course also means there is no help text associated
 *  with these parameters, so the user can't get any documentation about them from
 *  the command line.
 *
 * Custom keywords may be included in the instruction script using syntax similar to
 * the following examples:
 *
 * \code
 *     param "co2" (num 340)
 *     param "file_gridlist" (str "gridlist.txt")
 * \endcode
 *
 * To retrieve the values associated with the "param" strings in the above examples,
 * use the following function calls (may appear anywhere in this file; instruction
 * script must have been read in first):
 *
 * \code
 *     param["co2"].num
 *     param["file_gridlist"].str
 * \endcode
 *
 * Each "param" item can store EITHER a number (int or double) OR a string, but not
 * both types of data. Function fail is called to terminate output if a "param" item
 * with the specified identifier was not read in.
 */
class Paramlist : public ListArray<Paramtype> {

public:
	/// Adds a parameter with a numeric value, overwriting if it already existed
	void addparam(xtring name,xtring value);

	/// Adds a parameter with a string value, overwriting if it already existed
	void addparam(xtring name,double value);

	/// Fetches a parameter from the list, aborts the program if it didn't exist
	Paramtype& operator[](xtring name);

private:
	/// Tries to find the parameter in the list
	/** \returns 0 if it wasn't there. */
	Paramtype* find(xtring name);
};

/// The global Paramlist object
/** Contains all the custom parameters after reading in the instruction file */
extern Paramlist param;

/// Reads in the instruction file
/** Uses PLIB library functions to read instructions from file specified by
 * 'insfilename'.
 */
void read_instruction_file(const char* insfilename);

/// Displays documentation about the instruction file parameters to the user
void printhelp();


///////////////////////////////////////////////////////////////////////////////////////
// Interface for declaring parameters from other modules

/// Declares an xtring parameter
/** \param name     The name of the parameter
 *  \param param    Pointer to variable where the value of the parameter is to be placed
 *  \param maxlen   Maximum allowed length of the parameter in the ins file
 *  \param help     Documentation describing the parameter to the user
 */
void declare_parameter(const char* name, xtring* param, int maxlen, const char* help = "");

/// Declares a std:string parameter
/** \param name     The name of the parameter
 *  \param param    Pointer to variable where the value of the parameter is to be placed
 *  \param maxlen   Maximum allowed length of the parameter in the ins file
 *  \param help     Documentation describing the parameter to the user
 */
void declare_parameter(const char* name, std::string* param, int maxlen, const char* help = "");

/// Declares an int parameter
/** \param name     The name of the parameter
 *  \param param    Pointer to variable where the value of the parameter is to be placed
 *  \param min      Minimum allowed value of the parameter in the ins file
 *  \param max      Maximum allowed value of the parameter in the ins file
 *  \param help     Documentation describing the parameter to the user
 */
void declare_parameter(const char* name, int* param, int min, int max, const char* help = "");

/// Declares a double parameter
/** \param name     The name of the parameter
 *  \param param    Pointer to variable where the value of the parameter is to be placed
 *  \param min      Minimum allowed value of the parameter in the ins file
 *  \param max      Maximum allowed value of the parameter in the ins file
 *  \param help     Documentation describing the parameter to the user
 */
void declare_parameter(const char* name, double* param, double min, double max, const char* help = "");

/// Declares a bool parameter
/** \param name     The name of the parameter
 *  \param param    Pointer to variable where the value of the parameter is to be placed
 *  \param help     Documentation describing the parameter to the user
 */
void declare_parameter(const char* name, bool* param, const char* help = "");

#endif // LPJ_GUESS_PARAMETERS_H
