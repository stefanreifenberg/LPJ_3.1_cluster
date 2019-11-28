///////////////////////////////////////////////////////////////////////////////////////
/// \file parameters.cpp
/// \brief Implementation of the parameters module
///
/// New instructions (PLIB keywords) may be added (this would require addition of a
/// declareitem call in function plib_declarations, and possibly some additional code in
/// function plib_callback).
///
/// \author Joe Siltberg
/// $Date: 2015-04-08 13:56:22 +0200 (Mi, 08 Apr 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "parameters.h"
#include "guess.h"
#include "plib.h"
#include <map>

// Definitions of parameters defined globally in parameters.h,
// for documentation, see parameters.h

xtring title;
vegmodetype vegmode;
int npatch;
double patcharea;
bool ifbgestab;
bool ifsme;
bool ifstochestab;
bool ifstochmort;
bool iffire;
bool ifdisturb;
bool ifcalcsla;
bool ifcalccton;
int estinterval;
double distinterval;
bool ifcdebt;

bool ifcentury;
bool ifnlim;
int freenyears;
double nrelocfrac;
double nfix_a;
double nfix_b;

bool ifsmoothgreffmort;
bool ifdroughtlimitedestab;
bool ifrainonwetdaysonly;

bool ifbvoc;

wateruptaketype wateruptake;

bool run_landcover;
bool run[NLANDCOVERTYPES];
bool lcfrac_fixed;
bool all_fracs_const;
bool ifslowharvestpool;
int nyear_spinup;
bool textured_soil;


xtring state_path;
bool restart;
bool save_state;
int state_year;

//SR_20161123
int plant_year;

// f_js_20170126
int ndiamclass;
double diamclasswidth;

///////////////////////////////////////////////////////////////////////////////////////
// Implementation of the Paramlist class

Paramlist param;

void Paramlist::addparam(xtring name,xtring value) {
	Paramtype* p = find(name);
	if (p == 0) {
		p = &createobj();
	}
	p->name=name.lower();
	p->str=value;
}

void Paramlist::addparam(xtring name,double value) {
	Paramtype* p = find(name);
	if (p == 0) {
		p = &createobj();
	}
	p->name=name.lower();
	p->num=value;
}

Paramtype& Paramlist::operator[](xtring name) {
	Paramtype* param = find(name);

	if (param == 0) {
		fail("Paramlist::operator[]: parameter \"%s\" not found",(char*)name);
	}

	return *param;
}

Paramtype* Paramlist::find(xtring name) {
	name = name.lower();
	firstobj();
	while (isobj) {
		Paramtype& p=getobj();
		if (p.name==name) return &p;
		nextobj();
	}
	// nothing found
	return 0;
}


///////////////////////////////////////////////////////////////////////////////////////
// ENUM DECLARATIONS OF INTEGER CONSTANTS FOR PLIB INTERFACE

enum {BLOCK_GLOBAL,BLOCK_PFT,BLOCK_PARAM};
enum {CB_NONE,CB_VEGMODE,CB_CHECKGLOBAL,CB_LIFEFORM,CB_LANDCOVER,CB_PHENOLOGY,CB_LEAFPHYSIOGNOMY,
	CB_PATHWAY,	CB_ROOTDIST,CB_EST,CB_CHECKPFT,CB_STRPARAM,CB_NUMPARAM,CB_WATERUPTAKE,CB_THINNING_REGIME};
// SR_20161205 added CB_THINNING_REGIME above

// File local variables
namespace {

Pft* ppft; // pointer to Pft object currently being assigned to

xtring paramname;
xtring strparam;
double numparam;
bool ifhelp=false;

// 'include' parameter for currently scanned PFT
bool includepft;

// 'include' parameter per PFT
std::map<xtring, bool> includepft_map;

// Whether each PFT has had their parameters checked.
// We only check a PFT:s parameters (in plib_callback) the first time the PFT is
// parsed. If the same PFT occurs again (probably in a different file with a few
// minor modifications) we don't check again (because plib's itemparsed() function
// doesn't remember the old parsed parameters).
std::map<xtring, bool> checked_pft;

}

void initsettings() {

	// Initialises global settings
	// Parameters not initialised here must be set in instruction script

	iffire=true;
	ifcalcsla=true;
	ifdisturb=false;
	ifcalcsla=false;
	ifcalccton=true;
	ifcdebt=false;
	distinterval=1.0e10;
	npatch=1;
	vegmode=COHORT;
	run_landcover = false;

	save_state = false;
	restart = false;
	textured_soil = true;

	//SR_20161123
	plant_year = 0;

  // f_js_20170126
  ndiamclass = 28;
  diamclasswidth = 0.05;
}

void initpft(Pft& pft,xtring& setname) {

	// Initialises a PFT object
	// Parameters not initialised here must be set in instruction script

	pft.name=setname;
	pft.lifeform=NOLIFEFORM;
	pft.phenology=NOPHENOLOGY;

	// Set bioclimatic limits so that PFT can establish and survive under all
	// conditions (may be overridden by settings in instruction script)

	pft.tcmin_surv=-1000.0;
	pft.tcmin_est=-1000.0;
	pft.tcmax_est=1000.0;
	pft.twmin_est=-1000.0;
	pft.gdd5min_est=1000.0;
	pft.twminusc=0.0;

	// Set chilling parameters so that no chilling period required for budburst

	pft.k_chilla=0.0;
	pft.k_chillb=0.0;
	pft.k_chillk=0.0;
}


///////////////////////////////////////////////////////////////////////////////////////
// Data structures for storing parameters declared from other modules
// (see the declare_parameter functions)
//

namespace {

struct xtringParam {

	xtringParam(const char* name, xtring* param, int maxlen, const char* help)
		: name(name),
		  param(param),
		  maxlen(maxlen),
		  help(help) {}

	const char* name;
	xtring* param;
	int maxlen;
	const char* help;
};

std::vector<xtringParam> xtringParams;

struct stringParam {

	stringParam(const char* name, std::string* param, int maxlen, const char* help)
		: name(name),
		  param(param),
		  maxlen(maxlen),
		  help(help) {}

	const char* name;
	std::string* param;
	int maxlen;
	const char* help;
};

std::vector<stringParam> stringParams;

struct intParam {

	intParam(const char* name, int* param, int min, int max, const char* help)
		: name(name),
		  param(param),
		  min(min),
		  max(max),
		  help(help) {}

	const char* name;
	int* param;
	int min;
	int max;
	const char* help;
};

std::vector<intParam> intParams;

struct doubleParam {

	doubleParam(const char* name, double* param, double min, double max, const char* help)
		: name(name),
		  param(param),
		  min(min),
		  max(max),
		  help(help) {}

	const char* name;
	double* param;
	double min;
	double max;
	const char* help;
};

std::vector<doubleParam> doubleParams;

struct boolParam {

	boolParam(const char* name, bool* param, const char* help)
		: name(name),
		  param(param),
		  help(help) {}

	const char* name;
	bool* param;
	const char* help;
};

std::vector<boolParam> boolParams;

} // namespace

///////////////////////////////////////////////////////////////////////////////////////
// The following declare_parameter allow other modules to declare instruction file
// parameters. The information is simply stored and then sent to plib in
// plib_declarations with calls to declareitem().

void declare_parameter(const char* name, xtring* param, int maxlen, const char* help) {
	xtringParams.push_back(xtringParam(name, param, maxlen, help));
}

void declare_parameter(const char* name, std::string* param, int maxlen, const char* help) {
	stringParams.push_back(stringParam(name, param, maxlen, help));
}

void declare_parameter(const char* name, int* param, int min, int max, const char* help) {
	intParams.push_back(intParam(name, param, min, max, help));
}

void declare_parameter(const char* name, double* param, double min, double max, const char* help) {
	doubleParams.push_back(doubleParam(name, param, min, max, help));
}

void declare_parameter(const char* name, bool* param, const char* help) {
	boolParams.push_back(boolParam(name, param, help));
}


///////////////////////////////////////////////////////////////////////////////////////
// The following code uses functionality from the PLIB library to process an
// instruction script (ins) file containing simulation settings and PFT parameters.
// Function read_instruction_file() is called by the framework to initiate parsing of the script.
// Function printhelp() is called if GUESS is run with '-help' instead of an ins file
// name as a command line argument. Functions plib_declarations, plib_callback and
// plib_receivemessage comprise part of the interface to PLIB.

void plib_declarations(int id,xtring setname) {

	switch (id) {

	case BLOCK_GLOBAL:

		declareitem("title",&title,80,CB_NONE,"Title for run");
		declareitem("nyear_spinup",&nyear_spinup,1,10000,1,CB_NONE,"Number of simulation years to spinup for");
		declareitem("vegmode",&strparam,16,CB_VEGMODE,
			"Vegetation mode (\"INDIVIDUAL\", \"COHORT\", \"POPULATION\")");
		declareitem("ifbgestab",&ifbgestab,1,CB_NONE,
			"Whether background establishment enabled (0,1)");
		declareitem("ifsme",&ifsme,1,CB_NONE,
			"Whether spatial mass effect enabled for establishment (0,1)");
		declareitem("ifstochmort",&ifstochmort,1,CB_NONE,
			"Whether mortality stochastic (0,1)");
		declareitem("ifstochestab",&ifstochestab,1,CB_NONE,
			"Whether establishment stochastic (0,1)");
		declareitem("estinterval",&estinterval,1,10,1,CB_NONE,
			"Interval for establishment of new cohorts (years)");
		declareitem("distinterval",&distinterval,1.0,1.0e10,1,CB_NONE,
			"Generic patch-destroying disturbance interval (years)");
		declareitem("iffire",&iffire,1,CB_NONE,
			"Whether fire enabled (0,1)");
		declareitem("ifdisturb",&ifdisturb,1,CB_NONE,
			"Whether generic patch-destroying disturbance enabled (0,1)");
		declareitem("ifcalcsla",&ifcalcsla,1,CB_NONE,
			"Whether SLA calculated from leaf longevity");
		declareitem("ifcalccton",&ifcalccton,1,CB_NONE,
			"Whether leaf C:N min calculated from leaf longevity");
		declareitem("ifcdebt",&ifcdebt,1,CB_NONE,
			"Whether to allow C storage");
		declareitem("npatch",&npatch,1,1000,1,CB_NONE,
			"Number of patches simulated");
		declareitem("patcharea",&patcharea,1.0,1.0e4,1,CB_NONE,
			"Patch area (m2)");
		declareitem("wateruptake", &strparam, 20, CB_WATERUPTAKE,
			"Water uptake mode (\"WCONT\", \"ROOTDIST\", \"SMART\", \"SPECIESSPECIFIC\")");

		declareitem("nrelocfrac",&nrelocfrac,0.0,0.99,1,CB_NONE,
			"Fractional nitrogen relocation from shed leaves & roots");
		declareitem("nfix_a",&nfix_a,0.0,0.4,1,CB_NONE,
			"first term in nitrogen fixation eqn");
		declareitem("nfix_b",&nfix_b,-10.0,10.,1,CB_NONE,
			"second term in nitrogen fixation eqn");

		declareitem("ifcentury",&ifcentury,1,CB_NONE,
			"Whether to use CENTURY SOM dynamics (default standard LPJ)");
		declareitem("ifnlim",&ifnlim,1,CB_NONE,
			"Whether plant growth limited by available nitrogen");
		declareitem("freenyears",&freenyears,0,1000,1,CB_NONE,
			"Number of years to spinup without nitrogen limitation");

		declareitem("ifsmoothgreffmort",&ifsmoothgreffmort,1,CB_NONE,
			"Whether to vary mort_greff smoothly with growth efficiency (0,1)");
		declareitem("ifdroughtlimitedestab",&ifdroughtlimitedestab,1,CB_NONE,
			"Whether establishment drought limited (0,1)");
		declareitem("ifrainonwetdaysonly",&ifrainonwetdaysonly,1,CB_NONE,
			"Whether it rains on wet days only (1), or a little every day (0);");
		//SR_20161123
		declareitem("plant_year", &plant_year, 0, 10000, 1, CB_NONE,
			"year of first plantation");

    // f_js_20170126
    declareitem("ndiamclass", &ndiamclass, 0, 10000, 1, CB_NONE,
			"number of dbh classes for output");
    declareitem("diamclasswidth", &diamclasswidth, 0.01, 1.0, 1, CB_NONE,
			"step size in m for dbh class output");

		// bvoc
		declareitem("ifbvoc",&ifbvoc,1,CB_NONE,
			"Whether or not BVOC calculations are performed (0,1)");
		declareitem("run_landcover",&run_landcover,1,CB_NONE,"Landcover version");
		declareitem("run_urban",&run[URBAN],1,CB_NONE,"Whether urban land is to be simulated");
		declareitem("run_crop",&run[CROPLAND],1,CB_NONE,"Whether crop-land is to be simulated");
		declareitem("run_pasture",&run[PASTURE],1,CB_NONE,"Whether pasture is to be simulated");
		declareitem("run_forest",&run[FOREST],1,CB_NONE,"Whether managed forest is to be simulated");
		declareitem("run_natural",&run[NATURAL],1,CB_NONE,"Whether natural vegetation is to be simulated");
		declareitem("run_peatland",&run[PEATLAND],1,CB_NONE,"Whether peatland is to be simulated");
		declareitem("ifslowharvestpool",&ifslowharvestpool,1,CB_NONE,"If a slow harvested product pool is included in patchpft.");
		declareitem("lcfrac_fixed",&lcfrac_fixed,1,CB_NONE,"Whether static landcover fractions are set in the ins-file (0,1)");

		declareitem("textured_soil",&textured_soil,1,CB_NONE,"Use silt/sand fractions specific to soiltype");

		declareitem("state_path", &state_path, 300, CB_NONE, "State files directory (for restarting from, or saving state files)");
		declareitem("restart", &restart, 1, CB_NONE, "Whether to restart from state files");
		declareitem("save_state", &save_state, 1, CB_NONE, "Whether to save new state files");
		declareitem("state_year", &state_year, 1, 20000, 1, CB_NONE, "Save/restart year. Unspecified means just after spinup");

		declareitem("pft",BLOCK_PFT,CB_NONE,"Header for block defining PFT");
		declareitem("param",BLOCK_PARAM,CB_NONE,"Header for custom parameter block");

		for (size_t i = 0; i < xtringParams.size(); ++i) {
			const xtringParam& p = xtringParams[i];
			declareitem(p.name, p.param, p.maxlen, 0, p.help);
		}

		for (size_t i = 0; i < stringParams.size(); ++i) {
			const stringParam& p = stringParams[i];
			declareitem(p.name, p.param, p.maxlen, 0, p.help);
		}

		for (size_t i = 0; i < intParams.size(); ++i) {
			const intParam& p = intParams[i];
			declareitem(p.name, p.param, p.min, p.max, 1, 0, p.help);
		}

		for (size_t i = 0; i < doubleParams.size(); ++i) {
			const doubleParam& p = doubleParams[i];
			declareitem(p.name, p.param, p.min, p.max, 1, 0, p.help);
		}

		for (size_t i = 0; i < boolParams.size(); ++i) {
			const boolParam& p = boolParams[i];
			declareitem(p.name, p.param, 1, 0, p.help);
		}

		callwhendone(CB_CHECKGLOBAL);


		break;

	case BLOCK_PFT:

		if (!ifhelp) {

			ppft = 0;

			// Was this pft already created?
			for (size_t p = 0; p < pftlist.nobj; ++p) {
				if (pftlist[p].name == setname) {
					ppft = &pftlist[p];
				}
			}

			if (ppft == 0) {
				// Create and initialise a new Pft object and obtain a reference to it

				ppft=&pftlist.createobj();
				initpft(*ppft,setname);
				includepft_map[setname] = true;
			}
		}

		declareitem("include",&includepft,1,CB_NONE,"Include PFT in analysis");
		declareitem("lifeform",&strparam,16,CB_LIFEFORM,
			"Lifeform (\"TREE\" or \"GRASS\")");
		declareitem("landcover",&strparam,16,CB_LANDCOVER,
			"Landcovertype (\"URBAN\", \"CROP\", \"PASTURE\", \"FOREST\", \"NATURAL\" or \"PEATLAND\")");
		declareitem("phenology",&strparam,16,CB_PHENOLOGY,
			"Phenology (\"EVERGREEN\", \"SUMMERGREEN\", \"RAINGREEN\" or \"ANY\")");
		declareitem("leafphysiognomy",&strparam,16,CB_LEAFPHYSIOGNOMY,
			"Leaf physiognomy (\"NEEDLELEAF\" or \"BROADLEAF\")");
		declareitem("phengdd5ramp",&ppft->phengdd5ramp,0.0,1000.0,1,CB_NONE,
			"GDD on 5 deg C base to attain full leaf cover");
		declareitem("wscal_min",&ppft->wscal_min,0.0,1.0,1,CB_NONE,
			"Water stress threshold for leaf abscission (raingreen PFTs)");
		declareitem("pathway",&strparam,16,CB_PATHWAY,
			"Biochemical pathway (\"C3\" or \"C4\")");
		declareitem("pstemp_min",&ppft->pstemp_min,-50.0,50.0,1,CB_NONE,
			"Approximate low temp limit for photosynthesis (deg C)");
		declareitem("pstemp_low",&ppft->pstemp_low,-50.0,50.0,1,CB_NONE,
			"Approx lower range of temp optimum for photosynthesis (deg C)");
		declareitem("pstemp_high",&ppft->pstemp_high,0.0,60.0,1,CB_NONE,
			"Approx higher range of temp optimum for photosynthesis (deg C)");
		declareitem("pstemp_max",&ppft->pstemp_max,0.0,60.0,1,CB_NONE,
			"Maximum temperature limit for photosynthesis (deg C)");
		declareitem("lambda_max",&ppft->lambda_max,0.1,0.99,1,CB_NONE,
			"Non-water-stressed ratio of intercellular to ambient CO2 pp");
		declareitem("rootdist",ppft->rootdist,0.0,1.0,NSOILLAYER,CB_ROOTDIST,
			"Fraction of roots in each soil layer (first value=upper layer)");
		declareitem("gmin",&ppft->gmin,0.0,1.0,1,CB_NONE,
			"Canopy conductance not assoc with photosynthesis (mm/s)");
		declareitem("emax",&ppft->emax,0.0,50.0,1,CB_NONE,
			"Maximum evapotranspiration rate (mm/day)");
		// guess2008 - increased the upper limit to possible respcoeff values (was 1.2)
		declareitem("respcoeff",&ppft->respcoeff,0.0,3,1,CB_NONE,
			"Respiration coefficient (0-1)");

		declareitem("cton_root",&ppft->cton_root,1.0,1.0e4,1,CB_NONE,
			"Reference Fine root C:N mass ratio");
		declareitem("cton_sap",&ppft->cton_sap,1.0,1.0e4,1,CB_NONE,
			"Reference Sapwood C:N mass ratio");
		declareitem("nuptoroot",&ppft->nuptoroot,0.0,1.0,1,CB_NONE,
			"Maximum nitrogen uptake per fine root");
		declareitem("km_volume",&ppft->km_volume,0.0,10.0,1,CB_NONE,
			"Michaelis-Menten kinetic parameters for nitrogen uptake");
		declareitem("fnstorage",&ppft->fnstorage,0.0,10.0,1,CB_NONE,
			"fraction of sapwood (root for herbaceous pfts) that can be used as a nitrogen storage scalar");

		declareitem("reprfrac",&ppft->reprfrac,0.0,1.0,1,CB_NONE,
			"Fraction of NPP allocated to reproduction");
		declareitem("turnover_leaf",&ppft->turnover_leaf,0.0,1.0,1,CB_NONE,
			"Leaf turnover (fraction/year)");
		declareitem("turnover_root",&ppft->turnover_root,0.0,1.0,1,CB_NONE,
			"Fine root turnover (fraction/year)");
		declareitem("turnover_sap",&ppft->turnover_sap,0.0,1.0,1,CB_NONE,
			"Sapwood turnover (fraction/year)");
		declareitem("wooddens",&ppft->wooddens,10.0,1000.0,1,CB_NONE,
			"Sapwood and heartwood density (kgC/m3)");
		declareitem("crownarea_max",&ppft->crownarea_max,1.0,1000.0,1,CB_NONE,
			"Maximum tree crown area (m2)");
		declareitem("k_allom1",&ppft->k_allom1,10.0,1000.0,1,CB_NONE,
			"Constant in allometry equations");
		// guess2008 - changed lower limit for k_allom2 to 1 from 10. This is needed
		// for the shrub allometries.
		declareitem("k_allom2",&ppft->k_allom2,1.0,1.0e4,1,CB_NONE,
			"Constant in allometry equations");
		declareitem("k_allom3",&ppft->k_allom3,0.1,1.0,1,CB_NONE,
			"Constant in allometry equations");
		declareitem("k_rp",&ppft->k_rp,1.0,2.0,1,CB_NONE,
			"Constant in allometry equations");
		declareitem("k_latosa",&ppft->k_latosa,100.0,1.0e5,1,CB_NONE,
			"Tree leaf to sapwood xs area ratio");
		declareitem("sla",&ppft->sla,1.0,1000.0,1,CB_NONE,
			"Specific leaf area (m2/kgC)");
		declareitem("cton_leaf_min",&ppft->cton_leaf_min,1.0,1.0e4,1,CB_NONE,
			"Minimum leaf C:N mass ratio");
		declareitem("ltor_max",&ppft->ltor_max,0.1,10.0,1,CB_NONE,
			"Non-water-stressed leaf:fine root mass ratio");
		declareitem("litterme",&ppft->litterme,0.0,1.0,1,CB_NONE,
			"Litter moisture flammability threshold (fraction of AWC)");
		declareitem("fireresist",&ppft->fireresist,0.0,1.0,1,CB_NONE,
			"Fire resistance (0-1)");
		declareitem("tcmin_surv",&ppft->tcmin_surv,-1000.0,50.0,1,CB_NONE,
			"Min 20-year coldest month mean temp for survival (deg C)");
		declareitem("tcmin_est",&ppft->tcmin_est,-1000.0,50.0,1,CB_NONE,
			"Min 20-year coldest month mean temp for establishment (deg C)");
		declareitem("tcmax_est",&ppft->tcmax_est,-50.0,1000.0,1,CB_NONE,
			"Max 20-year coldest month mean temp for establishment (deg C)");
		declareitem("twmin_est",&ppft->twmin_est,-1000.0,50.0,1,CB_NONE,
			"Min warmest month mean temp for establishment (deg C)");
		declareitem("twminusc",&ppft->twminusc,0,100,1,CB_NONE,
			"Stupid larch parameter");
		declareitem("gdd5min_est",&ppft->gdd5min_est,0.0,5000.0,1,CB_NONE,
			"Min GDD on 5 deg C base for establishment");
		declareitem("k_chilla",&ppft->k_chilla,0.0,5000.0,1,CB_NONE,
			"Constant in equation for budburst chilling time requirement");
		declareitem("k_chillb",&ppft->k_chillb,0.0,5000.0,1,CB_NONE,
			"Coefficient in equation for budburst chilling time requirement");
		declareitem("k_chillk",&ppft->k_chillk,0.0,1.0,1,CB_NONE,
			"Exponent in equation for budburst chilling time requirement");
		declareitem("parff_min",&ppft->parff_min,0.0,1.0e7,1,CB_NONE,
			"Min forest floor PAR for grass growth/tree estab (J/m2/day)");
		declareitem("alphar",&ppft->alphar,0.01,100.0,1,CB_NONE,
			"Shape parameter for recruitment-juv growth rate relationship");
		declareitem("est_max",&ppft->est_max,1.0e-4,1.0,1,CB_NONE,
			"Max sapling establishment rate (indiv/m2/year)");
		declareitem("kest_repr",&ppft->kest_repr,1.0,1000.0,1,CB_NONE,
			"Constant in equation for tree estab rate");
		declareitem("kest_bg",&ppft->kest_bg,0.0,1.0,1,CB_NONE,
			"Constant in equation for tree estab rate");
		declareitem("kest_pres",&ppft->kest_pres,0.0,1.0,1,CB_NONE,
			"Constant in equation for tree estab rate");
		declareitem("longevity",&ppft->longevity,0.0,3000.0,1,CB_NONE,
			"Expected longevity under lifetime non-stressed conditions (yr)");
		declareitem("greff_min",&ppft->greff_min,0.0,1.0,1,CB_NONE,
			"Threshold for growth suppression mortality (kgC/m2 leaf/yr)");
		declareitem("leaflong",&ppft->leaflong,0.1,100.0,1,CB_NONE,
			"Leaf longevity (years)");
		declareitem("intc",&ppft->intc,0.0,1.0,1,CB_NONE,"Interception coefficient");

		
		declareitem("rotation_interv", &ppft->rotation_interv, 0, 10000, 1, CB_NONE,
			"final harvesting interval");
		declareitem("thinning_interv", &ppft->thinning_interv, 0, 10000, 1, CB_NONE,
			"thinning interval in years"); // SR 20161124
		declareitem("thinning_frac", &ppft->thinning_frac, 0, 1.0, 1, CB_NONE,
			"fraction of BA removed");
		declareitem("thinning_regime", &strparam, 16, CB_THINNING_REGIME,
			"Thinning regime (\"NONE\",\"COMBINED\", \"ABOVE\" or \"BELOW\")");
		declareitem("min_diam", &ppft->min_diam, 0, 1.0, 1, CB_NONE,
			"Minimum Diameter");
		declareitem("max_diam", &ppft->max_diam, 0, 1.0, 1, CB_NONE,
			"Maximum Diameter");
		declareitem("harvest_frac", &ppft->harvest_frac, 0, 1.0, 1, CB_NONE,
			"harvested fraction"); //SR 20170523
		declareitem("est_indiv", &ppft->est_indiv, 0, 1.0, 1, CB_NONE,
			"number of planted individuals");

		
		// guess2008 - DLE
		declareitem("drought_tolerance",&ppft->drought_tolerance,0.0,1.0,1,CB_NONE,
			"Drought tolerance level (0 = very -> 1 = not at all) (unitless)");

		// bvoc
		declareitem("ga",&ppft->ga,0.0,1.0,1,CB_NONE,
			"aerodynamic conductance (m/s)");
		declareitem("eps_iso",&ppft->eps_iso,0.,100.,1,CB_NONE,
			"isoprene emission capacity (ug C g-1 h-1)");
		declareitem("seas_iso",&ppft->seas_iso,1,CB_NONE,
			"whether (1) or not (0) isoprene emissions show seasonality");
		declareitem("eps_mon",&ppft->eps_mon,0.,100.,1,CB_NONE,
			"monoterpene emission capacity (ug C g-1 h-1)");
		declareitem("storfrac_mon",&ppft->storfrac_mon,0.,1.,1,CB_NONE,
			"fraction of monoterpene production that goes into storage pool (-)");

		declareitem("harv_eff",&ppft->harv_eff,0.0,1.0,1,CB_NONE,"Harvest efficiency");
		declareitem("harvest_slow_frac",&ppft->harvest_slow_frac,0.0,1.0,1,CB_NONE,
			"Fraction of harvested products that goes into carbon depository for long-lived products like wood");
		declareitem("turnover_harv_prod",&ppft->turnover_harv_prod,0.0,1.0,1,CB_NONE,"Harvested products turnover (fraction/year)");
		declareitem("res_outtake",&ppft->res_outtake,0.0,1.0,1,CB_NONE,"Fraction of residue outtake at harvest");

		callwhendone(CB_CHECKPFT);

		break;

	case BLOCK_PARAM:

		paramname=setname;
		declareitem("str",&strparam,300,CB_STRPARAM,
			"String value for custom parameter");
		declareitem("num",&numparam,-1.0e38,1.0e38,1,CB_NUMPARAM,
			"Numerical value for custom parameter");

		break;
	}
}

void badins(xtring missing) {

	xtring message=(xtring)"Missing mandatory setting: "+missing;
	sendmessage("Error",message);
	plibabort();
}

void plib_callback(int callback) {

	xtring message;
	int i;
	double numval;

	switch (callback) {

	case CB_VEGMODE:
		if (strparam.upper()=="INDIVIDUAL") vegmode=INDIVIDUAL;
		else if (strparam.upper()=="COHORT") vegmode=COHORT;
		else if (strparam.upper()=="POPULATION") vegmode=POPULATION;
		else {
			sendmessage("Error",
				"Unknown vegetation mode (valid types: \"INDIVIDUAL\",\"COHORT\", \"POPULATION\")");
			plibabort();
		}
		break;
	case CB_WATERUPTAKE:
		if (strparam.upper() == "WCONT") wateruptake = WR_WCONT;
		else if (strparam.upper() == "ROOTDIST") wateruptake = WR_ROOTDIST;
		else if (strparam.upper() == "SMART") wateruptake = WR_SMART;
		else if (strparam.upper() == "SPECIESSPECIFIC") wateruptake = WR_SPECIESSPECIFIC;
		else {
			sendmessage("Error",
				"Unknown water uptake mode (valid types: \"WCONT\", \"ROOTDIST\", \"SMART\", \"SPECIESSPECIFIC\")");
		}
		break;
	case CB_LIFEFORM:
		if (strparam.upper()=="TREE") ppft->lifeform=TREE;
		else if (strparam.upper()=="GRASS") ppft->lifeform=GRASS;
		else {
			sendmessage("Error",
				"Unknown lifeform type (valid types: \"TREE\", \"GRASS\")");
			plibabort();
		}
		break;
	case CB_LANDCOVER:
		if (strparam.upper()=="NATURAL") ppft->landcover=NATURAL;
		else if (strparam.upper()=="URBAN") ppft->landcover=URBAN;
		else if (strparam.upper()=="CROPLAND") ppft->landcover=CROPLAND;
		else if (strparam.upper()=="PASTURE") ppft->landcover=PASTURE;
		else if (strparam.upper()=="FOREST") ppft->landcover=FOREST;
		else if (strparam.upper()=="PEATLAND") ppft->landcover=PEATLAND;
		else {
			sendmessage("Error",
				"Unknown landcover type (valid types: \"URBAN\", \"CROPLAND\", \"PASTURE\", \"FOREST\", \"NATURAL\" or \"PEATLAND\")");
			plibabort();
		}
		break;
	case CB_PHENOLOGY:
		if (strparam.upper()=="SUMMERGREEN") ppft->phenology=SUMMERGREEN;
		else if (strparam.upper()=="RAINGREEN") ppft->phenology=RAINGREEN;
		else if (strparam.upper()=="EVERGREEN") ppft->phenology=EVERGREEN;
		else if (strparam.upper()=="ANY") ppft->phenology=ANY;
		else {
			sendmessage("Error",
				"Unknown phenology type\n  (valid types: \"EVERGREEN\", \"SUMMERGREEN\", \"RAINGREEN\" or \"ANY\")");
			plibabort();
		}
		break;
	case CB_LEAFPHYSIOGNOMY:
		if (strparam.upper()=="NEEDLELEAF") ppft->leafphysiognomy=NEEDLELEAF;
		else if (strparam.upper()=="BROADLEAF") ppft->leafphysiognomy=BROADLEAF;
		else {
			sendmessage("Error",
				"Unknown leaf physiognomy (valid types: \"NEEDLELEAF\", \"BROADLEAF\")");
			plibabort();
		}
		break;
	case CB_PATHWAY:
		if (strparam.upper()=="C3") ppft->pathway=C3;
		else if (strparam.upper()=="C4") ppft->pathway=C4;
		else {
			sendmessage("Error",
				"Unknown pathway type\n  (valid types: \"C3\" or \"C4\")");
			plibabort();
		}
		break;
	case CB_ROOTDIST:
		numval=0.0;
		for (i=0;i<NSOILLAYER;i++) numval+=ppft->rootdist[i];
		if (numval<0.99 || numval>1.01) {
			sendmessage("Error","Specified root fractions do not sum to 1.0");
			plibabort();
		}
		ppft->rootdist[NSOILLAYER-1]+=1.0-numval;
		break;
	case CB_STRPARAM:
		param.addparam(paramname,strparam);
		break;
	case CB_NUMPARAM:
		param.addparam(paramname,numparam);
		break;
	case CB_CHECKGLOBAL:
		if (!itemparsed("title")) badins("title");
		if (!itemparsed("nyear_spinup")) badins("nyear_spinup");
		if (!itemparsed("vegmode")) badins("vegmode");
		if (!itemparsed("iffire")) badins("iffire");
		if (!itemparsed("ifcalcsla")) badins("ifcalcsla");
		if (!itemparsed("ifcalccton")) badins("ifcalccton");
		if (!itemparsed("ifcdebt")) badins("ifcdebt");
		if (!itemparsed("wateruptake")) badins("wateruptake");

		if (!itemparsed("nrelocfrac")) badins("nrelocfrac");
		if (!itemparsed("nfix_a")) badins("nfix_a");
		if (!itemparsed("nfix_b")) badins("nfix_b");

		if (!itemparsed("ifcentury")) badins("ifcentury");
		if (!itemparsed("ifnlim")) badins("ifnlim");
		if (!itemparsed("freenyears")) badins("freenyears");

		if (nyear_spinup <= freenyears) {
			sendmessage("Error", "freenyears must be smaller than nyear_spinup");
			plibabort();
		}
		//SR_20161123
		if (!itemparsed("plant_year")) badins("plant_year");

    // f_js_20170126
    if (!itemparsed("ndiamclass")) badins("ndiamclass");
    if (!itemparsed("diamclasswidth")) badins("diamclasswidth");

		if (!itemparsed("outputdirectory")) badins("outputdirectory");
		if (!itemparsed("ifsmoothgreffmort")) badins("ifsmoothgreffmort");
		if (!itemparsed("ifdroughtlimitedestab")) badins("ifdroughtlimitedestab");
		if (!itemparsed("ifrainonwetdaysonly")) badins("ifrainonwetdaysonly");
		// bvoc
		if (!itemparsed("ifbvoc")) badins("ifbvoc");

		if (!itemparsed("run_landcover")) badins("run_landcover");
		if (run_landcover) {
			if (!itemparsed("lcfrac_fixed")) badins("lcfrac_fixed");
			if (!itemparsed("equal_landcover_area")) badins("equal_landcover_area");
			if (!itemparsed("lc_fixed_urban")) badins("lc_fixed_urban");
			if (!itemparsed("lc_fixed_cropland")) badins("lc_fixed_cropland");
			if (!itemparsed("lc_fixed_pasture")) badins("lc_fixed_pasture");
			if (!itemparsed("lc_fixed_forest")) badins("lc_fixed_forest");
			if (!itemparsed("lc_fixed_natural")) badins("lc_fixed_natural");
			if (!itemparsed("lc_fixed_peatland")) badins("lc_fixed_peatland");
			if (!itemparsed("run_natural")) badins("run_natural");
			if (!itemparsed("run_crop")) badins("run_crop");
			if (!itemparsed("run_forest")) badins("run_forest");
			if (!itemparsed("run_urban")) badins("run_urban");
			if (!itemparsed("run_pasture")) badins("run_pasture");
			if (!itemparsed("ifslowharvestpool")) badins("ifslowharvestpool");
		}

		if (!itemparsed("pft")) badins("pft");
		if (vegmode==COHORT || vegmode==INDIVIDUAL) {
			if (!itemparsed("ifbgestab")) badins("ifbgestab");
			if (!itemparsed("ifsme")) badins("ifsme");
			if (!itemparsed("ifstochmort")) badins("ifstochmort");
			if (!itemparsed("ifstochestab")) badins("ifstochestab");
			if (itemparsed("ifdisturb") && !itemparsed("distinterval"))
				badins("distinterval");
			if (!itemparsed("npatch")) badins("npatch");
			if (!itemparsed("patcharea")) badins("patcharea");
			if (!itemparsed("estinterval")) badins("estinterval");
		}
		else if (vegmode==POPULATION && npatch!=1) {
			sendmessage("Information",
				"Value specified for npatch ignored in population mode");
			npatch=1;
		}

		if (save_state && restart) {
			sendmessage("Error",
			            "Can't save state and restart at the same time");
			plibabort();
		}

		if (!itemparsed("state_year")) {
			state_year = nyear_spinup;
		}

		if (state_path == "" && (save_state || restart)) {
			badins("state_path");
		}

		//	delete unused pft:s from pftlist

		pftlist.firstobj();
		while (pftlist.isobj) {
			Pft& pft = pftlist.getobj();
			bool include = includepft_map[pft.name];

			if (pft.landcover!=NATURAL) {
				if (!run_landcover || !run[pft.landcover])
					include = false;
			}
			else if (run_landcover && !run[NATURAL]) {
				if (pft.landcover==NATURAL)
					include = false;
			}

			if (!include) {
				// Remove this PFT from list
				pftlist.killobj();
			}
			else {
				pftlist.nextobj();
			}
		}

		// Set ids and npft variable after removing unused pfts
		npft = 0;
		pftlist.firstobj();
		while (pftlist.isobj) {
			pftlist.getobj().id = npft++;
			pftlist.nextobj();
		}

		// Call various init functions on each PFT now that all settings have been read
		pftlist.firstobj();
		while (pftlist.isobj) {
			ppft = &pftlist.getobj();

			if (ifcalcsla) {
				// Calculate SLA
				ppft->initsla();
			}

			if (ifcalccton) {
				// Calculate leaf C:N ratio minimum
				ppft->init_cton_min();
			}

			// Calculate C:N ratio limits
			ppft->init_cton_limits();

			// Calculate nitrogen uptake strength dependency on root distribution
			ppft->init_nupscoeff();

			// Calculate regeneration characteristics for population mode
			ppft->initregen();

			pftlist.nextobj();
		}


		break;
	case CB_CHECKPFT:
		if (!checked_pft[ppft->name]) {
			checked_pft[ppft->name] = true;

			if (!itemparsed("lifeform")) badins("lifeform");
			if (!itemparsed("phenology")) badins("phenology");
			if (ppft->phenology==SUMMERGREEN || ppft->phenology==ANY)
				if (!itemparsed("phengdd5ramp")) badins("phengdd5ramp");
			if (ppft->phenology==RAINGREEN || ppft->phenology==ANY)
				if (!itemparsed("wscal_min")) badins("wscal_min");
			if (!itemparsed("leafphysiognomy")) badins("leafphysiognomy");
			if (!itemparsed("pathway")) badins("pathway");
			if (!itemparsed("pstemp_min")) badins("pstemp_min");
			if (!itemparsed("pstemp_low")) badins("pstemp_low");
			if (!itemparsed("pstemp_high")) badins("pstemp_high");
			if (!itemparsed("pstemp_max")) badins("pstemp_max");
			if (!itemparsed("lambda_max")) badins("lambda_max");
			if (!itemparsed("rootdist")) badins("rootdist");
			if (!itemparsed("gmin")) badins("gmin");
			if (!itemparsed("emax")) badins("emax");
			if (!itemparsed("respcoeff")) badins("respcoeff");
			if (!itemparsed("sla") && !ifcalcsla) badins("sla");
			if (!itemparsed("cton_leaf_min") && !ifcalccton) badins("cton_leaf_min");

			if (!itemparsed("cton_root")) badins("cton_root");
			if (!itemparsed("nuptoroot")) badins("nuptoroot");
			if (!itemparsed("km_volume")) badins("km_volume");
			if (!itemparsed("fnstorage")) badins("fnstorage");

			if (!itemparsed("reprfrac")) badins("reprfrac");
			if (!itemparsed("turnover_leaf")) badins("turnover_leaf");
			if (!itemparsed("turnover_root")) badins("turnover_root");
			if (!itemparsed("ltor_max")) badins("ltor_max");
			if (!itemparsed("intc")) badins("intc");

			if (run_landcover)
				{
					if (!itemparsed("landcover")) badins("landcover");
					if (!itemparsed("turnover_harv_prod")) badins("turnover_harv_prod");
					if (!itemparsed("harvest_slow_frac")) badins("harvest_slow_frac");
					if (!itemparsed("harv_eff")) badins("harv_eff");
					if (!itemparsed("res_outtake")) badins("res_outtake");
				}

			// guess2008 - DLE
			if (!itemparsed("drought_tolerance")) badins("drought_tolerance");

			// SR_20161125
			if (!itemparsed("rotation_interv")) badins("rotation_interv");
			if (!itemparsed("thinning_interv")) badins("thinning_interv");
			if (!itemparsed("thinning_frac")) badins("thinning_frac");
			if (!itemparsed("thinning_regime")) badins("thinning_regime");
			if (!itemparsed("min_diam")) badins("min_diam");
			if (!itemparsed("max_diam")) badins("max_diam");
			if (!itemparsed("harvest_frac")) badins("harvest_frac");
			//if (!itemparsed("est_indiv")) badins("est_indiv");


			// bvoc
			if(ifbvoc){
				if (!itemparsed("ga")) badins("ga");
				if (!itemparsed("eps_iso")) badins("eps_iso");
				if (!itemparsed("seas_iso")) badins("seas_iso");
				if (!itemparsed("eps_mon")) badins("eps_mon");
				if (!itemparsed("storfrac_mon")) badins("storfrac_mon");
			}

			if (ppft->lifeform==TREE) {
				if (!itemparsed("cton_sap")) badins("cton_sap");
				if (!itemparsed("turnover_sap")) badins("turnover_sap");
				if (!itemparsed("wooddens")) badins("wooddens");
				if (!itemparsed("crownarea_max")) badins("crownarea_max");
				if (!itemparsed("k_allom1")) badins("k_allom1");
				if (!itemparsed("k_allom2")) badins("k_allom2");
				if (!itemparsed("k_allom3")) badins("k_allom3");
				if (!itemparsed("k_rp")) badins("k_rp");
				if (!itemparsed("k_latosa")) badins("k_latosa");
				if (vegmode==COHORT || vegmode==INDIVIDUAL) {
					if (!itemparsed("kest_repr")) badins("kest_repr");
					if (!itemparsed("kest_bg")) badins("kest_bg");
					if (!itemparsed("kest_pres")) badins("kest_pres");
					if (!itemparsed("longevity")) badins("longevity");
					if (!itemparsed("greff_min")) badins("greff_min");
					if (!itemparsed("alphar")) badins("alphar");
					if (!itemparsed("est_max")) badins("est_max");
				}
			}
			if (iffire) {
				if (!itemparsed("litterme")) badins("litterme");
				if (!itemparsed("fireresist")) badins("fireresist");
			}
			if (ifcalcsla) {
				if (!itemparsed("leaflong")) {
					sendmessage("Error",
					            "Value required for leaflong when ifcalcsla enabled");
					plibabort();
				}
				if (itemparsed("sla"))
					sendmessage("Warning",
					            "Specified sla value not used when ifcalcsla enabled");
			}
			if (vegmode==COHORT || vegmode==INDIVIDUAL) {
				if (!itemparsed("parff_min")) badins("parff_min");
			}

			if (ifcalccton) {
				if (!itemparsed("leaflong")) {
					sendmessage("Error",
					            "Value required for leaflong when ifcalccton enabled");
					plibabort();
				}
				if (itemparsed("cton_leaf_min"))
					sendmessage("Warning",
					            "Specified cton_leaf_min value not used when ifcalccton enabled");
			}
		}
		else {
			// This PFT has already been parsed once, don't allow changing parameters
			// which would have incurred different checks above.
			if (itemparsed("lifeform") ||
			    itemparsed("phenology")) {
				sendmessage("Error",
				            "Not allowed to redefine lifeform or phenology in second PFT definition");
				plibabort();
			}
		}

		if (itemparsed("include")) {
			includepft_map[ppft->name] = includepft;
		}

		break;

		// f_js_20160920
	case CB_THINNING_REGIME:
		if (strparam.upper() == "ABOVE") ppft->thinning_regime = ABOVE;
		else if (strparam.upper() == "BELOW") ppft->thinning_regime = BELOW;
		else if (strparam.upper() == "COMBINED") ppft->thinning_regime = COMBINED;
		else {
			ppft->thinning_regime = NONE;
		}
		break;
	}
}

void plib_receivemessage(xtring text) {

	// Output of messages to user sent by PLIB

	dprintf((char*)text);
}

void read_instruction_file(const char* insfilename) {
	if (!fileexists(insfilename)) {
		fail("Error: could not open %s for input",(const char*)insfilename);
	}

	// Initialise PFT count
	npft=0;

	initsettings();
	param.killall();

	// Initialise simulation settings and PFT parameters from instruction script
	if (!plib(insfilename)) {
		fail("Bad instruction file!");
	}
}

void printhelp() {

	// Calls PLIB to output help text

	ifhelp=true;
	plibhelp();
	ifhelp=false;
}
