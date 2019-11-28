///////////////////////////////////////////////////////////////////////////////////////
/// \file demoinput.cpp
/// \brief LPJ-GUESS input module for a toy data set (for demonstration purposes)
///
///
/// \author Ben Smith
/// $Date: 2014-03-13 19:00:46 +0100 (Thu, 13 Mar 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "asciiinput2.h"
#include "driver.h"
#include "config.h"

REGISTER_INPUT_MODULE("ascii2", AsciiInput2)

AsciiInput2::AsciiInput2()
	: ndep_timeseries("historic"),
	lc_fixed_frac(NLANDCOVERTYPES, 0),
	equal_landcover_area(false) {

	// Declare instruction file parameters

	declare_parameter("ndep_timeseries", &ndep_timeseries, 10, "Nitrogen deposition time series to use (historic, rcp26, rcp45, rcp60 or rcp85");

	// Not used by this input module currently, but included as parameters so
	// common ins files can be used.
	declare_parameter("equal_landcover_area", &equal_landcover_area, "Whether enforced static landcover fractions are equal-sized stands of all included landcovers (0,1)");
	declare_parameter("lc_fixed_urban", &lc_fixed_frac[URBAN], 0, 100, "% lc_fixed_urban");
	declare_parameter("lc_fixed_cropland", &lc_fixed_frac[CROPLAND], 0, 100, "% lc_fixed_cropland");
	declare_parameter("lc_fixed_pasture", &lc_fixed_frac[PASTURE], 0, 100, "% lc_fixed_pasture");
	declare_parameter("lc_fixed_forest", &lc_fixed_frac[FOREST], 0, 100, "% lc_fixed_forest");
	declare_parameter("lc_fixed_natural", &lc_fixed_frac[NATURAL], 0, 100, "% lc_fixed_natural");
	declare_parameter("lc_fixed_peatland", &lc_fixed_frac[PEATLAND], 0, 100, "% lc_fixed_peatland");
}



// First and last historical year obtained from the input meteo file
int first_hist_year, last_hist_year;

void AsciiInput2::init() {

		
	double dlon,dlat;
	bool eof=false;
	xtring descrip;

	// Read list of grid coordinates and store in global Coord object 'gridlist'

	// Retrieve name of grid list file as read from ins file
	xtring file_gridlist=param["file_gridlist"].str;

	FILE* in_grid=fopen(file_gridlist,"r");
	if (!in_grid) fail("initio: could not open %s for input",(char*)file_gridlist);

	ngridcell=0;
	while (!eof) {
		
		// Read next record in file
		eof=!readfor(in_grid,"f,f,a#",&dlon,&dlat,&descrip);

		if (!eof && !(dlon==0.0 && dlat==0.0)) { // ignore blank lines at end (if any)
			Coord& c=gridlist.createobj(); // add new coordinate to grid list

			c.lon=dlon;
			c.lat=dlat;
			c.descrip=descrip;
			ngridcell++;
		}
	}

	fclose(in_grid);	
	
	

	// Read CO2 data from file
	co2.load_file(param["file_co2"].str);

	file_cru = param["file_cru"].str;

	// Set timers
	tprogress.init();
	tmute.init();

	tprogress.settimer();
	tmute.settimer(MUTESEC);
}

void AsciiInput2::site_specific_met(double lon, double lat) {

	xtring latstring, lonstring;
	int lines=0;
	sprintf(lonstring,"%d",(int)(lon*10));
	sprintf(latstring,"%d",(int)(lat*10));
	xtring fname = param["file_met"].str + lonstring + latstring + ".txt";
	std::ifstream ifs(fname, std::ifstream::in);

	if (!ifs.good()) {
		fail("Could not open %s for input", (char*)fname);
	}

	int i = 0;
	std::string line;
	while (getline(ifs, line)) {

		std::istringstream iss(line);

		int year, doy;
		double radiation, prec, temp;
		if (iss >> year >> doy >> radiation >> prec >> temp) {
			if (!i) {
				date.set_first_calendar_year(year - nyear_spinup);
			}

			vrad.push_back(radiation);
			vprec.push_back(prec);
			vtemp.push_back(temp);
			i++;
		}
	}
	ifs.close();
	std::vector<double>::size_type days = vrad.size();
	if (days % 365) {
		fail("Given timeseries doesn't extend for a full number of years (length: %d)\n", days);
	}
	nyear_hist = days / 365;

}


/// Called by the framework at the start of the simulation for a particular grid cell
bool AsciiInput2::getgridcell(Gridcell& gridcell) {

	static bool first_call = true;

	if (first_call) {
		gridlist.firstobj();
		first_call = false;
	}
	else gridlist.nextobj();

	if (gridlist.isobj) {

		Coord& c = gridlist.getobj();

		// Find nearest CRU grid cell in order to get the soilcode
		double cru_lon = c.lon;
		double cru_lat = c.lat;
		double dummy[CRU_TS30::NYEAR_HIST][12];

		const double searchradius = 1;

		int soilcode;
		if (!CRU_TS30::findnearestCRUdata(searchradius, file_cru, cru_lon, cru_lat,
					soilcode, dummy, dummy, dummy)) {
			fail("Failed to find soil code from CRU archive, close to coordinates (%g,%g).\n",
					cru_lon, cru_lat);
		}

		// Get nitrogen deposition, using the found CRU coordinates
		ndep.getndep(param["file_ndep"].str, cru_lon, cru_lat,
						Lamarque::parse_timeseries(ndep_timeseries));

		site_specific_met(c.lon, c.lat);


		dprintf("\nCommencing simulation for stand at (%g, %g)", c.lon, c.lat);
		if (c.descrip!="") {
			dprintf(" (%s)", c.descrip.c_str());
		}
		dprintf("\n\nUsing soil code and Nitrogen deposition for (%3.1f, %3.1f)\n", c.lon, c.lat);

		// Tell framework the coordinates of this grid cell
		gridcell.set_coordinates(gridlist.getobj().lon, gridlist.getobj().lat);

		// The insolation data will be sent (in function getclimate, below)
		// as percentage sunshine

		gridcell.climate.instype=SWRAD;

		// Tell framework the soil type of this grid cell
		soilparameters(gridcell.soiltype,soilcode);

		// For Windows shell - clear graphical output
		// (ignored on other platforms)

		clear_all_graphs();

		return true; // simulate this stand
	}
	return false; // no more stands
}

///	Gets gridcell.landcoverfrac from landcover input file(s) for one year or from ins-file .
void AsciiInput2::getlandcover(Gridcell& gridcell) {

	// Only 100% natural land cover is supported by this input module for now

	for (int i = 0; i < NLANDCOVERTYPES; ++i) {
		gridcell.landcoverfrac[i] = 0;
	}
	gridcell.landcoverfrac[NATURAL] = 1;
}


/// Called by the framework each simulation day before any process modelling is performed for this day
/** Obtains climate data (including atmospheric CO2 and insolation) for this day. */
bool AsciiInput2::getclimate(Gridcell& gridcell) {

	// DESCRIPTION
	// The function should returns false if the simulation is complete for this grid cell,
	// otherwise true. This will normally require querying the year and day member
	// variables of the global class object date:
	//
	// if (date.day==0 && date.year==nyear) return false;
	// // else
	// return true;
	//
	// Currently the following member variables of the climate member of gridcell must be
	// initialised: co2, temp, prec, insol.

	Climate& climate = gridcell.climate;

	const int NYEAR_SPINUP = 30; //Maximum number of years in climate data to repeat during spin-up

	// First day of year
	if (date.day == 0) {

		if (date.year == nyear_spinup + nyear_hist) {
			vrad.clear();
			vprec.clear();
			vtemp.clear();
			return false;
		}

		// Get nitrogen data and distribute it to daily array
		double mndrydep[12];
		double mnwetdep[12];
		year_in_file = date.year < nyear_spinup ? date.year % min(NYEAR_SPINUP,nyear_hist) : date.year - nyear_spinup;
		dprec = &vprec[year_in_file * 365];
		ndep.get_one_calendar_year(date.get_calendar_year(), mndrydep, mnwetdep);
		distribute_ndep(mndrydep, mnwetdep, dprec, dndep);
	}


	// Send climate data
	size_t daily_index = year_in_file * 365 + date.day;
	climate.prec = vprec[daily_index]; // Normal precipitation
	climate.insol = vrad[daily_index];
	climate.temp = vtemp[daily_index];


	// Send nitrogen and CO2 data
	climate.dndep = dndep[date.day];
	climate.dnfert = 0.0;
	climate.co2 = co2[date.get_calendar_year()];

	
	// bvoc
	//climate.dtr = dmax_temp[date.day] - dmin_temp[date.day];

	// Progress report to user and update timer, first day of year only
	if (date.day == 0) {

		if (tmute.getprogress()>=1.0) {
			int years_to_simulate = nyear_spinup + last_hist_year - first_hist_year + 1;
			double progress=(double)(gridlist.getobj().id*years_to_simulate
				+date.year)/(double)(ngridcell*years_to_simulate);

			tprogress.setprogress(progress);
			dprintf("%3d%% complete, %s elapsed, %s remaining\n",(int)(progress*100.0),
				tprogress.elapsed.str,tprogress.remaining.str);
			tmute.settimer(MUTESEC);
		}
	}

	return true;
}


AsciiInput2::~AsciiInput2() {

	// Performs memory deallocation, closing of files or other "cleanup" functions.

	// Clean up
	gridlist.killall();
}
