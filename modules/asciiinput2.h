///////////////////////////////////////////////////////////////////////////////////////
/// \file demoinput.h
/// \brief Input module for demo data set
///
/// \author Joe Siltberg
/// $Date: 2015-03-20 11:18:13 +0100 (Fri, 20 Mar 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_ASCIIINPUT2_H
#define LPJ_GUESS_ASCIIINPUT2_H

#include "guess.h"
#include "guessstring.h"
#include "cruinput.h"

#include <fstream>
#include <sstream>
#include "inputmodule.h"
#include <vector>
#include "gutil.h"

class AsciiInput2 : public InputModule {
public:

	AsciiInput2();

	~AsciiInput2();

	void init();

	bool getgridcell(Gridcell& gridcell);

	bool getclimate(Gridcell& gridcell);

	void getlandcover(Gridcell& gridcell);

private:

	struct Coord {

		// Type for storing grid cell longitude, latitude and description text

		int id;
		double lon;
		double lat;
		std::string descrip;
	} ;

	/// Landcover fractions read from ins-file (% area).
	/** One entry for each land cover type */
	std::vector<int> lc_fixed_frac;

	/// Whether gridcell is divided into equal active landcover fractions.
	bool equal_landcover_area;

	/// A list of Coord objects containing coordinates of the grid cells to simulate
	ListArray_id<Coord> gridlist;

	/// The number of grid cells to simulate
	int ngridcell;

	// Timers for keeping track of progress through the simulation
	Timer tprogress,tmute;
	static const int MUTESEC=20; // minimum number of sec to wait between progress messages

	/// Yearly CO2 data read from file
	/**
	 * This object is indexed with calendar years, so to get co2 value for
	 * year 1990, use co2[1990]. See documentation for GlobalCO2File for
	 * more information.
	 */
	GlobalCO2File co2;

	/// Path to CRU binary archive
	xtring file_cru;

	/// Temperature for current gridcell and current year (deg C)
	double dtemp[365];

	/// Minimum temperature for current gridcell and current year (deg C)
	double dmin_temp[365];

	/// Maximum temperature for current gridcell and current year (deg C)
	double dmax_temp[365];

	/// Precipitation for current gridcell and current year (mm/day)
	double* dprec;

	/// Insolation for current gridcell and current year (\see instype)
	double dinsol[365];

	/// Daily N deposition for one year
	double dndep[365];

	/// /// Mean PPFD (micro mol m-2 s-1, 24h based) for this month for use in dark respiration response. FL 141107
	double mppfd[12];

	int year_in_file;

	/// Nitrogen deposition forcing for current gridcell
	Lamarque::NDepData ndep;

	/// Nitrogen deposition time series to use (historic,rcp26,...)
	std::string ndep_timeseries;

	void site_specific_met(double lon, double lat);
	std::vector<double> vrad, vprec, vtemp;
	int nyear_hist;

};

#endif // LPJ_GUESS_DEMOINPUT_H
