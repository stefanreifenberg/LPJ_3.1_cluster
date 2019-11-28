///////////////////////////////////////////////////////////////////////////////////////
/// \file commandlinearguments.h
/// \brief Takes care of the command line arguments to LPJ-GUESS
///
/// $Date: 2013-07-17 09:22:52 +0200 (Mi, 17 Jul 2013) $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_COMMAND_LINE_ARGUMENTS_H
#define LPJ_GUESS_COMMAND_LINE_ARGUMENTS_H

#include <string>

/// Parses and stores the user's command line arguments
class CommandLineArguments {
public:
	/// Send in the arguments from main() to this constructor
	/** If the arguments are malformed, this constructor will
	 *  print out usage information and exit the program.
	 */
	CommandLineArguments(int argc, char** argv);

	/// Returns the instruction filename
	const char* get_instruction_file() const;

	/// Returns true if the user has specified the help option
	bool get_help() const;

	/// Returns true if the user has specified the parallel option
	bool get_parallel() const;

	/// Returns the chosen (or default) input module
	const char* get_input_module() const;

private:
	/// Does the actual parsing of the arguments
	bool parse_arguments(int argc, char** argv);

	/// Prints out usage information and exits the program
	void print_usage(const char* command_name) const;

	/// Instruction filename specified on the command line
	std::string insfile;

	/// Whether the user wants help on how to run the LPJ-GUESS command
	bool help;

	/// Whether the user requested a parallel run
	bool parallel;

	/// The chosen (or default) input module
	std::string input_module;
};

#endif // LPJ_GUESS_COMMAND_LINE_ARGUMENTS_H
