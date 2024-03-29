///////////////////////////////////////////////////////////////////////////////////////
/// \file inputmodule.cpp
/// \brief Implemenation file for inputmodule.h
///
/// \author Joe Siltberg
/// $Date: 2013-07-17 09:22:52 +0200 (Mi, 17 Jul 2013) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "inputmodule.h"
#include "guess.h"

///////////////////////////////////////////////////////////////////////////////////////
/// InputModuleRegistry
///

InputModuleRegistry& InputModuleRegistry::get_instance() {
	static InputModuleRegistry instance;
	return instance;
}

void InputModuleRegistry::register_input_module(const char* name, 
                                                InputModuleCreator imc) {
	modules.insert(make_pair(std::string(name), imc));
}
            
InputModule* InputModuleRegistry::create_input_module(const char* name) const {
	std::map<std::string, InputModuleCreator>::const_iterator itr = modules.find(name);
	
	if (itr != modules.end()) {
		return (itr->second)();
	}
	else {
		fail("Couldn't find input module %s\n", name);
		return 0;
	}
}
