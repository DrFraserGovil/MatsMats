#include "SettingGroups.h"


bool SystemSettings::Validate()
{
	if (Verbosity > 3)
	{
		Verbosity.SetValue(3,true);
	}
	
	GlobalLog::Config.Initialise(Verbosity,true);	
	
	return true;
}
