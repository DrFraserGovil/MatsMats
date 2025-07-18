#pragma once
#include "HelpMessages.h"
#define SETTINGS_CATEGORY SystemSettings
#define SETTINGS_FILE "definitions/system.def"
#define SETTINGS_VALIDATE
#include "SettingsConstructor.h"




//this here defines the global set which will be inserted into the application settings using a familiar X-macro pattern
//why not put into its own file? a) overkill -- parameters are far more likely to change than entire settings, and b) it's very obvious when a setting is missing 
//(i.e. it won't compile, because you can't access Settings.NonExistantGroup) 
#define SETTINGS_GROUPS \
	S_GROUP(SystemSettings, System) \
	

// SGroupDefinition.h
#ifndef S_GROUP_DEFINITION_H
#define S_GROUP_DEFINITION_H

#define S_GROUP(type,name) type name;

#endif // S_GROUP_DEFINITION_H