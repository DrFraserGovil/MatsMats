#include <iostream>

#include "settings/Settings.h"

#include "tools/tools.h"
#include "parallel/parallel.h"
#include "biology/DNASequence.h"

int main(int argc, char**argv)
{
	Settings.Initialise(argc,argv);
	LOG(WARN) << sizeof(dnabits);
	return 0;
}