#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defaults.h"
#include "options.h"

void resetCommandLineOptions(CommandLineOptions *opts)
{
	memset(opts, 0, sizeof(CommandLineOptions));
}

CommandLineOptions *newCommandLineOptions(int argc, char **argv)
{
	CommandLineOptions *opts;
	int a;

	opts = (CommandLineOptions *)malloc(sizeof(CommandLineOptions));

	if(!opts)
	{
		fprintf(stderr, "Error: newCommandLineOptions: cannot allocate %d bytes\n", (int)sizeof(CommandLineOptions));

		exit(EXIT_FAILURE);
	}

	/* trivial defaults */
	resetCommandLineOptions(opts);

	for(a = 1; a < argc; ++a)
	{
		if(argv[a][0] == '-')
		{
			if(strcmp(argv[a], "--verbose") == 0 ||
			   strcmp(argv[a], "-v") == 0)
			{
				++opts->verbose;
			}
			else if(strcmp(argv[a], "--quiet") == 0 ||
			   strcmp(argv[a], "-q") == 0)
			{
				--opts->verbose;
			}
			else if(strcmp(argv[a], "--help") == 0 ||
			   strcmp(argv[a], "-h") == 0)
			{
				++opts->usage;
			}
			else if(strcmp(argv[a], "--pband") == 0 ||
			   strcmp(argv[a], "-p") == 0)
			{
				opts->allowPband = 1;
			}
			else if(strcmp(argv[a], "--test") == 0 ||
			   strcmp(argv[a], "-t") == 0)
			{
				opts->testMode = 1;
			}
			else if(strcmp(argv[a], "--sniff") == 0 ||
			   strcmp(argv[a], "-s") == 0)
			{
				opts->doSniff = 1;
			}
			//else if(a+1 < argc)
			//{
			//	
			//}
			else
			{
				fprintf(stderr, "Unknown parameter %s\n", argv[a]);
				deleteCommandLineOptions(opts);

				return 0;
			}
		}
		else
		{
			fprintf(stderr, "Unexpected value %s\n", argv[a]);
			deleteCommandLineOptions(opts);

			return 0;
		}
	}

	/* non-trivial defaults */
	opts->program = argv[0];

	if(opts->testMode > 0)
	{
		/* This set is for test mode */

		opts->alertPort = TEST_ALERT_PORT;
		strcpy(opts->alertGroup, TEST_ALERT_GROUP);
		opts->scanInfoPort = TEST_SCANINFO_PORT;
		strcpy(opts->scanInfoGroup, TEST_SCANINFO_GROUP);
		opts->subarrayPort = TEST_SUBARRAY_PORT;
		strcpy(opts->subarrayGroup, TEST_SUBARRAY_GROUP);
	}
	else
	{
		/* This set corresponds to normal operation */

		opts->alertPort = EVLA_ALERT_PORT;
		strcpy(opts->alertGroup, EVLA_ALERT_GROUP);
		opts->scanInfoPort = EVLA_SCANINFO_PORT;
		strcpy(opts->scanInfoGroup, EVLA_SCANINFO_GROUP);
		opts->subarrayPort = EVLA_SUBARRAY_PORT;
		strcpy(opts->subarrayGroup, EVLA_SUBARRAY_GROUP);
	}

	return opts;
}

void deleteCommandLineOptions(CommandLineOptions *opts)
{
	if(opts)
	{
		free(opts);
	}
	else
	{
		fprintf(stderr, "Warning: deleteCommandLineOptions: null pointer\n");
	}
}

void printCommandLineOptions(const CommandLineOptions *opts)
{
	printf("Vlite Command Line Options\n\n");
	printf("  verbose         = %d\n", opts->verbose);
	printf("  usage           = %d\n", opts->usage);
	printf("  program         = %s\n", opts->program);
	printf("  scan info port  = %d\n", opts->scanInfoPort);
	printf("  scan info group = %s\n", opts->scanInfoGroup);
	printf("  alert port      = %d\n", opts->alertPort);
	printf("  alert group     = %s\n", opts->alertGroup);
	printf("  do sniff        = %d\n", opts->doSniff);
	printf("  allow P-band    = %d\n", opts->allowPband);
	printf("  test mode       = %d\n", opts->testMode);
	printf("\n\n");
}
