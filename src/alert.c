#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <expat.h>
#include <difxmessage.h>
#include "vlite_xml.h"
#include "alert.h"

/* Function that replaces illegal XML string characters with the
 * official "entity" replacements.  
 *
 * Code stolen from difxsend.c from difxmessage
 *
 * Returns:
 *   Increse of string size on success, or -1 on error.
 */
int expandEntityReferences(char *dest, const char *src, int maxLength)
{
	int i, j;

	for(i = j = 0; src[i]; i++)
	{
		if(j >= maxLength-7)
		{
			dest[maxLength-1] = 0;

			return -1;
		}

		if(src[i] == '>')
		{
			strcpy(dest+j, "&gt;");
			j += 4;
		}
		else if(src[i] == '<')
		{
			strcpy(dest+j, "&lt;");
			j += 4;
		}
		else if(src[i] == '&')
		{
			strcpy(dest+j, "&amp;");
			j += 5;
		}
		else if(src[i] == '"')
		{
			strcpy(dest+j, "&quot;");
			j += 6;
		}
		else if(src[i] == '\'')
		{
			strcpy(dest+j, "&apos;");
			j += 6;
		}
		else if(src[i] < 32)	/* ascii chars < 32 are not allowed */
		{
			sprintf(dest+j, "[[%3d]]", src[i]);
			j += 7;
		}
		else
		{
			dest[j] = src[i];
			j++;
		}
	}

	dest[j] = 0;

	return j - i;
}

int alertReceiveOpen(const CommandLineOptions *opts)
{
	return openMultiCastSocket(opts->alertGroup, opts->alertPort);
}

int alertReceiveClose(int sock)
{
	return closeMultiCastSocket(sock);
}

static void XMLCALL startElement(void *userData, const char *name, const char **atts)
{
	AlertDocument *A;
	VliteXML *X;

	X = (VliteXML *)userData;
	A = (AlertDocument *)X->userData;
	
	if(strcmp(name, "EVLAMessage") == 0)
	{
		int i;

		for(i = 0; atts[i]; i += 2)
		{
			if(strcmp(atts[i], "location") == 0)
			{
				snprintf(A->location, VLA_LOCATION_SIZE, "%s", atts[i+1]);
			}
			else if(strcmp(atts[i], "timestamp") == 0)
			{
				A->timeStamp = atof(atts[i+1]);	
			}
		}
	}
	else if(strcmp(name, "device") == 0)
	{
		int i;

		for(i = 0; atts[i]; i += 2)
		{
			if(strcmp(atts[i], "name") == 0)
			{
				snprintf(A->deviceName, VLA_DEVICE_SIZE, "%s", atts[i+1]);
			}
		}
	}
	else if(strcmp(name, "monitor") == 0)
	{
		int i;

		for(i = 0; atts[i]; i += 2)
		{
			if(strcmp(atts[i], "name") == 0)
			{
				snprintf(A->monitorName, VLA_MONITOR_NAME_SIZE, "%s", atts[i+1]);
			}
			else if(strcmp(atts[i], "alert") == 0)
			{
				A->alertState = atoi(atts[i+1]);
			}
		}
	}
}

static void XMLCALL endElement(void *userData, const char *name)
{
}

int parseAlertDocument(AlertDocument *A, const char *message)
{
	int v;

	memset(A, 0, sizeof(AlertDocument));

	v = parseVliteXML(A, startElement, endElement, message);

	if(A->location[0] == 'e' && A->location[1] == 'a' && isdigit(A->location[2]) && isdigit(A->location[3]) && A->location[4] == '-')
	{
		A->locationSuffix = A->location + 5;
		A->vlaAnt = 10*(A->location[2]-'0') + (A->location[3]-'0');
	}

	return v;
}

void printAlertDocument(const AlertDocument *A)
{
	printf("Alert Document:\n");
	printf("  location = %s\n", A->location);
	printf("  vlaAnt = %d\n", A->vlaAnt);
	printf("  locationSuffix = %s\n", A->locationSuffix);
	printf("  device = %s\n", A->deviceName);
	printf("  monitorName = %s\n", A->monitorName);
	printf("  timeStamp = %14.8f\n", A->timeStamp);
	printf("  alertState = %d\n", A->alertState);
}

void fprintAlertDocument(const AlertDocument *A, FILE *fd)
{
  fprintf(fd,"Alert Document:\n");
  fprintf(fd,"  location = %s\n", A->location);
  fprintf(fd,"  vlaAnt = %d\n", A->vlaAnt);
  fprintf(fd,"  locationSuffix = %s\n", A->locationSuffix);
  fprintf(fd,"  device = %s\n", A->deviceName);
  fprintf(fd,"  monitorName = %s\n", A->monitorName);
  fprintf(fd,"  timeStamp = %14.8f\n", A->timeStamp);
  fprintf(fd,"  alertState = %d\n", A->alertState);
}
