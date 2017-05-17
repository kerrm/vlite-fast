#ifndef __VLITE_XML_H__
#define __VLITE_XML_H__

#include <expat.h>

#define VLITE_XML_MAX_DEPTH		5
#define VLITE_XML_ELEMENT_SIZE		32
#define VLITE_XML_CHAR_SIZE		1024

typedef struct
{
	void *userData;

	XML_Parser parser;

	int level;
	char element[VLITE_XML_MAX_DEPTH][VLITE_XML_ELEMENT_SIZE];
	int error_count;
	char string[VLITE_XML_CHAR_SIZE];
	int antId;
	int eopId;
} VliteXML;

int parseVliteXML(void *userData, XML_StartElementHandler start, XML_EndElementHandler end, const char *message);

#endif
