#include <string.h>
#include "vlite_xml.h"

static void XMLCALL charHandler(void *userData, const XML_Char *str, int len)
{
	VliteXML *X;
	int l;

	X = (VliteXML *)userData;

	l = strlen(X->string);

	if(len + l >= VLITE_XML_CHAR_SIZE)
	{
		len = VLITE_XML_CHAR_SIZE - l - 1;
		if(len <= 0)
		{
			return;
		}
	}
	strncpy(X->string+l, str, len);
	X->string[len+l] = 0;
}

int parseVliteXML(void *userData, XML_StartElementHandler start, XML_EndElementHandler end, const char *message)
{
	VliteXML X;

	memset(&X, 0, sizeof(VliteXML));

	X.userData = userData;
	X.level = -1;
	X.parser = XML_ParserCreate(0);
	
	XML_ParserReset(X.parser, 0);
	XML_SetElementHandler(X.parser, start, end);
	XML_SetCharacterDataHandler(X.parser, charHandler);
	XML_SetUserData(X.parser, &X);
	XML_Parse(X.parser, message, strlen(message), 0);
	XML_ParserFree(X.parser);

	return X.error_count;
}
