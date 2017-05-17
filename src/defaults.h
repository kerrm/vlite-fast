#ifndef __DEFAULTS_H__
#define __DEFAULTS_H__

/* Some VLA specific values */
#define VLA_ANTENNA_COUNT		28		/* how many antennas hath the VLA */
#define VLA_PROJECT_CODE_SIZE		32
#define VLA_CENTER_X			-1601185.4	/* [m] in VLBI coordinate system (same as used by VLA Executor 20140207) */
#define VLA_CENTER_Y			-5041977.5	/* [m] in VLBI coordinate system (same as used by VLA Executor 20140207) */
#define VLA_CENTER_Z			3554875.9	/* [m] in VLBI coordinate system (same as used by VLA Executor 20140207) */

/* Some VLITE specific values */
#define VLITE_ANTENNA_COUNT		10		/* how many antennas in VLITE */
#define VLITE_MIN_SUBARRAY		3		/* minimum number of antennas in a VLITE subarray */
#define VLITE_TEST_SUBARRAY_SIZE	2		/* minimum number of antennas in a test subarray */
#define VLITE_SUBARRAYS			(VLITE_ANTENNA_COUNT/VLITE_MIN_SUBARRAY)	/* maximum number of subarrays to support */
#define STATE_N_EOP			5
#define VLITE_MAX_DIFX_CORE		20

/* Default vlite system parameters */
#define VLITE_CONFIG_DIRECTORY		"../data"	/* where to find the system files for VLITE */
#define DEFAULT_STRIP_BYTES		42		/* number of bytes to lop of start of each VDIF packet */
	/* note: the four directories below must live on the same filesystem */
#define DEFAULT_PRIVATE_DIRECTORY	"./private"	/* place where difx files live and FITS files are assembled */
#define DEFAULT_STAGING_DIRECTORY	"./staging"	/* place to stuff complete FITS files for VLITE backend to find */
#define DEFAULT_COPYING_DIRECTORY	"./copying"	/* temporary place where VLITE backend works on the files */
#define DEFAULT_ARCHIVE_DIRECTORY	"./archive"	/* endpoint for FITS files; VLITE backend moves files here when it is done */
#define MAX_TEMP_FILE_AGE		600		/* [s] max age of temporary files before raising alert */
#define TEMP_TCAL_FILE			"vlite.tcal.tmp"
#define TEMP_CONFIG_FILE		"vlite.config.tmp"
#define TEMP_NODES_FILE			"vlite.nodes.tmp"
#define MIN_FREE_SPACE_MB		100000		/* [MB] An alert will be issued if remaining space drops below this amount */
#define MAX_START_DELAY			10		/* [s] max time to wait after trying to start job for subarray to be idle */
#define STATE_LOOKBACK			10000		/* how many jobs of history to keep */

/* Default difx parameters */
#define DEFAULT_SCAN_LENGTH		3600		/* [s] maximum length of a VLITE correlation */
#define DEFAULT_VIS_BUFFER_LENGTH	40
#define DEFAULT_DATA_BUFFER_FACTOR	32
#define DEFAULT_N_DATA_SEGMENTS		8
#define DEFAULT_DELAY_POLY_ORDER	5		/* 5 is usual for VLBI */
#define DEFAULT_DELAY_POLY_INTERVAL	600		/* [s] ; 120 is usual for VLBI */
#define DEFAULT_INTEGRATION_TIME	2		/* [s] ; 2 is the VLITE spec */
#define DEFAULT_SUBINT_NS		10000000
#define DEFAULT_GUARD_NS		400
#define DEFAULT_N_CHAN			640		/* 640 is the VLITE spec for 100 kHz channels */
#define DEFAULT_N_FFT_CHAN		640		/* minimum consistent with spec; double or quadruple may be better */

/* Some string length values */
#define VLITE_DIRECTORY_SIZE		256		/* length of a unix directory */
#define VLITE_FILENAME_SIZE		256		/* length of a unix file */
#define DATABASE_ERROR_MESSAGE_SIZE	512
#define EXECUTOR_DATASETID_SIZE		64		/* size in bytes of the executor dataset identifier (used to connect an observation document to an antenna properties document */
#define DIFXNODE_SIZE			64
#define ETHDEV_SIZE			32
#define COORD_STRING_SIZE		24
#define SCAN_INTENT_SIZE		512

/* Some network parameters */
#define MULTICAST_GROUP_SIZE		16		/* max lenth + 1 of group strings */
#define EVLA_ALERT_PORT			20011
#define EVLA_ALERT_GROUP		"239.192.2.3"
#define EVLA_SUBARRAY_PORT		53000
#define EVLA_SUBARRAY_GROUP		"239.192.3.1"
#define EVLA_SCANINFO_PORT		53001
#define EVLA_SCANINFO_GROUP		"239.192.3.2"

#define TEST_ALERT_PORT			20911
#define TEST_ALERT_GROUP		"239.199.2.3"
#define TEST_SUBARRAY_PORT		53900
#define TEST_SUBARRAY_GROUP		"239.199.3.1"
#define TEST_SCANINFO_PORT		53901
#define TEST_SCANINFO_GROUP		"239.199.3.2"

#define D306_SWITCHEDPOWER_PORT		54321

#define SERVICE_PORT			54322
#define MAX_SERVICE_CONNECTIONS		10

/* Timnig constants */
#define MJD2000				51544.5		/* 2000.0 occured at noon */
#define MJD_UNIX0			40587.0		/* MJD at beginning of unix time */

#endif
