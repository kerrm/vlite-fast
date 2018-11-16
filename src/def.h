/*
#define CMD_START 101
#define CMD_STOP 102
#define CMD_QUIT 103
#define CMD_EVENT 104
#define CMD_NONE 0
*/

#define CMD_START 'S'
#define CMD_STOP  'C'
#define CMD_QUIT  'Q'
#define CMD_EVENT 'E'
#define CMD_NONE  'N'

#define STATE_STARTED 201
#define STATE_STOPPED 202

#define VDIF_PKT_SIZE 5032
#define UDP_HDR_SIZE 42

#define SRCMAXSIZE 32

#define ETHDEV "eth0"

#define MAXFRAMENUM 25599

#define EVENTDIR "/mnt/ssd/dumps"
#define LOGDIR "/data/kerrm/test_vlite/logs"
#define DATADIR "/data/kerrm/test_vlite/data"
#define OBSINFODIR "/data/kerrm/test_vlite/logs/obsinfo-antprop"

