#ifndef __DEF_H__
#define __DEF_H__

#define CMD_START 'S'
#define CMD_STOP  'C'
#define CMD_QUIT  'Q'
#define CMD_EVENT 'E'
#define CMD_NONE  'N'
#define CMD_FAKE_START  'F'
#define CMD_FAKE_STOP  'G'

#define STATE_STARTED 201
#define STATE_STOPPED 202

#define VDIF_FRAME_SIZE 5032
#define UDP_HDR_SIZE 42

#define SRCMAXSIZE 32

#define ETHDEV "eth0"

#define MAXFRAMENUM 25599
#define FRAMESPERSEC 25600

#define EVENTDIR "/mnt/ssd/dumps"
#define LOGDIR "/home/vlite-master/mtk/logs"
//#define DATADIR "/home/vlite-master/mtk/data"
#define DATADIR "/mnt/ssd/fildata"
#define CANDDIR "/mnt/ssd/cands"
#define OBSINFODIR "/home/vlite-master/mtk/logs/obsinfo-antprop"


#define MAX_THREADIOS 100
#define MAX_TRIGGERS 10

#endif
