# this file holds all the commands to be used to run VLITE-Fast
# made watching tutorial #1 by Dr. Matthew Kerr
vncserver --list 
uptime
# they have a lot of uptime
ssh vlite-master
# controls the disks, NFS .. vlite-master leads the way
ssh vlite-difx(nn)
# vd nodes are worker ndoes
# XML broadcasting
/users/mkerr/mtk@ --> /home/vlite-master/mtk 
# where everything sits
# /data on vlite-nrl is /dev/sda1 RAID
# NFS : /home/vlite-master/mtk/data
# to merge candidates
python /home/vlite-master/vlite-fast/scripts/merge_candidates.py
# scripts directory
/home/vlite-master/vlite-fast/scripts
debug_launch
		Key to lauching the telescope is by running the messenger file by giving a configuration file
		messenger config file is 
		hostname interface gpu reader_port writer_port info_port baseband_dada_key fb_dada_key write_fb nbit
vliteantennas.in
		/home/vlite-master/vliteantennas.in
		File holds all the vlite-antennas <-> nodes mapping
		vlite-antenna_index VLA_antenna_number difx_node_mapping difx_iface_mapping clock_offset _pad LOFIBERLENGTH enable/disable
debug_lauch
		lauches a set of start_single
		see logs realtime, terminals open in real time 
start_single
		hostname iface gpu_dev_id reader_port writer_port info_port baseband_dada_key fb_dada_key write_fb nbit
logcheck
		/home/vlite-master/logs
