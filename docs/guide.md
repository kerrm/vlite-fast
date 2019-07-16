# The super fast way to dive into VLITE-Fast


## TL;DR
		If you know what you're doing, do this:
		``debug_launch``

## Theory
		`debug_launch` should be run on `vlite-nrl`. 
		
		It all starts with `messenger` which :
		- connects to every vlite-difx nodes on {reader_port|writer_port|info_port} control.
		- connects to VLA {obsinfo|antprop|alert} multicast and captures XML (AntPropDocument|ScanInfoDocument|AlertDocument)
		`messenger` is like the leader in that it _controls_ based on the captured XMLs which are multicast. 

		`debug_launch` also provides the `start_single` commands which are to pasted in every vlite-difx* nodes and run there.
		Each `start_single` does :
		- `start_writer`
				- `start_dada`
						Yet another script is called.
				- `writer`
				- `dada_db` ???
		- `start_process`
				- `start_dada`
						Yet another script is called.
				- `nvidia-smi` _gpu logging_
				- `process_baseband`
				- `dada_db` ???

		`writer` is pretty important C code, hence, dedicating a paragraph to it. It receives `baseband_dada_key`.
		`writer` connects to PSRDADA and sets up socket to receive a ScanInfoDocument (ObservationDocument). 
		It collects raw packets and writes the raw voltages to PSRDADA.
		It works with raw voltage data so it handles dumps after receiving triggers too.

		`process_baseband` is important GPU+C code which as the name suggests performs
		baseband --> filterbank processing. 
		It's input is PSRDADA `shmkey`. Output can be also be PSRDADA but is not default behavior. 
		Some straightforward CLI arguments which I am skipping.
		It does:
				- Temporal kurtosis
				- Bandpass normalization
				- Add Polarizations
				- Average in time domain
				- Zap channels and Digitize
		- Histogram|Kurto (if asked)
		- Launches heimdall
		This happens on a per-second data
		
What is VDIF Frame Number?? What is VDIFThreadID?? 

## Glossary

### debug_launch
		`hostname interface gpu_id reader_port writer_port info_port baseband_dada_key fb_dada_key write_fb nbit`

		_ssh into every `vlite-difx` nodes and run `start_single` with arguments._
		_runs messenger with configuration._

### start_single
		`hostname interface gpu_id reader_port writer_port info_port baseband_dada_key fb_dada_key write_fb nbit`
		
		_`start_writer`_
				`interface writer_port info_port baseband_dada_key`
		_`start_process`_
				`reader_port baseband_dada_key fb_dada_key write_fb nbit`

### start_writer
		`interface writer_port info_port baseband_dada_key`

		_`start_dada`_
				`baseband_dada_key 257638400 32`
		_`writer`_
				`interface writer_port info_port baseband_dada_key`
				Binary
		_`dada_db`_
				`baseband_dada_key`
				Binary (is part of psrdada)

### start_process
		`reader_port baseband_dada_key fb_dada_key write_fb nbit gpu_id`	

		_`start_dada`_
				If input `fb_dada_key` is greater than 0.
				Why?
		_`nvidia-smi`_
				Logging if not already present
				Binary (is part of nVidia)
		_`process_baseband`_
				`reader_port baseband_dada_key fb_dada_key write_fb nbit gpu_id`	
				Binary
		_`dada_db`
				`fb_dada_key` as daemon
				Binary (is part of psrdada)

### start_dada
		`dada_key size_in_bytes number_of_buffers`

		_`dada_db`
				`fb_dada_key` as daemon
				Binary (is part of psrdada)








