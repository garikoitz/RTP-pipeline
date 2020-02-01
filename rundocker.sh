dataDir="/black/localhome/glerma/TESTDATA/RTP-pipeline"

#docker run -ti --rm --entrypoint /bin/bash  \
docker run -ti --rm  $2 \
	       -v $dataDir/input:/flywheel/v0/input  \
   	       -v $dataDir/output:/flywheel/v0/output  \
           -v $(pwd)/example_config.json:/flywheel/v0/config.json \
   	       garikoitz/rtp-pipeline:$1
