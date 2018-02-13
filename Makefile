SHELL:=/bin/bash

install:
	mkdir -p bin
	rm -f bin/*
	cp src/*.py src/SVXplorer src/VERSION bin

test:
	./bin/SVXplorer data/bams/test.discordants.bam data/bams/test.splitters.bam \
	data/bams/test.bam -f -w testSVE
	if cmp -s testFiles/variants.vcf testSVE/results/variants.vcf
	then
	   echo "Test successful"
	else
	   echo "Test failed. Please make sure data folder was not changed."
	   echo "Report bugs to Kunal: kk7t@virginia.edu"
	fi
