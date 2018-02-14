SHELL:=/bin/bash

.PHONY: test

install:
	mkdir -p bin
	rm -f bin/*
	cp src/*.py src/SVXplorer src/VERSION bin

test:
	./bin/SVXplorer data/bams/test.discordants.bam data/bams/test.splitters.bam \
	data/bams/test.bam -f -w test
	@cmp -s testFiles/variants.vcf test/results/variants.vcf \
	|| (echo -e "Test failed. Please make sure data folder was not changed, and open an issue on Github."; exit 1)
	echo "Test successful"
