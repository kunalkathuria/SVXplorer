SHELL:=/bin/bash

.PHONY: test

install:
	mkdir -p bin
	rm -f bin/*
	cp src/*.py src/SVXplorer src/VERSION bin
test:
	./bin/SVXplorer testCases/discordants.bam testCases/splitters.bam \
	testCases/sample.bam testCases/ref/10kbp.random.ref.fa -x -f -w test -d
	@cmp -s testFiles/variants.vcf test/results/variants.vcf \
	|| (echo -e "Test failed. Please make sure data folder was not changed, and open an issue on Github."; exit 1)
	echo "Test successful"
