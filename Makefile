SHELL:=/bin/bash

.PHONY: test

install:
	mkdir -p bin
	rm -f bin/*
	cp src/*.py src/SVXplorer src/VERSION bin
test:
	./bin/SVXplorer testCases/discordants.bam testCases/splitters.bam \
	testCases/sample.bam testCases/ref/10kbp.random.ref.fa -x -f -w test -d
	cat test/results/variants.bedpe | cut -f1-6 > test/results/variants.6col.bedpe
	cat testFiles/variants.bedpe | cut -f1-6 > test/results/variants.testFile.6col.bedpe
	@cmp -s test/results/variants.6col.bedpe test/results/variants.testFile.6col.bedpe \
	|| (echo -e "Test failed. Please make sure data folder was not changed, and open an issue on Github."; exit 1)
	echo "Test successful"
	rm test/results/variants.6col.bedpe test/results/variants.testFile.6col.bedpe
