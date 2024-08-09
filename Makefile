EXAMPLESAM       :=  "data/210610_50 against HDV 1.sam"
EXAMPLEREFERENCE :=  data/substituted_reference.fa
BOWTIEDIR        :=  bowtie_dir
FASTQDIR         :=  output/fastq
FASTQFILES       :=  $(wildcard $(FASTQDIR)/*.fq)

all: clean_output dedup

clean_output:
	touch output/justafilesoidontgetanerroriftherearenofiles
	rm -r output/*

dedup:
	./sort_reads.py $(EXAMPLESAM) $(EXAMPLEREFERENCE) -t output

clean_bowtie_tmp:
	rm -rf /tmp/bt2-*

bt-clean:
	rm -rf $(BOWTIEDIR)

bt-index:
	mkdir -p $(BOWTIEDIR)
	bowtie2-build -c -f $(EXAMPLEREFERENCE) $(BOWTIEDIR)/index

bt-align:
	for fq in $(FASTQFILES); do \
		cat $$fq >> $(BOWTIEDIR)/all_sequences.fq ;\
	done
	bowtie2 -q --xeq --local --very-sensitive-local -N 1 --no-unal -x "$(BOWTIEDIR)/index" -U $(BOWTIEDIR)/all_sequences.fq > '$(BOWTIEDIR)/result.sam' ;\

bt-all: bt-clean bt-index bt-align
	samtools view $(BOWTIEDIR)/result.sam | less
