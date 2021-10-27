.SECONDEXPANSION:
.SECONDARY:
.DELETE_ON_ERROR:

### OPTIONS
SECOND_PASS = 0

### SYSTEM TOOLS
LN = /bin/ln
MV = /bin/mv
CAT = /bin/cat
MKDIR = /bin/mkdir
ECHO = /bin/echo
CP = /bin/cp
CD = cd

### USER TOOLS
STAR = /home/paz2005/zeno/bin/src/STAR-STAR_2.4.2a/source/STAR
TRIM_GALORE = /zenodotus/abc/scratch/paz2005/bin/trim_galore
SAMTOOLS = /zenodotus/abc/scratch/paz2005/bin/samtools
FEATURE_COUNTS = /home/paz2005/zeno/bin/src/subread-1.4.6-p4-source/bin/featureCounts
JAVA = /scratchLocal/paz2005/pgs_pipeline/jdk1.8.0_171/bin/java
FASTQC = /home/paz2005/zeno/bin/src/fastqc_v0.11.3/FastQC/fastqc
QORTS = /zenodotus/abc/scratch/paz2005/bin/src/QoRTs_1.3/QoRTs.jar
RSCRIPT = /home/paz2005/zeno/bin/src/R-3.2.1/pzR3.2.1/bin/Rscript
R = /home/paz2005/zeno/bin/src/R-3.2.1/pzR3.2.1/bin/R

### RSCRIPTS
MERGE_GENE_COUNTS = /scratchLocal/paz2005/mojdeh_shakiba/rna_seq/scripts//mergeGeneCounts.R

### REFERENCES
REFERENCE = /zenodotus/abc/scratch/paz2005/references/GRCm38.p5/starIndex/
REFERENCE_FA = /zenodotus/abc/scratch/paz2005/references/GRCm38.p5/starIndex/GRCm38.primary_assembly.genome.fa
ANNOTATION = /zenodotus/abc/scratch/paz2005/references/GRCm38.p5/starIndex/gencode.vM16.annotation.gtf


### PARAMETERS
READ_LENGTH = 100
STAR_OPTIONS = --runThreadN 8 \
--runMode alignReads \
--genomeDir $(REFERENCE) \
--readFilesCommand zcat \
--outSAMstrandField intronMotif  \
--outFilterIntronMotifs RemoveNoncanonicalUnannotated  \
--outFilterType BySJout \
--outReadsUnmapped None \
--outSAMtype BAM Unsorted \
--chimOutType SeparateSAMold \
--limitSjdbInsertNsj 10000000 \


lc = $(subst A,a,$(subst B,b,$(subst C,c,$(subst D,d,$(subst E,e,$(subst F,f,$(subst G,g,$(subst H,h,$(subst I,i,$(subst J,j,$(subst K,k,$(subst L,l,$(subst M,m,$(subst N,n,$(subst O,o,$(subst P,p,$(subst Q,q,$(subst R,r,$(subst S,s,$(subst T,t,$(subst U,u,$(subst V,v,$(subst W,w,$(subst X,x,$(subst Y,y,$(subst Z,z,$1))))))))))))))))))))))))))


### DO NOT EDIT BELOW THIS LINE UNLESS YOU KNOW WHAT YOU ARE DOING
FASTQFILES := $(wildcard *_R1*.fastq.gz)
#TRIMMED_FASTQFILES :=  $(wildcard *_R1_val_1.fastq.gz)
TRIMMED_FASTQFILES := $(patsubst %.fastq.gz,%_val_1.fastq.gz,$(wildcard *_R1_???.fastq.gz))
SAMPLES := $(sort $(foreach a,$(FASTQFILES),$(firstword $(subst _, ,$a))))
FLOWCELLS := $(sort $(join $(foreach a,$(FASTQFILES),$(firstword $(subst _, ,$a))), $(addprefix _, $(foreach a,$(FASTQFILES),$(word 2, $(subst _, ,$a))))))
BARCODES :=  $(sort $(join $(join $(foreach a,$(FASTQFILES),$(firstword $(subst _, ,$a))), $(addprefix _, $(foreach a,$(FASTQFILES),$(word 2, $(subst _, ,$a))))), $(addprefix _, $(foreach a,$(FASTQFILES),$(word 3, $(subst _, ,$a))))))
LANES := $(sort $(join $(join $(join $(foreach a,$(FASTQFILES),$(firstword $(subst _, ,$a))), $(addprefix _, $(foreach a,$(FASTQFILES),$(word 2, $(subst _, ,$a))))), $(addprefix _, $(foreach a,$(FASTQFILES),$(word 3, $(subst _, ,$a))))), $(addprefix _, $(foreach a,$(FASTQFILES),$(word 4, $(subst _, ,$a))))))


all: gene_counts qc_data 
fastqc: $(patsubst %.fastq.gz,%_fastqc.zip,$(wildcard *_R1_???.fastq.gz))
first_pass: $(FASTQFILES:.fastq.gz=.1p.Aligned.out.bam) $(FASTQFILES:.fastq.gz=.1p.SJ.out.tab)
second_pass: $(FASTQFILES:.fastq.gz=.2p.Aligned.out.bam) $(FASTQFILES:.fastq.gz=.2p.SJ.out.tab)
sorted_bam: $(FASTQFILES:.fastq.gz=.sorted.bam)
merged_bam: $(addsuffix .merged.bam, $(SAMPLES))
coordSorted_bam: $(addsuffix .coordSorted.bam, $(SAMPLES))
gene_counts: $(addsuffix .gene.counts, $(SAMPLES)) gene.counts.txt
qc_data: $(addprefix ./qc/, $(addsuffix /QC.summary.txt, $(SAMPLES)))
decoder: decoder.txt

test:
	@echo $(SAMPLES)

### MERGE FASTQ FILES BY SAMPLE NAME
find-fastq-files = $(sort $(filter $1_% , $(FASTQFILES)))

define merge-fastq-files
$1.fastq.gz: $(call find-fastq-files,$1)
	$$(if $$(findstring 1, $$(words $$^)),$$(LN) -fs $$^ $$@,$(CAT) $$^ >> $$@)
endef

$(foreach s,$(LANES),$(eval $(call merge-fastq-files,$s)))

	
### FASTQC 
%_fastqc.zip: %.fastq.gz $$(subst R1,R2,%.fastq.gz)
	$(FASTQC) $^

%_fastqc.zip: %.fastq.gz
	$(FASTQC) $<

### FIRST PASS ALIGNMENT
%.1p.Aligned.out.bam %.1p.SJ.out.tab: %.fastq.gz $$(subst R1,R2,%.fastq.gz)
	$(STAR) --genomeLoad LoadAndKeep --readFilesIn $^ --outFileNamePrefix $*.1p. $(STAR_OPTIONS)

%.1p.Aligned.out.bam %.1p.SJ.out.tab: %.fastq.gz
	$(STAR) --genomeLoad LoadAndKeep --readFilesIn $< --outFileNamePrefix $*.1p. $(STAR_OPTIONS)
	
### SECOND PASS ALIGNMENT
%.2p.Aligned.out.bam %.2p.SJ.out.tab: %.fastq.gz $$(subst R1,R2,%.fastq.gz) $(FASTQFILES:.fastq.gz=.1p.SJ.out.tab)
	$(STAR) --genomeLoad NoSharedMemory --readFilesIn $< $(word 2,$^) --outFileNamePrefix $*.2p. --sjdbFileChrStartEnd $(TRIMMED_FASTQFILES:.fastq.gz=.1p.SJ.out.tab) $(STAR_OPTIONS)

%.2p.Aligned.out.bam %.2p.SJ.out.tab: %.fastq.gz
	$(STAR) --genomeLoad NoSharedMemory --readFilesIn $< --outFileNamePrefix $*.2p. --sjdbFileChrStartEnd $(TRIMMED_FASTQFILES:.fastq.gz=.1p.SJ.out.tab) $(STAR_OPTIONS)

### SORT BAM FILES
ifdef SECOND_PASS
ifeq ($(SECOND_PASS), 1)
%.sorted.bam: %.2p.Aligned.out.bam
	$(SAMTOOLS) sort -n -O bam -T $<.tmp $< > $@
else
%.sorted.bam: %.1p.Aligned.out.bam
	$(SAMTOOLS) sort -n -O bam -T $<.tmp $< > $@
endif
endif

### MERGE BAM FILES BY SAMPLE NAME
find-bam-files = $(sort $(filter $1_% , $(FASTQFILES:.fastq.gz=.sorted.bam)))

define merge-bam-files
$1.merged.bam: $(call find-bam-files,$1)
	$$(if $$(filter 1, $$(words $$^)),$$(LN) -fs $$^ $$@,$(SAMTOOLS) merge -f $$@ $$^)
endef

$(foreach s,$(SAMPLES),$(eval $(call merge-bam-files,$s)))

### COORD SORT
%.coordSorted.bam: %.merged.bam
	$(SAMTOOLS) sort -O bam -T $<.tmp $< > $@ ; $(SAMTOOLS) index $@

### COUNT READS PER GENE
%.gene.counts : %.merged.bam
	file_arr=(./$**R2*.fastq.gz) && if [ -f $${file_arr[0]} ] ; then $(FEATURE_COUNTS) -p -g "gene_name" -a $(ANNOTATION) -o $@  $< ; else $(FEATURE_COUNTS) -g "gene_name" -a $(ANNOTATION) -o $@  $< ; fi

### MERGE GENE COUNTS
gene.counts.txt: $(addsuffix .gene.counts, $(SAMPLES))
	$(RSCRIPT) $(MERGE_GENE_COUNTS) $@


### QC
./qc/%/QC.summary.txt: %.merged.bam
	file_arr=(./$*_*R2*.fastq.gz) && $(MKDIR) -p ./qc/$* && raw=$$(echo $$(zcat $*_*R1*.fastq.gz | wc -l)/4 | bc) && if [ -f $${file_arr[0]}  ] ; then $(JAVA) -Xmx70G -jar $(QORTS) QC --seqReadCt $$(echo $$raw) $< $(ANNOTATION) ./qc/$* ; else $(JAVA) -Xmx70G -jar $(QORTS) QC  --singleEnded --seqReadCt $$(echo $$raw) $< $(ANNOTATION) ./qc/$* ; fi
