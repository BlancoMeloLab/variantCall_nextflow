#!/bin/sh

module load bedtools samtools igvtools python/2.7.14/
WD=/Users/qiya9811/DTRA/ext_ADAR-KO_IFN_A549/mapped
mkdir $WD/bed $WD/tdf
for name in $(ls $WD/bams/*.bam)
do
  fn=$(basename -- $name)
  fn="${fn%.*.*}"
  samtools flagstat $name > $name.flagstat
  samtools view -@ 32 -h -b -f 0x0040 ${name} > ${name}.R1
  samtools view -@ 32 -h -b -f 0x0080 ${name} > ${name}.R2
  genome=hg38
  genomeCoverageBed -bg -split -strand - -ibam ${name}.R1 -g $genome | sortBed -i - > $WD/bed/${fn}.r1.pos
  genomeCoverageBed -bg -split -strand + -ibam ${name}.R1 -g $genome | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }' | sortBed -i - > $WD/bed/${fn}.r1.neg
  genomeCoverageBed -bg -split -strand + -ibam ${name}.R2 -g $genome > $WD/bed/${fn}.r2.pos
  genomeCoverageBed -bg -split -strand - -ibam ${name}.R2 -g $genome | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }' | sortBed -i - > $WD/bed/${fn}.r2.neg
  unionBedGraphs -i $WD/bed/${fn}.r1.pos $WD/bed/${fn}.r2.pos | awk -F '\t' '{OFS="\t"; print $1,$2,$3,$4+$5;}' - > $WD/bed/${fn}.pos.BedGraph
  unionBedGraphs -i $WD/bed/${fn}.r1.neg $WD/bed/${fn}.r2.neg | awk -F '\t' '{OFS="\t"; print $1,$2,$3,$4+$5;}' - > $WD/bed/${fn}.neg.BedGraph
  cat $WD/bed/${fn}.pos.BedGraph $WD/bed/${fn}.neg.BedGraph | sortBed -i - > $WD/bed/${fn}.BedGraph
  python $WD/scripts/readcountcorrectBGv2.py $WD/bed/${fn}.BedGraph $name.flagstat $WD/bed/${fn}.rcc.BedGraph
  igvtools toTDF $WD/bed/${fn}.rcc.BedGraph $WD/tdf/${fn}.tdf /scratch/Users/qiya9811/opt/IGVTools/genomes/hg38.chrom.sizes
done
