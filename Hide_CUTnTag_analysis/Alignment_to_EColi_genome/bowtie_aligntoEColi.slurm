#! /bin/bash
#SBATCH --verbose
#SBATCH -J bowtie2_EColi
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

# sbatch -p risc --array=1-4 --export=outdir=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/Hide_CUTnTag_analysis/Alignment_to_EColi_genome,inputfiles_R1=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/Hide_CUTnTag_analysis/Alignment_to_EColi_genome/inputfiles_R1.txt,inputfiles_R2=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/Hide_CUTnTag_analysis/Alignment_to_EColi_genome/inputfiles_R2.txt,genomeIndexName=/rugpfs/fs0/risc_lab/store/risc_data/downloaded/ecoli_mg1655/genome/Sequence/Bowtie2Index/genome,samplelist=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/Hide_CUTnTag_analysis/Alignment_to_EColi_genome/samplenames.txt,inputbed=/lustre/fs4/risc_lab/scratch/jyeung/for_Hide/Hide_CUTnTag_analysis/inputbed.txt,chromSize=Xenopus_laevis_v10.1.chrom.size /lustre/fs4/risc_lab/scratch/jyeung/for_Hide/Hide_CUTnTag_analysis/Alignment_to_EColi_genome/bowtie_aligntoEColi.slurm
source activate fastq2bam
INPUT_R1=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_R1 | tail -n 1)
INPUT_R2=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_R2 | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID $samplelist | tail -n 1)
INPUTBED=$(head -n $SLURM_ARRAY_TASK_ID $inputbed | tail -n 1)
OUTPUT=$outdir\/$NAME

bowtie2 --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 8 -x $genomeIndexName -1 $INPUT_R1 -2 $INPUT_R2 -S $OUTPUT\_EColi_spikeIn.sam &> $OUTPUT\_bowtie2_summary_EColi_spikeIn.txt

# 0x04 corresponds to the "unmapped" flag in the SAM/BAM format. By using -F 0x04, you're telling Samtools to exclude any reads from the output where this flag is set, effectively filtering out reads that are not properly aligned to the reference genome.
seqDepthDouble=`samtools view -F 0x04 $OUTPUT\_EColi_spikeIn.sam | wc -l`
# save number of read pairs aligned to Ecoli genome to .seqDepth file. 
seqDepth=$((seqDepthDouble/2))
echo $seqDepth >$OUTPUT\_bowtie2_summary_EColi_spikeIn.seqDepth
    
# generate scaling factor and create bedgraph normalized to scaling factor based on number of E.Coli Spikein reads. 10000 is an arbitrary constant. 
    scale_factor=`echo "10000 / $seqDepth" | bc -l`
    echo "Scaling factor for $histName is: $scale_factor!"  >> $OUTPUT\_bowtie2_summary_EColi_spikeIn.seqDepth
    bedtools genomecov -bg -scale $scale_factor -i $INPUTBED -g $chromSize > $OUTPUT\_fragments_normalized_EColi_spikeIn.bedgraph
