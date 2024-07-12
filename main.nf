#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input_bams = "/orfeo/cephfs/scratch/burlo/nardone/vcf_SV_valentina/ont_benchmark/BAM/downsampling/{minimap,winnowmap,lra}/HG002*.bam{,.bai}"
params.outdir = "/orfeo/cephfs/scratch/burlo/nardone/vcf_SV_valentina/ont_benchmark/SVCALL"

bams = Channel.fromFilePairs(${params.input_bams})

process NANOSV {
    tag {sample}_nanosv

    publishDir("${params.outdir}/nanosv", mode: 'copy')

    input:
        tuple val(sample), path(bam)
    output:
        tuple val(sample), path(${sample}.nanosv.vcf.gz)
    script:
        """
        NanoSV ${bam.first()} -o ${params.outdir}/nanosv/NanoSV.tmp.vcf -b /opt2/NanoSV/human_hg38.bed -t 8 && \
        bcftools view -i 'INFO/SVTYPE == "DEL"' ${params.outdir}/nanosv/NanoSV.tmp.vcf | bcftools sort -Oz -o ${params.outdir}/nanosv/${sample}.nanosv.vcf.gz && bcftools index ${params.outdir}/nanosv/${sample}.nanosv.vcf.gz
        """
}

process DUET {
    tag {sample}_duet

    publishDir("${params.outdir}/duet", mode: 'copy')

    input:
        tuple val(sample), path(bam)
    output:
        tuple val(sample), path(${sample}.duet.vcf.gz)
    script:
        """
        duet ${bam.first()} /ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ${params.outdir}/duet -t 8 &&
        sed '5 i\##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' ${params.outdir}/duet/phased_sv.vcf | bcftools view -i 'INFO/SVTYPE == "DEL"' | bcftools sort -Oz -o ${params.outdir}/duet/${sample}.duet.vcf.gz &&
        bcftools index ${params.outdir}/duet/${sample}.duet.vcf.gz
        """
}

process SNIFFLES {
    tag {sample}_sniffles

    publishDir("${params.outdir}/sniffles", mode: 'copy')

    input:
        tuple val(sample), path(bam)
    output:
        tuple val(sample), path(${sample}.sniffles.vcf.gz)
    script:
        """
        sniffles -t 8 --tandem-repeats /orfeo/cephfs/scratch/burlo/nardone/vcf_SV_valentina/ont_benchmark/resources/human_GRCh38_no_alt_analysis_set.trf.bed --input ${bam.first()} --output ${params.outdir}/sniffles/${sample}.sniffles.vcf.gz --reference /ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
        """
}

process NANOVAR {
    tag {sample}_nanovar

    publishDir("${params.outdir}/nanovar", mode: 'copy')

    input:
        tuple val(sample), path(bam)
    output:
        tuple val(sample), path(${sample}.nanovar.vcf.gz)
    script:
        """
        nanovar ${bam.first()} /ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ${params.outdir}/nanovar -t 8 -x ont -l 50 && 
        bcftools view --write-index -Oz -o ${params.outdir}/nanovar/${sample}.nanovar.vcf.gz ${params.outdir}/nanovar/${sample}.nanovar.pass.vcf
        """
}


process CAMPHOR {
    tag {sample}_camphor

    publishDir("${params.outdir}/camphor", mode: 'copy')

    input:
        tuple val(sample), path(bam)
    output:
        tuple val(sample), path(${sample}.camphor_DEL_only_sorted.vcf.gz)
    script:        
        """
        #!/usr/bin/env python3

        from glob import glob
        import re
        import os

        def couple_fastq(sample):
            cov_str = sample.split("_")[-1]
            cov_list = re.findall(r'\d+', cov_str)
            cov = "".join(cov_list)
            pattern = f"/orfeo/cephfs/scratch/burlo/nardone/vcf_SV_valentina/ont_benchmark/FASTQS/downsampling/HG002_minimap_down_to_{cov}X.fastq.gz"
            search = glob(pattern)
            fq = "".join(search)
            return fq

        def main():
            fq_dir = "/orfeo/cephfs/scratch/burlo/nardone/vcf_SV_valentina/ont_benchmark/FASTQS/downsampling"
            fq = couple_fastq(${sample})
            gunzip = f"gunzip -c {fq} > {fq_dir}/${sample}.tmp.decompressed.fastq"
            os.system(gunzip)
            py = "/orfeo/cephfs/scratch/burlo/nardone/vcf_SV_valentina/ont_benchmark/test_tools/ont_benchmark/src/CAMPHOR.py"
            camphor = f"{py} ${bam.first()} {fq_dir}/${sample}.tmp.decompressed.fastq ${params.outdir}/camphor ${sample}.camphor"
            os.system(camphor)

        main()

        """
}

process CUTESV {
    tag {sample}_camphor

    publishDir("${params.outdir}/cutesv", mode: 'copy')

    input:
        tuple val(sample), path(bam)
    output:
        tuple val(sample), path(${sample}.cutesv.vcf.gz)
    script:
        """
        cuteSV ${bam.first()} /ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ${params.outdir}/cutesv/unsorted.vcf ${params.outdir}/cutesv -t 8 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 10 --diff_ratio_merging_DEL 0.3 &&
        bcftools sort -W -Oz -o ${params.outdir}/cutesv/${sample}.cutesv.vcf.gz ${params.outdir}/cutesv/unsorted.vcf && rm ${params.outdir}/cutesv/unsorted.vcf
        """
}

workflow {
    bams = Channel.fromFilePairs(${params.input_bams})

    DUET(bams)
    NANOSV(bams)
    SNIFFLES(bams)
    NANOVAR(bams)
    CAMPHOR(bams)
    CUTESV(bams)
}
