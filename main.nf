#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input_bams = "all"
params.outdir = "/orfeo/cephfs/scratch/burlo/nardone/vcf_SV_valentina/ont_benchmark/SVCALL"

//bams = Channel.fromFilePairs(${params.input_bams})

process NANOSV {
    tag {"${sample}_nanosv"}

    publishDir("${params.outdir}/nanosv", mode: 'copy')

    input:
        tuple val(sample), path(bam)
    output:
        tuple val(sample), path("${sample}.nanosv.vcf.gz")
    script:
        """
        NanoSV ${bam.first()} -o NanoSV.tmp.vcf -b /orfeo/LTS/burlo/LT_storage/nardone/software/NanoSV/human_hg38.bed -t 8 && \
        bcftools view -i 'INFO/SVTYPE == "DEL"' NanoSV.tmp.vcf | bcftools sort -Oz -o ${params.outdir}/nanosv/${sample}.nanosv.vcf.gz && bcftools index ${params.outdir}/nanosv/${sample}.nanosv.vcf.gz
        """
}

process DUET {
    tag {"${sample}_duet"}

    publishDir("${params.outdir}/duet", mode: 'copy')

    input:
        tuple val(sample), path(bam)
    output:
        tuple val(sample), path("${sample}.duet.vcf.gz")
    script:
        """
        duet \$(realpath ${bam.first()}) /ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \${PWD} -t 8 &&
        sed '5 i\##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' phased_sv.vcf | bcftools view -i 'INFO/SVTYPE == "DEL"' | bcftools sort -Oz -o ${sample}.duet.vcf.gz &&
        bcftools index ${sample}.duet.vcf.gz
        """
}

process SNIFFLES {
    tag {"${sample}_sniffles"}

    publishDir("${params.outdir}/sniffles", mode: 'copy')

    input:
        tuple val(sample), path(bam)
    output:
        tuple val(sample), path("${sample}.sniffles.vcf.gz")
    script:
        """
        sniffles -t 8 --tandem-repeats /orfeo/cephfs/scratch/burlo/nardone/vcf_SV_valentina/ont_benchmark/resources/human_GRCh38_no_alt_analysis_set.trf.bed --input ${bam.first()} --vcf ${sample}.sniffles.vcf.gz --reference /ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna && bcftools index -f -t ${sample}.sniffles.vcf.gz
        """
}

process NANOVAR {
    tag {"${sample}_nanovar"}

    publishDir("${params.outdir}/nanovar", mode: 'copy')

    input:
        tuple val(sample), path(bam)
    output:
        tuple val(sample), path("${sample}.nanovar.vcf.gz")
    script:
        """
        nanovar ${bam.first()} /ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna . -t 8 -x ont -l 50 && 
        bcftools view --write-index -Oz -o ${sample}.nanovar.vcf.gz ./${sample}.nanovar.pass.vcf
        """
}


process CAMPHOR {
    tag {"${sample}_camphor"}

    publishDir("${params.outdir}/camphor", mode: 'copy')

    input:
        tuple val(sample), path(bam)
    output:
        tuple val(sample), path("${sample}.camphor_DEL_only_sorted.vcf.gz")
    script:        
        """
        #!/usr/bin/env python3

        from glob import glob
        import re
        import os

        def couple_fastq(sample):
            cov_str = sample.split("_")[-1]
            cov_list = re.findall(r'\\d+', cov_str)
            cov = "".join(cov_list)
            pattern = f"/orfeo/cephfs/scratch/burlo/nardone/vcf_SV_valentina/ont_benchmark/FASTQS/downsampling/HG002_minimap_down_to_{cov}X.fastq.gz"
            search = glob(pattern)
            fq = "".join(search)
            return fq

        def main():
            fq_dir = "/orfeo/cephfs/scratch/burlo/nardone/vcf_SV_valentina/ont_benchmark/FASTQS/downsampling"
            fq = couple_fastq("${sample}")
            gunzip = f"gunzip -c {fq} > {fq_dir}/${sample}.tmp.decompressed.fastq"
            os.system(gunzip)
            bam = os.path.realpath("${bam.first()}")
            py = "/orfeo/cephfs/scratch/burlo/nardone/vcf_SV_valentina/ont_benchmark/test_tools/ont_benchmark/src/CAMPHOR.py"
            camphor = f"{py} {bam} {fq_dir}/${sample}.tmp.decompressed.fastq {os.getcwd()} ${sample}.camphor"
            os.system(camphor)

        main()

        """
}

process CUTESV {
    tag {"${sample}_cuteSV"}

    publishDir("${params.outdir}/cutesv", mode: 'copy')

    input:
        tuple val(sample), path(bam)
    output:
        tuple val(sample), path("${sample}.cutesv.vcf.gz")
    script:
        """
        cuteSV ${bam.first()} /ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \${PWD}/unsorted.vcf \${PWD} -t 8 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 10 --diff_ratio_merging_DEL 0.3 &&
        bcftools sort -Oz -o ${sample}.cutesv.vcf.gz \${PWD}/unsorted.vcf && bcftools index ${sample}.cutesv.vcf.gz && rm \${PWD}/unsorted.vcf
        """
}

process BENCHMARK {
    tag {"${sample}_benchmark"}

    publishDir("${params.outdir}/benchmark/${sample}", mode: 'copy')

    input:
        tuple val(sample), path(vcfs)
    output:
        tuple val(sample), path("${sample}_benchmark.csv")
    script:
        """
        wittyer -i ${vcfs} -t /orfeo/cephfs/scratch/burlo/nardone/vcf_SV_valentina/ont_benchmark/GIAB_reference/lifted/HG002_SVs_Tier1_v0.6.hg38.vcf.gz \
        -o ${params.outdir}/benchmark/${sample} -em CrossTypeAndSimpleCounting --configFile /orfeo/cephfs/scratch/burlo/nardone/vcf_SV_valentina/ont_benchmark/GIAB_reference/config_witty_SV.json &&
        printf "Method,SV_Type,Interval,TP,FP,FN,Precision,Recall,F1Score\\n" > ${sample}_benchmark.csv
        jq -r --arg folder ${sample} '.PerSampleStats[0].DetailedStats[] | select(.VariantType=="Deletion+CopyNumberLoss") | .PerBinStats[] | .Bin as \$bin | .Stats[] | select(.StatsType=="Event") | [ \$folder, "DEL", \$bin, .TruthTpCount, .QueryFpCount, .TruthFnCount, .Precision, .Recall, .Fscore] | @csv' ${params.outdir}/benchmark/${sample}/Wittyer.Stats.json >> ${sample}_benchmark.csv
        jq -r --arg folder ${sample} '.PerSampleStats[0].DetailedStats[] | select(.VariantType=="Deletion+CopyNumberLoss") | .OverallStats[0] | [ \$folder, "DEL", "ALL", .TruthTpCount, .QueryFpCount, .TruthFnCount, .Precision, .Recall, .Fscore] | @csv' ${params.outdir}/benchmark/${sample}/Wittyer.Stats.json >> ${sample}_benchmark.csv
        """
}

def baseDir = "/orfeo/cephfs/scratch/burlo/nardone/vcf_SV_valentina/ont_benchmark/BAM/downsampling"

//def generateFilePatterns(keywords) {
//    return keywords.split(',').collect { keyword ->
//        "${baseDir}/${keyword}/HG002*.bam{,.bai}"
//    }
//}

//if (params.input_bams == 'all') {
//    bams = Channel.fromFilePairs("${baseDir}/*/HG002*.bam{,.bai}")
//} else {
//    def filePatterns = generateFilePatterns(params.input_bams)
//    bams = Channel.empty()
//    filePatterns.each { pattern ->
//        bams = bams.mix(Channel.fromFilePairs(pattern))
//    }
//}


workflow {
    if (params.input_bams == 'all') {
    bams = Channel.fromFilePairs("${baseDir}/*/HG002*.bam{,.bai}")
    } else {
        //def filePatterns = generateFilePatterns(params.input_bams)
        //bams = Channel.empty()
        //filePatterns.each { pattern ->
        //    bams = bams.mix(Channel.fromFilePairs(pattern))
        //}
        bams = Channel.fromFilePairs("${baseDir}/{${params.input_bams}}/HG002*.bam{,.bai}")
    }

    //DUET(bams)
    NANOSV(bams)
    SNIFFLES(bams)
    NANOVAR(bams)
    CAMPHOR(bams)
    CUTESV(bams)

    results_to_bench = NANOSV.out
                      .join(SNIFFLES.out)
                      .join(NANOVAR.out)
                      .join(CAMPHOR.out)
                      .join(CUTESV.out)
                      
    BENCHMARK(results_to_bench)
}

