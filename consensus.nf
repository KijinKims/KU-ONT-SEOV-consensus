nextflow.enable.dsl=2

workflow {
    main:
        if( params.help ) {
            help = """consensus.nf: Nextflow script to generate the consensus sequence of SEOV from Nanopore sequenicng datset
                    |Required arguments:
                    |  --fastq  Location of the input FASTQ file.
                    |  --prefix  Prefix of the run.
                    |  --outdir  Location output consensus will be saved to.
                    |
                    |Optional arguments:
                    |  --L   L segment reference in FASTA format.
                    |  --M   M segment reference in FASTA format.
                    |  --S   S segment reference in FASTA format.""".stripMargin()
            // Print the help with the stripped margin and exit
            println(help)
            exit(0)
        }

        Channel.fromPath(params.fastq).set{fastq}
        Channel.empty().set{init}
        Channel.fromPath(params.L).set{lseg}
        Channel.fromPath(params.M).set{mseg}
        Channel.fromPath(params.S).set{sseg}
        init.concat(lseg,mseg,sseg).set{segments}
        consensus(segments, fastq)
}

workflow consensus {
    take:
        segments
        fastq
    
    main:
        segments.combine(fastq).set{ref_fastq}
        consensus_process(ref_fastq)
}

process consensus_process {
    conda 'seov_consensus'
    publishDir "${params.outdir}", mode: 'copy'
    tag "${params.prefix}:consensus"

    input:
        tuple path(ref), path(fastq)
    output:
        path "${params.prefix}_${ref.simpleName}.consensus.fasta"
        path "bams/${params.prefix}_${ref.simpleName}.bam"
    """
    # variant call
    mini_align -t 8 -p ${params.prefix}_${ref.simpleName} -i $fastq -r $ref -m -f -M 2 -S 4 -O 4,24 -E 2,1
    medaka consensus --model r941_prom_variant_g360 --threads 8 --batch_size 100 ${params.prefix}_${ref.simpleName}.bam ${params.prefix}_${ref.simpleName}.hdf
    medaka variant $ref ${params.prefix}_${ref.simpleName}.hdf ${params.prefix}_${ref.simpleName}.medaka.vcf
    medaka tools annotate --dpsp ${params.prefix}_${ref.simpleName}.medaka.vcf $ref ${params.prefix}_${ref.simpleName}.bam ${params.prefix}_${ref.simpleName}.medaka.annotated.vcf
    bcftools view -i 'QUAL>${params.variant_quality_threshold} & INFO/DP>${params.variant_depth_threshold}' ${params.prefix}_${ref.simpleName}.medaka.annotated.vcf > ${params.prefix}_${ref.simpleName}.medaka.filtered.vcf
    bgzip -f ${params.prefix}_${ref.simpleName}.medaka.filtered.vcf
    tabix -p vcf ${params.prefix}_${ref.simpleName}.medaka.filtered.vcf.gz

    # apply variants to reference to generate consensus
    cat $ref | bcftools consensus ${params.prefix}_${ref.simpleName}.medaka.filtered.vcf.gz > ${ref.simpleName}.consensus.fasta

    # drop low coverage region
    bedtools genomecov -bga -ibam ${params.prefix}_${ref.simpleName}.bam | awk '\$4 < ${params.low_cov_threshold}' | bedtools merge -i - > low_cov_${params.low_cov_threshold}.bed
    bedtools maskfasta -fi ${ref.simpleName}.consensus.fasta -bed low_cov_${params.low_cov_threshold}.bed -fo ${params.prefix}_${ref.simpleName}.consensus.fasta
    mkdir bams
    mv ${params.prefix}_${ref.simpleName}.bam bams/

    # change fasta header
    header=">${params.prefix}_${ref.simpleName}_consensus low_cov_thrs=${params.low_cov_threshold} var_qual_thrs=${params.variant_quality_threshold} var_depth_thrs=${params.variant_depth_threshold}"
    sed -i "1s/.*/\$header/" ${params.prefix}_${ref.simpleName}.consensus.fasta
    """
}