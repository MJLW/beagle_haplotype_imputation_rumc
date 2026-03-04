#!/usr/bin/env nextflow

plink_dir_grch37 = "/ifs/data/research/projects/juliet/projects/prs_educational_attainment/in/plink.GRCh37/"
plink_prefix = "plink.chr"
plink_postfix = ".GRCh37.map"

ref_dir_grch37 = "/ifs/data/research/projects/juliet/projects/prs_educational_attainment/in/1000genomes_reference"
ref_prefix = "chr"
ref_postfix = ".1kg.phase3.v5a.b37.bref3"

beagle_jar = "/ifs/data/research/projects/juliet/projects/prs_educational_attainment/bin/beagle.27Feb25.75f.jar"
unbref_jar = "/ifs/data/research/projects/juliet/projects/prs_educational_attainment/bin/unbref3.27Feb25.75f.jar"

parent_input_dir = "/ifs/data/research/projects/juliet/projects/prs_educational_attainment/in/test_vcfs"
trio_dir_pattern = "rumc_trio_*"
output_dir = "/ifs/data/research/projects/juliet/projects/prs_educational_attainment/test_out/"

chr_map_file = "/ifs/data/research/projects/juliet/projects/prs_educational_attainment/chr_map.txt"

process ConformVCF {
    tag "${trio_id}.${relationship}"
    publishDir "${output_dir}/conformed_vcfs/${trio_id}", mode: "link"

    input:
    tuple val(trio_id), val(relationship), path(vcf), path(vcf_csi)

    output:
    tuple val(trio_id), val(relationship), path("${trio_id}.${relationship}.vcf.gz"), path("${trio_id}.${relationship}.vcf.gz.csi")

    script:
    """
        export TMPDIR=/ifs/temp/
        module load bioinf/bcftools

        echo "${trio_id}.${relationship}" > new_name.txt

        # Filter to autosomal chromosomes only, then remove 'chr' prefix from regions
        bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 $vcf | bcftools annotate --rename-chrs $chr_map_file | bcftools reheader -s new_name.txt -Oz -o "${trio_id}.${relationship}.vcf.gz"
        bcftools index "${trio_id}.${relationship}.vcf.gz"
    """
}

process MergeVCFs {
    tag "$chunk_id"
    publishDir "${output_dir}/chunked_vcfs/${chunk_id}", mode: "link"

    input:
    tuple val(chunk_id), path(vcfs), path(vcf_csis)

    output:
    tuple val(chunk_id), path("${chunk_id}.vcf.gz"), path("${chunk_id}.vcf.gz.csi")

    script:
    """
        export TMPDIR=/ifs/temp/
        module load bioinf/bcftools

        bcftools merge -Oz -o ${chunk_id}.vcf.gz $vcfs
        bcftools index ${chunk_id}.vcf.gz
    """
}

process ChromosomeSplitVCF {
    tag "$chunk_id"
    publishDir "${output_dir}/chr_vcfs/${chunk_id}", mode: "link"

    input:
    tuple val(chunk_id), val(chr_nr), path(vcf), path(vcf_csi)

    output:
    tuple val(chunk_id), val(chr_nr), path("${chunk_id}.${chr_nr}.vcf.gz"), path("${chunk_id}.${chr_nr}.vcf.gz.csi")


    script:
    """
        export TMPDIR=/ifs/temp/
        module load bioinf/bcftools

        bcftools view -r $chr_nr -Oz -o ${chunk_id}.${chr_nr}.vcf.gz $vcf
        bcftools index ${chunk_id}.${chr_nr}.vcf.gz
    """
}

process Beagle {
    tag "$chunk_id"
    publishDir "${output_dir}/imputed_chr_brefs/${chunk_id}", mode: "link"

    input:
    tuple val(chunk_id), val(chr_nr), path(vcf), path(vcf_csi)

    output:
    tuple val(chunk_id), val(chr_nr), path("${chunk_id}.chr${chr_nr}.vcf.gz"), path("${chunk_id}.chr${chr_nr}.vcf.gz.csi")

    script:
    """
        export TMPDIR=/ifs/temp/
        module load bioinf/bcftools

        java -jar $beagle_jar ref=${ref_dir_grch37}/${ref_prefix}${chr_nr}${ref_postfix} map=${plink_dir_grch37}/${plink_prefix}${chr_nr}${plink_postfix} gt=$vcf out=${chunk_id}.chr${chr_nr}
        bcftools index ${chunk_id}.chr${chr_nr}.vcf.gz
    """
}

process MergeChromosomeChunks {
    tag "chr${chr_nr}"
    publishDir "${output_dir}/merged_chrs/${chr_nr}", mode: "link"

    input:
    tuple val(chr_nr), path(vcfs), path(vcf_csis)

    output:
    tuple val(chr_nr), path("chr${chr_nr}.vcf.gz"), path("chr${chr_nr}.vcf.gz.csi")

    script:
    """
        export TMPDIR=/ifs/temp/
        module load bioinf/bcftools

        bcftools merge -Oz -o chr${chr_nr}.vcf.gz $vcfs
        bcftools index chr${chr_nr}.vcf.gz
    """
}

# process CombineChromosomeVCFs {
#     tag "$chunk_id"
#     publishDir "${output_dir}/combined_vcfs/${chunk_id}", mode: "link"
#
#     input:
#     tuple val(chunk_id), val(chr_nrs), path(vcfs), path(vcfs_csi)
#
#     output:
#     tuple val(chunk_id), path("${chunk_id}.vcf.gz"), path("${chunk_id}.vcf.gz.csi")
#
#     script:
#     """
#         export TMPDIR=/ifs/temp/
#         module load bioinf/bcftools
#
#         # Concat and sort
#         bcftools concat $vcfs | bcftools sort -Oz -o ${chunk_id}.vcf.gz
#         bcftools index ${chunk_id}.vcf.gz
#     """
# }

# process ExtractSampleVCFs {
#     tag "${trio_id}.${relationship}"
#     publishDir "${output_dir}/sample_vcfs/${trio_id}", mode: "link"
#
#     input:
#     tuple val(trio_id), val(relationship), path(chunked_vcf), path(chunked_vcf_csi)
#
#     output:
#     tuple val(trio_id), val(relationship), path("${trio_id}.${relationship}.vcf.gz"), path("${trio_id}.${relationship}.vcf.gz.csi")
#
#     script:
#     """
#         export TMPDIR=/ifs/temp/
#         module load bioinf/bcftools
#
#         # Extract sample from 
#     """
# }


workflow {
    // Prepare input channel
    trio_dirs = Channel.fromPath("/ifs/data/research/projects/juliet/projects/prs_educational_attainment/in/test_vcfs/rumc_trio_*", type: 'dir')

    // trio_dirs = Channel.fromPath("${parent_input_dir}/${trio_dir_pattern}", type: 'dir')
    relationships = Channel.of('proband', 'mother', 'father')
    input_vcfs = trio_dirs.combine(relationships).map { trio_dir, relationship -> 
        def trio_id = trio_dir.getName()
        def vcf = file("${trio_dir}/${trio_id}.rescue.${relationship}.vcf.gz")
        def vcf_csi = file("${trio_dir}/${trio_id}.rescue.${relationship}.vcf.gz.csi")
        tuple(trio_id, relationship, vcf, vcf_csi)
    }

    // Conform Radboud VCF to 1000Genomes format
    // TODO: Split unphased variants into files kept in separate channel
    input_vcfs | ConformVCF | set { conformed_vcfs }

    def enumerator = 0
    chunks = conformed_vcfs.collate(500).map { it -> tuple(enumerator++, *data) }
    # chunk_mapping = chunked_vcfs.map { chunk_id, trio_ids, relationships, vcfs, vcf_csis -> tuple(chunk_id, trio_ids, relationships) }
    chunked_vcfs = chunks.map { chunk_id, trio_ids, relationships, vcfs, vcf_csis -> chunk_id, vcfs, vcf_csis }

    chunked_vcfs | MergeVCFs | set { merged_vcfs }

    // Add chromosome information
    chromosomes = Channel.of(1..22)
    chr_in_vcfs = conformed_vcfs.combine(chromosomes).map { trio_id, relationship, vcf, vcf_csi, chr_nr -> tuple(trio_id, relationship, chr_nr, vcf, vcf_csi) }

    // Split to chromosome level, run beagle, format and combine output
    chr_in_vcfs | ChromosomeSplitVCF | Beagle | set { chr_out_vcfs }

    // Combine chromosome level output VCFs into whole exome output VCFs
    grouped_by_chr_vcfs = chr_out_vcfs.groupTuple(by: 1)
    grouped_by_chr_vcfs | MergeChromosomeChunks | set { merged_chr_vcfs }


    // TODO: Merge unphased variants back into output vcfs
}

