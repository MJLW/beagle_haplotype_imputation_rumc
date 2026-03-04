#!/usr/bin/env nextflow

plink_grch37 = "in/plink.GRCh37.map.zip"
ref_dir_grch37 = "in/b37.bref3/"
#ref_prefix = "chr"
#ref_postfix = ".1kg.phase3.v5a.b37.bref3"

beagle_jar = "bin/beagle.27Feb25.75f.jar"
unbref_jar = "bin/unbref3.27Feb25.75f.jar"

input_vcf_dir = "in/vcfs/"
input_vcf_pattern = "test*.vcf.gz"
output_dir = "out/"

process ConformVCF {
    tag "$trio_id"
    publishDir "${output_dir}/conformed_vcfs", mode: "link"

    input:
    tuple val(trio_id), val(relationship), path(vcf), path(vcf_csi)

    output:
    tuple val(trio_id), val(relationship), path("${trio_id}.${relationship}.vcf.gz"), path("${trio_id}.${relationship}.vcf.gz.csi")

    script:
    """
        # Mapping for removing 'chr' prefixes from regions
        awk 'BEGIN { for(i=1;i<=22;i++) print "chr"i"\t"i } END' > chr_map.txt

        # Filter to autosomal chromosomes only, then remove 'chr' prefix from regions
        bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 $vcf | bcftools annotate --rename-chrs chr_map.txt -Oz -o "${trio_id}.${relationship}.vcf.gz"
    """
}

process ChromosomeSplitVCF {
    tag "$trio_id"
    publishDir "${output_dir}/chr_vcfs", mode: "link"

    input:
    tuple val(trio_id), val(chr_nr), path(vcf), path(vcf_csi)

    output:
    tuple val(trio_id), val(relationship), val(chr_nr), path("${trio_id}.${relationship}.${chr_nr}.vcf.gz"), path("${trio_id}.${relationship}.${chr_nr}.vcf.gz.csi")


    script:
    """
        bcftools view -r $chr_nr -Oz -o ${trio_id}.${relationship}.${chr_nr}.vcf.gz $vcf
        bcftools index ${trio_id}.${relationship}.${chr_nr}.vcf.gz
    """
}

process Beagle {
    tag "$trio_id"
    publishDir "${output_dir}/imputed_chr_brefs", mode: "link"

    input:
    tuple val(trio_id), val(chr_nr), path(vcf), path(vcf_csi)

    output:
    tuple val(trio_id), val(relationship), val(chr_nr), path("{trio_id}.${relationship}.chr${chr_nr}.bref3")

    script:
    """
        java -jar $beagle_jar ref=$ref_dir_grch37 map=$plink_grch37 gt=$vcf out=${trio_id}.${relationship}.chr${chr_nr}.bref3
    """
}

process BRefToVCF {
    tag "$trio_id"
    publishDir "${output_dir}/imputed_chr_vcfs", mode: "link"

    input:
    tuple val(trio_id), val(relationship), val(chr_nr), path(bref3)

    output:
    tuple val(trio_id), val(relationship), val(chr_nr), path("${trio_id}.${relationship}.chr${chr_nr}.vcf.gz"), path("${trio_id}.${relationship}.chr${chr_nr}.vcf.gz.csi")

    script:
    """
        java -jar $unbref_jar $bref3 | bgzip > ${trio_id}.${relationship}.chr${chr_nr}.vcf.gz
        bcftools index ${trio_id}.${relationship}.chr${chr_nr}.vcf.gz
    """
}

process CombineVCF {
    tag "$trio_id"
    publishDir "${output_dir}/combined_vcfs", mode: "link"

    input:
    tuple val(trio_id), val(relationship), val(chr_nrs), path(vcfs), path(vcfs_csi)

    output:
    tuple val(trio_id), path("${trio_id}.${relationship}.vcf.gz"), path("${trio_id}.${relationship}.vcf.gz.csi")

    script:
    """
        # Concat and sort
        bcftools concat $vcfs | bcftools sort -Oz -o ${trio_id}.${relationship}.vcf.gz
        bcftools index ${trio_id}.${relationship}.vcf
    """
}


workflow {
    # Prepare input channel
    trio_dirs = Channel.fromPath("${input_vcf_dir}/rumc_trio_*")
    relationships = Channel.fromList(['proband', 'mother', 'father'])
    input_vcfs = trio_dirs.cross(relationships).map { trio_dir, relationship -> 
        def trio_id = trio_dir.getName()
        def vcf = file("${trio_dir}/${trio_id}.rescue.${relationship}.vcf.gz")
        def vcf_csi = file("${trio_dir}/${trio_id}.rescue.${relationship}.vcf.gz.csi")
        tuple(trio_id, relationship, vcf, vcf_csi)
    }

    # Conform Radboud VCF to 1000Genomes format
    # TODO: Split unphased variants into files kept in separate channel
    input_vcfs | ConformVCF | set { conformed_vcfs }

    # Add chromosome information
    chromosomes = Channel.of(1..22)
    chr_in_vcfs = conformed_vcfs.cross(chromosomes).map { trio_id, relationship, vcf, vcf_csi, chr_nr -> tuple(trio_id, relationship, chr_nr, vcf, vcf_csi) }

    # Split to chromosome level, run beagle, format and combine output
    chr_in_vcfs | ChromosomeSplitVCF | Beagle | BRefToVCF | set { chr_out_vcfs }

    # Combine chromosome level output VCFs into whole exome output VCFs
    grouped_chr_vcfs = chr_out_vcfs.groupTuple(by=[0,1], size=22) | CombineVCF | set { out_vcfs }

    # TODO: Merge unphased variants back into output vcfs
}

