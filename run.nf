#!/usr/bin/env nextflow

plink_dir_grch37 = "/ifs/data/research/projects/juliet/projects/prs_educational_attainment/in/plink.GRCh37/"
plink_prefix = "plink.chr"
plink_postfix = ".GRCh37.map"

ref_dir_grch37 = "/ifs/data/research/projects/juliet/projects/prs_educational_attainment/in/1000genomes_reference"
ref_prefix = "chr"
ref_postfix = ".1kg.phase3.v5a.b37.bref3"

beagle_jar = "/ifs/data/research/projects/juliet/projects/prs_educational_attainment/bin/beagle.27Feb25.75f.jar"
unbref_jar = "/ifs/data/research/projects/juliet/projects/prs_educational_attainment/bin/unbref3.27Feb25.75f.jar"

//parent_input_dir = "/ifs/data/research/projects/juliet/projects/prs_educational_attainment/in/test_vcfs"
parent_input_dir = "/ifs/data/research/projects/juliet/projects/inherited_ptvs/rerun/out"
trio_dir_pattern = "rumc_trio_*"
output_dir = "/ifs/data/research/projects/juliet/projects/prs_educational_attainment/test_out/"

chr_map_file = "/ifs/data/research/projects/juliet/projects/prs_educational_attainment/chr_map.txt"
summary_stats = "/ifs/data/research/projects/juliet/projects/prs_educational_attainment/snp_weights/EA4_additive_excl_23andMe.txt.gz"

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
        bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 $vcf | bcftools annotate --rename-chrs $chr_map_file -Oz -o "${trio_id}.${relationship}.temp.vcf.gz"
        bcftools reheader -s new_name.txt "${trio_id}.${relationship}.temp.vcf.gz" > "${trio_id}.${relationship}.vcf.gz"
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

process ImputeGenotypesBeagle {
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

// process CombineChromosomeVCFs {
//     tag "$chunk_id"
//     publishDir "${output_dir}/combined_vcfs/${chunk_id}", mode: "link"

//     input:
//     tuple val(chunk_id), val(chr_nrs), path(vcfs), path(vcfs_csi)

//     output:
//     tuple val(chunk_id), path("${chunk_id}.vcf.gz"), path("${chunk_id}.vcf.gz.csi")

//     script:
//     """
//         export TMPDIR=/ifs/temp/
//         module load bioinf/bcftools

//         # Concat and sort
//         bcftools concat $vcfs | bcftools sort -Oz -o ${chunk_id}.vcf.gz
//         bcftools index ${chunk_id}.vcf.gz
//     """
// }

// process ExtractSampleVCFs {
//     tag "${trio_id}.${relationship}"
//     publishDir "${output_dir}/sample_vcfs/${trio_id}", mode: "link"

//     input:
//     tuple val(trio_id), val(relationship), path(chunked_vcf), path(chunked_vcf_csi)

//     output:
//     tuple val(trio_id), val(relationship), path("${trio_id}.${relationship}.vcf.gz"), path("${trio_id}.${relationship}.vcf.gz.csi")

//     script:
//     """
//         export TMPDIR=/ifs/temp/
//         module load bioinf/bcftools

//         # Extract sample from
//     """
 // }

process MakeLDReference {
    tag "chr$chr_nr"
    publishDir "${output_dir}/1000G_plink", mode: "link"

    input:
    tuple val(chr_nr), val(ref_name), path(vcf), path(vcf_csi)

    output:
    tuple val(chr_nr), val("1000G_chr${chr_nr}"), path("1000G_chr${chr_nr}.bed"), path("1000G_chr${chr_nr}.bim"), path("1000G_chr${chr_nr}.fam")

    script:
    """
        eval "\$(conda shell.bash hook)"
        conda init
        conda activate /ifs/home/juliet/mamba/envs/plink_env

        plink --vcf $vcf --make-bed --out $ref_name
    """
}

process ConvertToPlink {
    tag "chr$chr_nr"
    publishDir "${output_dir}/plink", mode: "link"

    input:
    tuple val(chr_nr), path(vcf), path(vcf_csi)

    output:
    tuple val(chr_nr), path("rumc_chr${chr_nr}.bed"), path("rumc_chr${chr_nr}.bim"), path("rumc_chr${chr_nr}.fam")

    script:
    """
        eval "\$(conda shell.bash hook)"
        conda init
        conda activate /ifs/home/juliet/mamba/envs/plink_env

        plink --vcf $vcf --make-bed --out rumc_chr${chr_nr}
    """
}

process AdjustSNPWeights {
    tag "chr$chr_nr"
    publishDir "${output_dir}/adjusted_weights", mode: "link"
    memory '125 GB'

    input:
    tuple val(chr_nr), val(ref_name), path(ref_bed), path(ref_bim), path(ref_fam)

    output:
    tuple val(chr_nr), path("chr${chr_nr}_ldradius405.pkl"), path("snp_weights.chr${chr_nr}_LDpred_*.txt")

    script:
    """
        zcat $summary_stats | head -n 1 > summary_stats_${chr_nr}.txt
        awk -F '\\t' '\$2 == ${chr_nr}' $summary_stats >> summary_stats_${chr_nr}.txt

        ldpred coord --gf $ref_name --ssf summary_stats_${chr_nr}.txt --ssf-format CUSTOM --rs rsID --A1 Effect_allele --A2 Other_allele --pos BP --chr Chr --pval P --eff Beta --se SE --N 765283 --out sumstats_coords.chr${chr_nr}.hdf5
        ldpred gibbs --cf sumstats_coords.chr${chr_nr}.hdf5 --ldr 405 --ldf chr${chr_nr} --out snp_weights.chr${chr_nr} --N 765283 --no-ld-compression
    """
}

process ProcessSNPs {
    // errorStrategy "ignore"
    tag "ProcessSNPs"

    input:
    path(snp_weights)

    output:
    path("EA_exome_weights.fixed.txt")

    publishDir "snp_weights/pruned_snps", mode: "link"

    script:
    """
    awk 'BEGIN {OFS="\\t"; print "CHROM","POS","EFFECT","OTHER","AF","BETA"} NR>1 {print "chr"\$2,\$3,\$4,\$5,\$6,\$8}' $snp_weights | (head -n 1 && tail -n +2 | sort -k1.4n -k2,2n) > snps_in.txt

    /ifs/data/research/projects/juliet/tools/haplotype_prs_rs/target/release/fix_snps -f ${params.ref_fa} -s snps_in.txt -o EA_exome_weights.fixed.txt
    """
}


workflow {
    // Prepare input channel
    trio_dirs = Channel.fromPath("/ifs/data/research/projects/juliet/projects/inherited_ptvs/rerun/out/rumc_trio_*", type: 'dir')

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

    conformed_vcfs
        .branch {
            proband: it[1] == 'proband'
            non_proband: it[1] != 'proband'
        }
        .set { split_vcfs }

    proband_conformed_vcfs = split_vcfs.proband
    parent_conformed_vcfs = split_vcfs.non_proband

    def enumerator = 0
    // toSortedList{}.flatMap{} is required for resume to work, it's stupid but it works
    proband_chunks = proband_conformed_vcfs.toSortedList { a, b -> a[0] <=> b[0] }.flatMap { it }.collate(1000).map { it -> tuple(enumerator++, *(it.transpose())) }
    parent_chunks = parent_conformed_vcfs.toSortedList { a, b -> a[0] <=> b[0] }.flatMap { it }.collate(1000).map { it -> tuple(enumerator++, *(it.transpose())) }

    // chunk_mapping = chunked_vcfs.map { chunk_id, trio_ids, relationships, vcfs, vcf_csis -> tuple(chunk_id, trio_ids, relationships) }
    proband_chunked_vcfs = proband_chunks.map { chunk_id, trio_ids, relationships, vcfs, vcf_csis -> tuple(chunk_id, vcfs, vcf_csis) }
    parent_chunked_vcfs = parent_chunks.map { chunk_id, trio_ids, relationships, vcfs, vcf_csis -> tuple(chunk_id, vcfs, vcf_csis) }

    proband_chunked_vcfs | MergeVCFs | set { proband_merged_vcfs }
    parent_chunked_vcfs | MergeVCFs | set { parent_merged_vcfs }

    // Add chromosome information
    chromosomes = Channel.of(1..22)
    proband_chr_in_vcfs = proband_merged_vcfs.combine(chromosomes).map { chunk_id, vcf, vcf_csi, chr_nr -> tuple(chunk_id, chr_nr, vcf, vcf_csi) }
    parent_chr_in_vcfs = parent_merged_vcfs.combine(chromosomes).map { chunk_id, vcf, vcf_csi, chr_nr -> tuple(chunk_id, chr_nr, vcf, vcf_csi) }

    // Split to chromosome level, run beagle, format and combine output
    proband_chr_in_vcfs | ChromosomeSplitVCF | ImputeGenotypesBeagle | set { proband_chr_out_vcfs }
    parent_chr_in_vcfs | ChromosomeSplitVCF | ImputeGenotypesBeagle | set { parent_chr_out_vcfs }

    // Combine chunks
    proband_grouped_by_chr_vcfs = proband_chr_out_vcfs.groupTuple(by: 1).map { chunk_ids, chr_nr, vcfs, vcf_csis -> tuple(chr_nr, vcfs, vcf_csis) }
    proband_grouped_by_chr_vcfs | MergeChromosomeChunks | set { proband_merged_chr_vcfs }

    parent_grouped_by_chr_vcfs = parent_chr_out_vcfs.groupTuple(by: 1).map { chunk_ids, chr_nr, vcfs, vcf_csis -> tuple(chr_nr, vcfs, vcf_csis) }
    parent_grouped_by_chr_vcfs | MergeChromosomeChunks | set { parent_merged_chr_vcfs }

    // Adjust GWAS summary SNP weights using parents as LD reference
    reference_input = parent_merged_chr_vcfs.map { chr_nr, vcf, vcf_csi -> tuple(chr_nr, "RUMC_parents_chr${chr_nr}", vcf, vcf_csi) }
    reference_input | MakeLDReference | AdjustSNPWeights | set { ld_reference }

    // Adjust SNP weights using 1000G dataset
    // reference_input = chromosomes.map { it ->
    //     def chr_nr = it
    //     def ref_name = "1000G_chr${chr_nr}"
    //     def vcf = "/ifs/data/research/projects/juliet/projects/prs_educational_attainment/in/1000genomes_genotypes/ALL.chr${chr_nr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    //     tuple(chr_nr, ref_name, vcf)
    // }
    // reference_input | MakeLDReference | AdjustSNPWeights | set { ld_reference }

    // Convert RUMC to Plink
    // merged_chr_vcfs | ConvertToPlink | set { ld_target }


    // TODO: Subset RUMC VCF to parents to create dataset of unrelated samples


    // TODO: Merge unphased variants back into output vcfs
}
