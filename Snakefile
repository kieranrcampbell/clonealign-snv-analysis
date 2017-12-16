"""
Variant calling for clones and 10x RNA-seq reads

"""


import pandas as pd

ref_genome = "/shahlab/kicampbell/reference/GRCh37-lite.fa"

clones = ["A", "B", "C"]

rna_sample_names = "SA501X2B"
dna_sample_names = "SA501X3F"

cell_barcode_df = pd.read_csv("data/{}/rna/{}_cell_barcodes_for_em_Dec_15th.csv".format(rna_sample_names, rna_sample_names))
cell_barcodes = cell_barcode_df['barcode'].tolist()



dna_clone_variant_bcfs = expand("data/{dna_sample_name}/dna/variants/bcf/{dna_sample_name}-clone_{clone}_variants.bcf",
dna_sample_name = dna_sample_names, clone = clones)
dna_clone_variant_vcfs = expand("data/{dna_sample_name}/dna/variants/vcf/{dna_sample_name}-clone_{clone}_variants.vcf",
dna_sample_name = dna_sample_names, clone = clones)
dna_initial_snvs = expand("data/{dna_sample_name}/dna/variants/bed_initial_snvs/{dna_sample_name}-clone_{clone}_variants.bed",
dna_sample_name = dna_sample_names, clone = clones)
dna_clone_specific_beds_unfiltered = expand("data/{dna_sample_name}/dna/variants/clone_specific_beds/{dna_sample_name}-clone_{clone}_SNVs.bed", dna_sample_name = dna_sample_names, clone = clones)

dna_depths = expand("data/{dna_sample_name}/dna/variants/clone_variant_depths/{dna_sample_name}_SNVs_from_clone_{snv_clone}_depth_in_clone_{bam_clone}.tsv", dna_sample_name = dna_sample_names, snv_clone = clones, bam_clone = clones)

dna_clone_specific_beds_filtered = expand("data/{dna_sample_name}/dna/variants/clone_specific_beds/{dna_sample_name}-clone_{clone}_SNVs_filtered.bed", dna_sample_name = dna_sample_names, clone = clones)



rna_cell_bams = expand("data/rna/bam/{rna_sample_name}_cell_{barcode}.bam", 
                        rna_sample_name = rna_sample_names, barcode = cell_barcodes)
rna_cell_pileup = expand("data/{rna_sample_name}/rna/variants/bcf_by_clone/{rna_sample_name}-cell_{barcode}_variants_for_clone_{clone}_in_{dna_sample_name}.bcf", rna_sample_name = rna_sample_names, barcode = cell_barcodes, clone = clones, dna_sample_name = dna_sample_names)
rna_cell_variants = expand("data/{rna_sample_name}/rna/variants/vcf_by_clone/{rna_sample_name}-cell_{barcode}_variants_for_clone_{clone}_in_{dna_sample_name}.vcf", rna_sample_name = rna_sample_names, barcode = cell_barcodes, clone = clones, dna_sample_name = dna_sample_names)
collate_logs = expand("data/{rna_sample_name}/rna/variants/logs/{rna_sample_name}-cell_{barcode}_variants_for_clone_{clone}_in_{dna_sample_name}.txt",
rna_sample_name = rna_sample_names, barcode = cell_barcodes, clone = clones, dna_sample_name = dna_sample_names)


# rna_cell_bams = expand("data/rna/bam/{sample_name}_cell_{barcode}.bam", 
#                         sample_name = rna_sample_name, barcode = cell_barcodes)
# rna_cell_fastq = expand("data/rna/fastq/{sample_name}_cell_{barcode}.fastq", 
#                         sample_name = rna_sample_name, barcode = cell_barcodes)
# rna_cell_bams_realigned = expand("data/rna/bam_realigned/{sample_name}_cell_{barcode}.bam", 
#                         sample_name = rna_sample_name, barcode = cell_barcodes)

# rna_clone_variant_bcfs = expand("data/rna/variants/bcf/{sample_name}-cell_{barcode}_variants.bcf",
#                                 sample_name = rna_sample_name, barcode = cell_barcodes)
# rna_clone_variant_vcfs = expand("data/rna/variants/vcf/{sample_name}-cell_{barcode}_variants.vcf",
                                # sample_name = rna_sample_name, barcode = cell_barcodes)


rule all:
    input:
        dna_initial_snvs,
        dna_clone_specific_beds_unfiltered,
        dna_depths,
        dna_clone_specific_beds_filtered,
        # rna_cell_pileup#,
        collate_logs



# DNA specific analysis here ------

rule dna_pileup:
    input:
        bamfile = "data/{dna_sample_name}/dna/bam/{dna_sample_name}-cluster_{clone}.sorted.realigned.rmdups.bam",
        ref = ref_genome
    output:
        "data/{dna_sample_name}/dna/variants/bcf/{dna_sample_name}-clone_{clone}_variants.bcf"
    shell:
        "samtools mpileup -R -g -f {input.ref} {input.bamfile} > {output}"

rule dna_bcf_to_vcf:
    input:
        "data/{dna_sample_name}/dna/variants/bcf/{dna_sample_name}-clone_{clone}_variants.bcf"
    output:
        "data/{dna_sample_name}/dna/variants/vcf/{dna_sample_name}-clone_{clone}_variants.vcf"
    shell:
        "bcftools call -c -v {input} > {output}"

rule dna_vcf_to_bed:
    input:
        "data/{dna_sample_name}/dna/variants/vcf/{dna_sample_name}-clone_{clone}_variants.vcf"
    output:
        "data/{dna_sample_name}/dna/variants/bed_initial_snvs/{dna_sample_name}-clone_{clone}_variants.bed"
    shell:
        "vcf2bed < {input} > {output}"

rule temp_intersect:
    input:
        bedA = "data/{dna_sample_name}/dna/variants/bed_initial_snvs/{dna_sample_name}-clone_A_variants.bed",
        bedB = "data/{dna_sample_name}/dna/variants/bed_initial_snvs/{dna_sample_name}-clone_B_variants.bed"
    output:
        "data/{dna_sample_name}/dna/variants/tmp_intersect/SA501X3F-intersect_AB.bed"
    shell:
        "bedtools intersect -a {input.bedA} -b {input.bedB} > {output}"

rule germline_intersect:
    input:
        bedAB = "data/{dna_sample_name}/dna/variants/tmp_intersect/SA501X3F-intersect_AB.bed",
        bedC = "data/{dna_sample_name}/dna/variants/bed_initial_snvs/{dna_sample_name}-clone_A_variants.bed"
    output:
        "data/{dna_sample_name}/dna/variants/germline_intersect/{dna_sample_name}-germline_SNVs.bed"
    shell:
        "bedtools intersect -a {input.bedAB} -b {input.bedC} > {output}"

rule clone_specific_beds:
    input:
       bed = "data/{dna_sample_name}/dna/variants/bed_initial_snvs/{dna_sample_name}-clone_{clone}_variants.bed",
       germline =  "data/{dna_sample_name}/dna/variants/germline_intersect/{dna_sample_name}-germline_SNVs.bed"
    output:
        "data/{dna_sample_name}/dna/variants/clone_specific_beds/{dna_sample_name}-clone_{clone}_SNVs.bed"
    shell:
        "bedtools subtract -a {input.bed} -b {input.germline} > {output}"

rule dna_depth:
    input:
        snvs = "data/{dna_sample_name}/dna/variants/clone_specific_beds/{dna_sample_name}-clone_{snv_clone}_SNVs.bed",
        bam = "data/{dna_sample_name}/dna/bam/{dna_sample_name}-cluster_{bam_clone}.sorted.realigned.rmdups.bam"
    output:
        "data/{dna_sample_name}/dna/variants/clone_variant_depths/{dna_sample_name}_SNVs_from_clone_{snv_clone}_depth_in_clone_{bam_clone}.tsv"
    shell:
        "samtools depth -b {input.snvs} {input.bam} > {output}"

rule filter_beds:
    input:
        bed="data/{dna_sample_name}/dna/variants/clone_specific_beds/{dna_sample_name}-clone_{clone}_SNVs.bed",
        depth=dna_depths
    output:
        "data/{dna_sample_name}/dna/variants/clone_specific_beds/{dna_sample_name}-clone_{clone}_SNVs_filtered.bed"
    shell:
        """
        rm -f {output}
        python3 scripts/filter_snv.py --input_bed {input.bed} --depth_template data/{wildcards.dna_sample_name}/dna/variants/clone_variant_depths/{wildcards.dna_sample_name}_SNVs_from_clone_{wildcards.clone}_depth_in_clone_CLONE.tsv --clone {wildcards.clone} --output_bed {output}
        """


# python scripts/filter_snv.py --input_bed data/SA501X3F/dna/variants/clone_specific_beds/SA501X3F-clone_C_SNVs.bed --depth_template data/SA501X3F/dna/variants/clone_variant_depths/SA501X3F_SNVs_from_clone_C_depth_in_clone_CLONE.tsv --clone C --output_bed tmp.bed

# RNA here -----------------



rule split_rna_sam_on_barcode:
    input:
        "data/{}/rna/bam_from_10x/{}_possorted_genome_bam.bam".format(rna_sample_names, rna_sample_names)
    output:
        "data/{rna_sample_name}/rna/bam/{rna_sample_name}_cell_{barcode}.bam"
    shell:
        "samtools view -h {input} | grep '^\\@\\|{wildcards.barcode}' | samtools view -Sb -o {output} -"


rule rna_pileup:
    input:
        clone_bedfile="data/{dna_sample_name}/dna/variants/clone_specific_beds/{dna_sample_name}-clone_{clone}_SNVs_filtered.bed",
        bam="data/{rna_sample_name}/rna/bam/{rna_sample_name}_cell_{barcode}.bam",
        ref = ref_genome
    output:
        "data/{rna_sample_name}/rna/variants/bcf_by_clone/{rna_sample_name}-cell_{barcode}_variants_for_clone_{clone}_in_{dna_sample_name}.bcf"
    shell:
        "samtools mpileup -l {input.clone_bedfile} -R -g -f {input.ref} {input.bam} > {output}"

rule rna_call_variants:
    input:
        "data/{rna_sample_name}/rna/variants/bcf_by_clone/{rna_sample_name}-cell_{barcode}_variants_for_clone_{clone}_in_{dna_sample_name}.bcf"
    output:
        "data/{rna_sample_name}/rna/variants/vcf_by_clone/{rna_sample_name}-cell_{barcode}_variants_for_clone_{clone}_in_{dna_sample_name}.vcf"
    shell:
        "bcftools call -c -v {input} > {output}"

rule rna_collate_snvs:
    input:
        "data/{rna_sample_name}/rna/variants/vcf_by_clone/{rna_sample_name}-cell_{barcode}_variants_for_clone_{clone}_in_{dna_sample_name}.vcf"
    output:
        "data/{rna_sample_name}/rna/variants/logs/{rna_sample_name}-cell_{barcode}_variants_for_clone_{clone}_in_{dna_sample_name}.txt"
    shell:
        """
        Rscript scripts/parse_vcf.R --input_vcf {input} --cell_barcode {wildcards.barcode} --clone {wildcards.clone} >> data/{wildcards.rna_sample_name}/rna/variants/cell_snps.csv
        touch {output}
        """




# rule rna_pileup:
#     input: 
#         bamfile = "data/rna/bam/{sample_name}_cell_{barcode}.bam",
#         ref = ref_genome
#     output:
#         "data/rna/variants/bcf/{sample_name}-cell_{barcode}_variants.bcf"
#     shell:
#         "samtools mpileup -R -g -f {input.ref} {input.bamfile} > {output}"

# rule rna_bcf_to_vcf:
#     input:
#         "data/rna/variants/bcf/{sample_name}-cell_{barcode}_variants.bcf"
#     output:
#         "data/rna/variants/vcf/{sample_name}-cell_{barcode}_variants.vcf"
#     shell:
#         "bcftools call -c -v {input} > {output}"



# rule split_rna_sam_on_barcode:
#     input:
#         "data/rna/bam_from_10x/{}_possorted_genome_bam.bam".format(rna_sample_name)
#     output:
#         rna_cell_bams
#     shell:
#         """
#         samtools view data/rna/bam_from_10x/SA501X2B_possorted_genome_bam.bam | head  | awk -F "\t" 'match($0, /CB:Z:[ACGT]-1/) {if (RSTART>0){OUTPUT=substr($0,RSTART+5,18); print $0 >> "data/rna/bam/"OUTPUT".sam"}}'
#         samtools view -H data/rna/bam_from_10x/SA501X2B_possorted_genome_bam.bam > tmp_header.sam
#         for f in data/rna/sam/*.sam
#         do
#             b=`basename $f .sam`
#             cat tmp_header.sam $f.sam | samtools view -b - > $b.bam
#         done
#         """

# rule rna_bam_to_fastq:
#     input:
#         "data/rna/bam/{sample_name}_cell_{barcode}.bam"
#     output:
#         "data/rna/fastq/{sample_name}_cell_{barcode}.fastq"
#     shell:
#         "bedtools bamtofastq -i {input} -fq {output}"

# rule star_index:
#     input:
#         ref_genome
#     output:
#         "data/rna/star_index/genomeParameters.txt"
#     shell:
#         "STAR --runMode genomeGenerate --genomeDir data/rna/star_index --genomeFastaFiles {input}"

# rule rna_realign:
#     input:
#         fastq = "data/rna/fastq/{sample_name}_cell_{barcode}.fastq",
#         star_index = "data/rna/star_index/genomeParameters.txt"
#     output:
#         "data/rna/bam_realigned/{sample_name}_cell_{barcode}.bam"
#     shell:
#         "STAR --genomeDir data/rna/star_index --readFilesIn {input.fastq}"


