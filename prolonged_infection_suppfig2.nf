#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// NOTE that helpful runtime parameters are in the file `nextflow.config`
// We recommend that you modify any file paths, inputs, or outputs there.


// Workflow specification
workflow {

	ch_patient_counts = Channel
		.fromPath( params.patient_variants )

	ch_palette = Channel
		.fromPath( params.color_palette )
	
	ch_include_check = Channel
		.fromPath( params.include )
		.ifEmpty( 'Accessions not yet selected' )
	
	
	PULL_METADATA (	
		ch_include_check
	)
	
	REFORMAT_METADATA (
		ch_include_check,
		PULL_METADATA.out
	)
	
	SELECT_SUBSAMPLE (
		ch_include_check,
		REFORMAT_METADATA.out
	)
	
	ARRANGE_INCLUDE_LIST ( )
	
	PULL_FASTAS (
		ARRANGE_INCLUDE_LIST.out
			.splitCsv ( header: true )
			.map {row -> tuple(row.accession, row.date, row.pango)}
	)
	
	SUBSAMPLE_ALIGNMENT (
		PULL_FASTAS.out
	)
	
	SUBSAMPLE_VARIANT_CALLING (
		SUBSAMPLE_ALIGNMENT.out
	)

	SUPP_FIGURE_2_PLOTTING (
		SUBSAMPLE_VARIANT_CALLING.out.collect(),
		ch_patient_counts,
		ch_palette
	)

}


// Defining each process in the workflow
process PULL_METADATA {
	
	when:
	include == 'Accessions not yet selected'
	
	input:
	val(include_check)
	
	output:
	path("sarscov2-metadata.jsonl")
	
	script:
	"""
	datasets summary virus genome taxon sars-cov-2 \
	--complete-only \
	--as-json-lines \
	> sarscov2-metadata.jsonl
	"""
	
}


process REFORMAT_METADATA {
	
	publishDir params.results_data_files
	
	when:
	include == 'Accessions not yet selected'
	
	input:
	val(include_check)
	path(jsonl)
	
	output:
	path("*.tsv")
	
	script:
	"""
	dataformat tsv virus-genome \
	--fields accession,bioprojects,biosample-acc,geo-location,geo-region,host-organism-name,isolate-collection-date,length,nucleotide-completeness,release-date,sourcedb,sra-accs,virus-pangolin,virus-strain,virus-tax-id --inputfile ${jsonl} > sarscov2-metadata.tsv
	"""
	
}


process SELECT_SUBSAMPLE {
	
	publishDir params.refdir, mode: 'copy'
	
	when:
	include == 'Accessions not yet selected'
	
	input:
	val(include_check)
	path(tsv)
	
	output:
	path("include_list.csv")
	
	script:
	"""
	select_subsample.R ${tsv} ${params.subsample_size} \
	${params.min_date} ${params.max_date}
	"""
	
}


process ARRANGE_INCLUDE_LIST { 
	
	output:
	path("include_list.csv")
	
	script:
	"""
	cp ${params.include} .
	"""
	
}


process PULL_FASTAS {
	
	tag "${accession}"
	
	cpus 1
	
	input:
	tuple val(accession), val(date), val(pango)
	
	output:
	tuple val(accession), path("*.fasta")
	
	script:
	"""
	datasets download virus genome accession "${accession}" \
	--exclude-cds --exclude-protein
	unzip ncbi_dataset.zip
	mv ncbi_dataset/data/genomic.fna ./"${accession}".fasta
	rm -rf ncbi_dataset/
	"""
	
}


process SUBSAMPLE_ALIGNMENT {

	// Aligning consensus sequences so that they are in sam format

	tag "${accession}"

	input:
	tuple val(accession), file(fasta)

	output:
	tuple val(accession), path("*.sam")

	script:
	"""
	minimap2 -a ${params.refseq} ${fasta} > ${accession}.sam
	"""

}


process SUBSAMPLE_VARIANT_CALLING {

	// calling single-nucleotide variants from each genbank sequence

	tag "${accession}"

	publishDir params.results_data_files, mode: "copy"

	cpus 4

	input:
	tuple val(accession), path(sam)

	output:
	path("*.vcf")

	script:
	"""

	callvariants.sh \
	in=${sam} out=${accession}.vcf.gz \
	ref=${params.refseq} samstreamer=t ss=4 clearfilters \
	ploidy=1 mincov=0 callsub=t calldel=f callins=f overwrite=t

	gunzip ${accession}.vcf.gz

	"""

}


process SUPP_FIGURE_2_PLOTTING {

	// counts variants for each genbank sequence and plots the figure

	publishDir params.visuals, pattern: "*.pdf", mode: "move"
	publishDir params.results, pattern: "*.csv", mode: "move"

	input:
	path(vcf_list)
	file(mutations_counts)
	path(palette)

	output:
	file("*.pdf")
	file("*.csv")

	script:
	"""
	SupplementalFigure2_global_roottotip_plot.R "." "${params.include}"
	"""

}


