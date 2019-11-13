#!/usr/bin/env nextflow
/*
========================================================================================
                                PHW PrepareRef Pipeline
========================================================================================
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    Usage:
    nextflow run connor-lab/prepareref_nexflow [ options ] --outdir [output_directory] --scheme [pubmlst_scheme] ( --taxid [NCBI TaxID] | --fasta_dir [fasta_directory] )
    
    Mandatory arguments:
      --taxid                       NCBI TaxID. Used to download complete- and chromosomal-level assemblies
        AND/OR 
      --chromosome_dir              Path to directory containing complete genome fasta files
        AND/OR
      --contig_dir                  Path to directory containing contig-level assembly fasta files

      --scheme                      PubMLST scheme name 
      
      --outdir                      The output directory where the results will be saved
    
    Other options:      
      --ncbi_section                NCBI section to download genomes from
                                    Allowed values: ${params.ncbi_section_av}
                                    Default: ${params.ncbi_section}
      
      --minstsize                   Minimum number of genomes in a ST to retain the ST in final dataset
                                    Allowed values: >3
                                    Default: ${params.minstsize}

      --maxstsize                   Maximum number of genomes in a ST to retain in final dataset
                                    Default: ${params.maxstsize}

      --mapper                      PHEnix mapper
                                    Allowed values: ${params.mapper_av}
                                    Default: ${params.mapper}

      --variant                     PHEnix variant caller
                                    Allowed values: ${params.variant_av}
                                    Default: ${params.variant}

      --force                       Force overwrite of output directory
                                    Default: ${params.force}
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}


// Required options, at least one must be set
chromosome_dir = params.chromosome_dir
contig_dir = params.contig_dir
taxid = params.taxid
if ( !chromosome_dir && !contig_dir && !taxid ){
    println("Set either \"--taxid\" or \"--contig_dir\" \"--chromosome_dir\" or any combination!\n")
    helpMessage()
    exit 1
}


// Output directory
if (params.outdir) {
    // Trim trailing slash
    outputDir = "${params.outdir}".replace(/\/$/, "")
    
    // Exit if output dir exists and --force not set, otherwise delete output dir
    def folder = new File("$outputDir")
    if( folder.exists() && !params.force ) {
        println("Output directory already exists and \"--force\" not set!\n")
        helpMessage()
        exit 1
    } else if ( folder.exists() && params.force ) {
        folder.deleteDir()
        println("Output directory already exists and \"--force\" set - DELETING!\n")
    }

} else {
    println("Please supply an output directory\n")
    helpMessage()
    exit 1
}

// PubMLST scheme
scheme = params.scheme
if ( !params.scheme ) {
    println("Please choose a MLST scheme supported by https://github.com/tseemann/mlst\n")
    helpMessage()
    exit 1
}

// NCBI section for ncbi-genome-download
ncbi_section = params.ncbi_section
if ( !params.ncbi_section_av.contains(ncbi_section) ) {
    println("Please choose a supported NCBI section\n")
    helpMessage()
    exit 1
} 

// Minimum size of an ST
minstsize = params.minstsize
if (minstsize < 3 ){
    println("Please choose a minstsize > 3\n")
    helpMessage()
    exit 1
}

// Maximum ST size
maxstsize = params.maxstsize

// PHEnix mapper
mapper = params.mapper

if ( !params.mapper_av.contains(mapper) ){
    println("Please choose a supported mapper\n")
    helpMessage()
    exit 1
}

//PHEnix variant
variant = params.variant
if ( !params.variant_av.contains(variant) ){
    println("Please choose a supported variant caller\n")
    helpMessage()
    exit 1
}


if ( taxid ) {
    Channel.from( "${taxid}" )
           .into{ ch_taxid_getCompleteSeqs ; ch_taxid_getContigSeqs }
}

if ( chromosome_dir ) {
    Channel.fromPath( "${chromosome_dir}/*.{fas,fasta,fa,fsa,fna}" )
           .map{ tuple( it, "COMPLETE" ) }
           .set{ ch_chromosome_dir }
} 

if ( contig_dir ) {
    Channel.fromPath( "${contig_dir}/*.{fas,fasta,fa,fsa,fna}" )
           .map{ tuple( it, "CONTIG" ) }
           .set{ ch_contig_dir }
}



if ( taxid ) {
    process GETCOMPLETESEQS_NCBIGENOMEDOWNLOAD {
        tag "${taxonid}"

        publishDir "${outputDir}/downloaded_genomes", mode: 'copy', pattern: "*.fna"

        cpus 2

        input:
        val taxonid from ch_taxid_getCompleteSeqs

        output:
        tuple file("*.fna"), val("COMPLETE") into ch_getCompleteSeqs

        script:
        """
        ncbi-genome-download -F fasta \
        -l complete,chromosome \
        -p ${task.cpus} \
        -m TaxID_${taxonid}.genome_download.tab \
        -s ${ncbi_section}\
        -t ${taxonid} \
        bacteria
        
        for seqfile in \$(find ${ncbi_section} -name "*.fna.gz"); do
            acc=\$(echo "\${seqfile}" | rev | cut -f2 -d "/" | rev)
            zcat \${seqfile} > \${acc}.fna
        done
        """
    }

    process GETCONTIGSEQS_NCBIGENOMEDOWNLOAD {
        tag "${taxonid}"

        publishDir "${outputDir}/downloaded_genomes", mode: 'copy', pattern: "*.fna"

        cpus 2

        input:
        val taxonid from ch_taxid_getContigSeqs

        output:
        tuple file("*.fna"), val("CONTIG") into ch_getContigSeqs

        script:
        """
        ncbi-genome-download -F fasta \
        -l scaffold,contig \
        -p ${task.cpus} \
        -m TaxID_${taxonid}.genome_download.tab \
        -s ${ncbi_section}\
        -t ${taxonid} \
        bacteria
        
        for seqfile in \$(find ${ncbi_section} -name "*.fna.gz"); do
            acc=\$(echo "\${seqfile}" | rev | cut -f2 -d "/" | rev)
            zcat \${seqfile} > \${acc}.fna
        done
        """
    }
}

if ( taxid ) {
    ch_getContigSeqs.transpose()
                    .concat(ch_getCompleteSeqs.transpose())
                    .branch {
                        complete: it[1] == "COMPLETE"
                        contig: it[1] == "CONTIG"
                    }
                    .set{ ch_downloadedFasta }
} else {
    Channel.empty()
           .branch {
               complete: it[1] == "COMPLETE"
               contig: it[1] == "CONTIG"
               }
               .set{ ch_downloadedFasta }
}

if ( contig_dir && chromosome_dir ) {
    ch_downloadedFasta.complete.concat( ch_chromosome_dir )
                               .concat( ch_contig_dir )
                               .concat( ch_downloadedFasta.contig )
                               .branch {
                                   complete: it[1] == "COMPLETE"
                                   contig: it[1] == "CONTIG"
                                   }
                               .set{ ch_allFastaSeqs }
} else if ( !contig_dir && chromosome_dir ) {
    ch_downloadedFasta.complete.concat( ch_chromosome_dir )
                               .concat( ch_downloadedFasta.contig )
                               .branch {
                                   complete: it[1] == "COMPLETE"
                                   contig: it[1] == "CONTIG"
                                   }
                               .set{ ch_allFastaSeqs }
} else if ( contig_dir && !chromosome_dir ) {
    ch_downloadedFasta.complete.concat( ch_contig_dir )
                               .concat( ch_downloadedFasta.contig )
                               .branch {
                                   complete: it[1] == "COMPLETE"
                                   contig: it[1] == "CONTIG"
                                   }
                               .set{ ch_allFastaSeqs }
} else if ( !contig_dir && !chromosome_dir ) {
    ch_downloadedFasta.complete.concat( ch_downloadedFasta.contig )
                               .branch {
                                   complete: it[1] == "COMPLETE"
                                   contig: it[1] == "CONTIG"
                               }
                               .set{ ch_allFastaSeqs }
} else { 
    error "Error combining local and downloaded genome sequences"
}


process REMOVEEXTRACHROMOSOMALSEQS_BIOPYTHON {
    tag "${genome_fasta}"

    input:
    tuple file(genome_fasta), assembly_status from ch_allFastaSeqs.complete

    output:
    tuple file("chr_${genome_fasta}"), assembly_status into ch_removeExtrachromosomalSeqs_callSequenceType
    file("chr_${genome_fasta}") into ch_removeExtrachromosomalSeqs_clusterGenomes
    file("chr_${genome_fasta}") into ch_removeExtrachromosomalSeqs_mapGenomes

    script:
    """
    #!/usr/bin/env python3

    import re

    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    max_len = 0

    for seq_record in SeqIO.parse("${genome_fasta}", "fasta"):
        if len(seq_record) > max_len:
            max_len = len(seq_record)
            max_seq_record_str = re.sub("[^ACGTN]", "N", str(seq_record.seq).upper())

            output_seq = SeqRecord(Seq(max_seq_record_str), id=seq_record.id, description=seq_record.description)

    SeqIO.write(output_seq, "chr_${genome_fasta}", "fasta")
    """
}

process CALLSEQUENCETYPE_MLST {
    tag "${genome_fasta}"

    input:
    tuple file(genome_fasta), assembly_status from ch_removeExtrachromosomalSeqs_callSequenceType.concat(ch_allFastaSeqs.contig)

    output:
    tuple assembly_status, stdout into ch_callSequenceType1
    file("${genome_fasta_basename}_mlst.tab") into ch_callSequenceType_sequenceTypeSummary
    file("${genome_fasta}") into ch_callSequenceType2

    script:
    genome_fasta_basename = genome_fasta.getBaseName()
    """
    mlst --scheme ${scheme} --nopath ${genome_fasta} > ${genome_fasta_basename}_mlst.tab
    grep "${genome_fasta}" ${genome_fasta_basename}_mlst.tab
    """
}

process SEQUENCETYPESUMMARY_LINUX {

    publishDir "${outputDir}", mode: 'copy', pattern: "mlst_summary.tab"

    input:
    file("*") from ch_callSequenceType_sequenceTypeSummary.collect()

    output:
    file("mlst_summary.tab")
        
    script:
    """
    grep -h SCHEME *.tab | uniq > header.txt
    grep -h "${scheme}" *.tab > body.txt
    cat header.txt body.txt > mlst_summary.tab
    """
}


/* This filters 
    - uncalled MLST STs
    - STs without a COMPLETE genome
    - STs smaller than $params.minstsize
    - Genomes from STs larger than $params.maxstsize (sorts first for determinism) */


ch_callSequenceType1.map{ tuple( it[0], it[1].trim().split("\t") )}
                   .map{ tuple(it[1][2], tuple( it[1][0] , it[0] ) ) }
                   .filter{ !( it[0] =~ /-/ ) }
                   .groupTuple()
                   .map{ [ it[0], it[1].take( maxstsize.toInteger() ) ] }
                   .filter{ it[1].size() >= minstsize }
                   .filter{ it[1][0][1] ==~ /COMPLETE/ }
                   .transpose()
                   .map{ tuple( it[1][0], it[0], it[1][1] )}
                   .join( ch_removeExtrachromosomalSeqs_clusterGenomes.concat(ch_callSequenceType2)
                        .map{ tuple( it.getName(), it ) }, by: 0 )
                   .map{ tuple( it[1], it[2], it[0], it[3] ) }
                   .into { ch_callSequenceType_clusterGenomes ; 
                           ch_callSequenceType_mapGenomes ;
                           ch_callSequenceType_copyGenomes }
                    //[ ST, assembly_status, seq_name, seq_path ]


process COPYGENOMES_LOCAL {
    tag "${name}"

    publishDir "${outputDir}/ST${ST}/fasta/${assembly_status}", mode: 'copy', pattern: "${fasta}"

    input:
    tuple ST, assembly_status, name, file(fasta) from ch_callSequenceType_copyGenomes

    output:
    file("${fasta}")

    script:
    """
    sleep 1
    """
}



process CLUSTERGENOMES_MASH {
    tag "ST${ST}"

    publishDir "${outputDir}/ST${ST}", mode: 'copy', pattern: "ST${ST}.mash_dist.tab"

    input:
    tuple ST, file("fasta/*") from ch_callSequenceType_clusterGenomes.filter{ it[1] ==~ /COMPLETE/ }.groupTuple().map{ tuple( it[0], it[3] ) }

    output:
    tuple ST, file("ST${ST}.mash_dist.tab") into ch_clusterGenomes_identifyRefSequence
    
    script:
    """
    mash sketch -o all.msh -s 5000 fasta/*

    mash dist all.msh all.msh | sed 's/fasta\\///g' > ST${ST}.mash_dist.tab
    """
}

process IDENTIFYSTREFSEQUENCE_PYTHON {
    tag "ST${ST}"

    input:
    tuple ST, file(mash_dist_tab) from ch_clusterGenomes_identifyRefSequence
    
    output:
    tuple ST, stdout into ch_identifyRefSequence_mapGenomes

    script:
    """
    #!/usr/bin/env python3
    import csv 
    import statistics

    distances = {}

    with open("${mash_dist_tab}") as tabfile:
        lines = csv.reader(tabfile, delimiter = "	")
        for line in lines:
            item = distances.get(line[0], dict())
            item[line[1]] = float(line[2])

            distances[line[0]] = item


    min_dist_med = 2

    for source, dest in distances.items():
        median_dist = statistics.median(list(dest.values()))

        if median_dist < min_dist_med:
            min_dist_med = median_dist
            lowest_dist_med = source

    print(lowest_dist_med, end = '')
    """
}


/* Use reference genome selection to link 
the reference genome filepath with its ST */


ch_identifyRefSequence_mapGenomes.map{ [ it[1], it[0] ] }
                                 .join(ch_removeExtrachromosomalSeqs_mapGenomes
                                        .map{ [it.getName() , it ] }, by: 0)
                                 .map{ [ it[1], it[2] ] }
                                 .into{ ch_mapGenomes ;
                                        ch_renameRef ;
                                        ch_makeAlignment ; 
                                        ch_findHgtSequence ; 
                                        ch_maskReferenceRecomb ; 
                                        ch_calculateNoSimReads }
                                 // [ ST, ST_ref_genome_file ]

process RENAMEREF_LOCAL {
    tag "ST${ST}-${ref_fasta}"

    input:
    tuple ST, file(ref_fasta) from ch_renameRef

    output:
    file("ST${ST}.fasta") into ch_renameRef_mashSketchRefs

    script:
    """
    mv ${ref_fasta} ST${ST}.fasta
    """
}

process MASHSKETCHREFS_MASH {

    publishDir "${outputDir}", mode: 'copy', pattern: "${scheme}.msh"

    input:
    file("*") from ch_renameRef_mashSketchRefs.collect()

    output:
    file("${scheme}.msh")

    script:
    """
    mash sketch -s 5000 -o ${scheme}.msh *
    """
}


process MAPGENOMES_SNIPPY {
    tag "ST${ST}-${query}"

    cpus 2

    input:
    tuple ST, file(query), file("ref*.fa") from ch_callSequenceType_mapGenomes.map{ [ it[0], it[3] ] }.combine( ch_mapGenomes, by: 0)

    output:
    tuple ST, file("${query}_snippy") into ch_mapGenomes_makeAlignment
    
    script:
    """
    snippy --cpus ${task.cpus} --cleanup --outdir ${query}_snippy --ref ref.fa --ctgs ${query}
    """
}

process MAKEALIGNMENT_SNIPPYCORE {
    tag "ST${ST}"

    publishDir "${outputDir}/ST${ST}", mode: 'copy', pattern: "ST${ST}*"

    input:
    tuple ST, file("*"), file(ref) from ch_mapGenomes_makeAlignment.groupTuple().combine(ch_makeAlignment, by: 0)

    output:
    tuple ST, file("ST${ST}.full.aln") into ch_makeAlignment_detectRecombination
    tuple ST, file("ST${ST}.full.aln") into ch_makeAlignment_calculatePhylogeny
    file("ST${ST}*")


    script:
    """
    snippy-core --prefix ST${ST} --ref ${ref} *_snippy
    """
}


process CALCULATEPHYLOGENY_IQTREE {
    tag "ST${ST}"

    publishDir "${outputDir}/ST${ST}", mode: 'copy', pattern: "ST${ST}.treefile"

    cpus 2

    input:
    tuple ST, file(aln) from ch_makeAlignment_calculatePhylogeny

    output:
    tuple ST, file("ST${ST}.treefile") into ch_calculatePhylogeny_detectRecombination
    
    
    script:
    """ 
    iqtree -nt ${task.cpus} -s ${aln} -bb 1000 -t PARS -ninit 2 -m GTR -pre ST${ST}
    """
}

process DETECTRECOMBINATION_GUBBINS {
    tag "ST${ST}"

    publishDir "${outputDir}/ST${ST}/gubbins_recombination_predictions", mode: 'copy', pattern: "ST${ST}.*"

    cpus 2

    input:
    tuple ST, file(tree), file(full_aln) from ch_calculatePhylogeny_detectRecombination.join(ch_makeAlignment_detectRecombination)

    output:
    tuple ST, file("ST${ST}.recombination_predictions.embl") into ch_detectRecombination_mergeHgtRegions
    file("ST${ST}.*")

    script:
    """
    run_gubbins.py --starting_tree ${tree} --threads ${task.cpus} --prefix ST${ST} ${full_aln}
    """
}

process FINDHGTSEQUENCE_ALIENHUNTER {
    tag "ST${ST}"

    publishDir "${outputDir}/ST${ST}/alienhunter_recombination_predictions", mode: 'copy', pattern: "ST${ST}_alienhunter.embl"

    input:
    tuple ST, file(fasta) from ch_findHgtSequence

    output:
    tuple ST, file("ST${ST}_alienhunter.embl") into ch_findHgtSequence_mergeHgtRegions

    script:
    """
    alien_hunter ${fasta} ST${ST}_alienhunter.embl
    """
}

process MERGEHGTREGIONS_PYTHON {
    tag "ST${ST}"

    publishDir "${outputDir}/ST${ST}", mode: 'copy', pattern: "ST${ST}.recombination.bed"
    
    input:
    tuple ST, file(gubbins_embl), file(ah_embl) from ch_detectRecombination_mergeHgtRegions.join( ch_findHgtSequence_mergeHgtRegions )

    output:
    tuple ST, file("ST${ST}.recombination.bed") into ch_mergeHgtRegions_maskReferenceRecomb

    script:
    """
    #!/usr/bin/env python3
    
    files = { "gubbins": "${gubbins_embl}", "alienhunter": "${ah_embl}" }
              
    features = []

    for file, path in files.items():
        with open(path, "r") as f:
            for line in f:
                if "misc_feature" in line:
                    feature = {"source" : file, 
                                "start" : int(line.split()[2].split("..")[0]), 
                                  "end" : int(line.split()[2].split("..")[1]) }
                    features.append(feature)

    features_sorted = sorted(features, key = lambda i: i['start']) 


    with open("ST${ST}.recombination.bed", 'w') as outfile:
        for feature in features_sorted:
            line = "SEQUENCE\\t{start!s}\\t{end!s}\\t{source}\\n".format(**feature)
            outfile.write(line)
    """
}

process MASKREFERENCERECOMB_BEDTOOLS {
    tag "ST${ST}-${fasta}"

    publishDir "${outputDir}/ST${ST}", mode: 'copy', pattern: "mask_${fasta}"

    input:
    tuple ST, file(bedfile), file(fasta) from ch_mergeHgtRegions_maskReferenceRecomb.join(ch_maskReferenceRecomb)

    output:
    tuple ST, file("mask_${fasta}") into ch_maskReferenceRecomb_makeSnapperDbFiles

    script:
    """
    head -n1 ${fasta} > tmp.header
    sed '1s/^.*\$/>SEQUENCE/' ${fasta} > tmp.fasta
 
    bedtools maskfasta -fi tmp.fasta -fo tmp1.fasta -bed ${bedfile}
    
    grep -v "SEQUENCE" tmp1.fasta > tmp2.fasta
    cat tmp.header tmp2.fasta > mask_${fasta}
    """
}

process CALCULATENOSIMREADS_PYTHON {

    input:
    tuple ST, file(ref_fasta) from ch_calculateNoSimReads

    output:
    tuple ST, file("${ref_fasta}"), stdout into ch_calculateNoSimReads_simulateFastq

    script:
    """
    #!/usr/bin/env python3

    from Bio import SeqIO

    readlength = int(${params.readlength})
    coverage = int(${params.coverage})

    for seq_record in SeqIO.parse("${ref_fasta}", "fasta"):
        no_reads = int( ( len(seq_record) * coverage ) / ( readlength * 2 ) )

    print(no_reads, end = '')
    """

}

process SIMULATEREFFASTQS_WGSIM {
    tag "ST${ST}-${ref_fasta}"

    publishDir "${outputDir}/ST${ST}", mode: 'copy', pattern: "${ref_basename}.R{1,2}.fq.gz"

    input:
    tuple ST, file(ref_fasta), no_reads from ch_calculateNoSimReads_simulateFastq
    
    output:
    tuple ST, file("${ref_fasta}"), file("${ref_basename}.R1.fq.gz"), file("${ref_basename}.R2.fq.gz") into ch_simulateFastqs_makeSnapperDbJson

    script:
    ref_basename = ref_fasta.getBaseName()
    """
    wgsim -d 600 -e 0 -r 0 -R 0 -X 0 \
    -N ${no_reads} \
    -1 ${params.readlength} -2 ${params.readlength} \
    ${ref_fasta} \
    ${ref_basename}.R1.fq \
    ${ref_basename}.R2.fq
    
    gzip *.fq
    """
}

process MAKESNAPPERDBFILES_BWAPICARDSAMTOOLS {
    tag "ST${ST}-${ref_fasta}"

    publishDir "${outputDir}/ST${ST}", mode: 'copy', pattern: "${ref_basename}*{fai,dict,bt2}"

    input:
    tuple ST, file(ref_fasta) from ch_maskReferenceRecomb_makeSnapperDbFiles

    output:
    tuple ST, file("${ref_fasta}"), file("${ref_basename}*{fai,dict,bt2,amb,ann,bwt,pac,sa}") into ch_makeSnapperDbFiles_makeSnapperDbJson
    
    script:
    ref_basename = ref_fasta.getBaseName()

    if ( mapper == 'bwa' && variant == 'gatk' )
        """
        samtools faidx ${ref_fasta}
        picard CreateSequenceDictionary R=${ref_fasta} O=${ref_basename}.dict URI=file:/${ref_fasta}
        bwa index ${ref_fasta}
        """
    
    else if ( mapper == 'minimap2' && variant == 'gatk' )
        """
        samtools faidx ${ref_fasta}
        picard CreateSequenceDictionary R=${ref_fasta} O=${ref_basename}.dict URI=file:/${ref_fasta}
        """

    else if ( mapper == 'bowtie2' && variant == 'gatk' )
        """
        samtools faidx ${ref_fasta}
        picard CreateSequenceDictionary R=${ref_fasta} O=${ref_basename}.dict URI=file:/${ref_fasta}
        bowtie2-build ${ref_fasta} ${ref_fasta}
        """

    else if ( mapper == 'bwa' && variant == 'mpileup' )
        """
        samtools faidx ${ref_fasta}
        bwa index ${ref_fasta}
        """

    else if ( mapper == 'minimap2' && variant == 'mpileup' )
        """
        samtools faidx ${ref_fasta}
        """

    else if ( mapper == 'bowtie2' && variant == 'mpileup' )
        """
        samtools faidx ${ref_fasta}
        bowtie2-build ${ref_fasta} ${ref_fasta}
        """
    
    else
        error "Invalid mapper/variant-caller combination (mapper: ${mapper} | variant: ${variant}"
}

process MAKESNAPPERDDBJSON_PHENIX {
    tag "ST${ST}-${ref_fasta}"

    publishDir "${outputDir}/ST${ST}", mode: 'copy', pattern: "snapperdb_${ref_basename}"

    input:
    tuple ST, 
    file(forward), 
    file(reverse), 
    file(ref_fasta), 
    file("*") from ch_simulateFastqs_makeSnapperDbJson.map{ [ it[0], it[2], it[3] ] }
                                                      .join( ch_makeSnapperDbFiles_makeSnapperDbJson )

    output:
    file("snapperdb_${ref_basename}")

    script:
    ref_basename = ref_fasta.getBaseName()
    """
    export TMPDIR=\$(pwd)

    phenix.py run_snp_pipeline \
    --json \
    --json-info \
    --annotators coverage \
    --filters mq_score:30,min_depth:10,ad_ratio:0.9 \
    --mapper ${mapper} \
    --variant ${variant} \
    --sample-name ${ref_basename} \
    --reference ${ref_fasta} \
    -r1 ${forward} \
    -r2 ${reverse} \
    --outdir snapperdb_${ref_basename}
    """
}
