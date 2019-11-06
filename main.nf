
ncbi_section = params.ncbi_section

fasta_dir = params.fasta_dir

taxid = params.taxid

minstsize = params.minstsize

scheme = params.scheme

if ( taxid ) {
    Channel.from( "${taxid}" )
           .set{ ch_taxid_getSeqs }
}

if ( fasta_dir ) {
    Channel.fromPath( "${fasta_dir}/*.{fas,fasta,fa,fsa,fna}" )
           .set{ ch_fasta_dir_removeExtrachromosomalSeqs }
}


if ( taxid ) {
    process GETSEQS_NCBIGENOMEDOWNLOAD {
        input:
        val taxonid from ch_taxid_getSeqs

        output:
        file("*.fna") into ch_getSeqs_removeExtrachromosomalSeqs

        script:
        """
        ncbi-genome-download -s ${ncbi_section} -F fasta -l complete,chromosome -t ${taxonid} -p ${task.cpus} bacteria
        for seqfile in \$(find ${ncbi_section} -name "*.fna.gz"); do
            acc = \$(echo "\${seqfile}" | rev | cut -f2 -d "/" | rev)
            zcat \${seqfile} > \${acc}.fna
        """
    }
}

if ( taxid && fasta_dir ) {
    ch_getSeqs_removeExtrachromosomalSeqs.flatten()
                                         .concat( ch_fasta_dir_removeExtrachromosomalSeqs )
                                         .unique()
                                         .set{ ch_removeExtrachromosomalSeqs }
} else if ( taxid && !fasta_dir ) {
    ch_getSeqs_removeExtrachromosomalSeqs.flatten()
                                         .set{ ch_removeExtrachromosomalSeqs }
} else if ( fasta_dir && !taxid ) {
    ch_fasta_dir_removeExtrachromosomalSeqs.set{ ch_removeExtrachromosomalSeqs }
}

process REMOVEEXTRACHROMOSOMALSEQS_BIOPYTHON {
    input:
    file(genome_fasta) from ch_removeExtrachromosomalSeqs

    output:
    file("chr_${genome_fasta}") into ch_removeExtrachromosomalSeqs_callSequenceType
    file("chr_${genome_fasta}") into ch_removeExtrachromosomalSeqs_clusterGenomes
    file("chr_${genome_fasta}") into ch_removeExtrachromosomalSeqs_mapGenomes

    script:
    """
    #!/usr/bin/env python3

    from Bio import SeqIO

    max_len = 0

    for seq_record in SeqIO.parse("${genome_fasta}", "fasta"):
        if len(seq_record) > max_len:
            max_len = len(seq_record)
            max_seq_record = seq_record

    SeqIO.write(max_seq_record, "chr_${genome_fasta}", "fasta")
    """
}

process CALLSEQUENCETYPE_MLST {
    input:
    file(genome_fasta) from ch_removeExtrachromosomalSeqs_callSequenceType

    output:
    stdout into ch_callSequenceType 
    
    script:
    """
    mlst --scheme ${scheme} --nopath ${genome_fasta} | grep "${genome_fasta}"
    """
}

/*
This filters uncalled MLST STs and filters STs with fewer 
than $params.minstsize representatives
*/

ch_callSequenceType.splitCsv(sep: "\t").map{ [ it[0], it [2] ] }
                   .filter{ !( it[1] =~ /-/ ) }
                   .join( ch_removeExtrachromosomalSeqs_clusterGenomes.map{ [ it.getName(), it ] } )
                   .map{ [ it[1], it[0], it[2] ] }
                   .groupTuple()
                   .filter{ it[1].size() >= minstsize }
                   .transpose()
                   .into { ch_callSequenceType_clusterGenomes ; ch_callSequenceType_mapGenomes }
                   // [ ST, seq_name, seq_path ]

           
process CLUSTERGENOMES_MASH {
    input:
    tuple ST, file("fasta/*") from ch_callSequenceType_clusterGenomes.groupTuple().map{ [ it[0], it[2] ] }

    output:
    tuple ST, file("${ST}.mash_dist.tab") into ch_clusterGenomes_identifyRefSequence
    
    script:
    """
    mash sketch -o all.msh -s 5000 fasta/*

    mash dist all.msh all.msh | sed 's/fasta\\///g' > ${ST}.mash_dist.tab
    """
}

process IDENTIFYREFSEQUENCE_PYTHON {
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


ch_identifyRefSequence_mapGenomes.map{ [ it[1], it[0] ] }
                                 .join(ch_removeExtrachromosomalSeqs_mapGenomes.map{ [it.getName() , it ] }, by: 0)
                                 .map{ [ it[1], it[2] ] }
                                 .into{ ch_mapGenomes ; ch_makeAlignment ; ch_findHgtSequence ; ch_maskReferenceRecomb ; ch_simulateFastq}
                                 // [ ST, ref_genome_file ]

process MAPGENOMES_SNIPPY {

    cpus 2

    input:
    tuple ST, file(query), file("ref*.fa") from ch_callSequenceType_mapGenomes.map{ [ it[0], it[2] ] }.combine( ch_mapGenomes, by: 0)

    output:
    tuple ST, file("${query}_snippy") into ch_mapGenomes_makeAlignment
    
    script:
    """
    snippy --cpus ${task.cpus} --cleanup --outdir ${query}_snippy --ref ref.fa --ctgs ${query}
    """
}

process MAKEALIGNMENT_SNIPPYCORE {
    input:
    tuple ST, file("*"), file(ref) from ch_mapGenomes_makeAlignment.groupTuple().combine(ch_makeAlignment, by: 0)

    output:
    tuple ST, file("ST${ST}.full.aln") into ch_makeAlignment_detectRecombination
    tuple ST, file("ST${ST}.aln") into ch_makeAlignment_calculatePhylogeny


    script:
    """
    snippy-core --prefix ST${ST} --ref ${ref} *_snippy
    """
}


process CALCULATEPHYLOGENY_IQTREE {
    cpus 2

    input:
    tuple ST, file(aln) from ch_makeAlignment_calculatePhylogeny

    output:
    tuple ST, file("ST${ST}.treefile") into ch_calculatePhylogeny_detectRecombination
    
    
    script:
    """ 
    iqtree -nt ${task.cpus} -s ${aln} -bb 1000 -m GTR+ASC -pre ST${ST}
    """
}

process DETECTRECOMBINATION_GUBBINS {

    cpus 2

    input:
    tuple ST, file(tree), file(full_aln) from ch_calculatePhylogeny_detectRecombination.join(ch_makeAlignment_detectRecombination)

    output:
    tuple ST, file("ST${ST}.recombination_predictions.embl") into ch_detectRecombination_mergeHgtRegions

    script:
    """
    run_gubbins.py --starting_tree ${tree} --threads ${task.cpus} --prefix ST${ST} ${full_aln}
    """
}

process FINDHGTSEQUENCE_ALIENHUNTER {

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

process SIMULATEREFFASTQS_WGSIM {

    input:
    tuple ST, file(ref) from ch_simulateFastq
    
    output:
    tuple ST, file("${ref}"), file("${ref_basename}.R1.fq.gz"), file("${ref_basename}.R2.fq.gz")

    script:
    ref_basename = ref.getBaseName()
    """
    wgsim -d 600 -e 0 -N 1500000 -1 250 -2 250 -r 0 -R 0 -X 0 ${ref} ${ref_basename}.R1.fq ${ref_basename}.R2.fq
    gzip *.fq
    """
}

/*
process MAKESNAPPERDBFILES_BWAPICARDSAMTOOLS {
    input:
    output:
    script:
    """

    """
}
*/