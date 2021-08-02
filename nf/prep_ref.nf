
workflow prep_ref {
    take:
        ref_fasta
        mode
    main:
        ref = ref_fasta ==~ '^(ftp|https?)://.+' ?
            wget(ref_fasta) :
            Channel.from(path(ref))
        if (mode == 'mmi') {
            mmi(ref)
            out = mmi.out.ccs_mmi | mix(mmi.out.subread_mmi)
        } else {
            fai_dict(ref)
            out = fai_dict.out
        }
    emit:
        out
}

process wget {
    label 'S_L'
    publishDir "progress/ref", mode: params.intermediate_pub_mode

    input:
        val url

    output:
        path name

    script:
        name_gz = new File(url).name
        is_gz = name_gz ==~ '.+\\.b?gz$'
        name = name_gz.toString().replaceAll('\\.b?gz$', '')
        if (is_gz)
            """
            wget $url -O $name_gz
            zcat $name_gz > $name
            """
        else
            """
            wget $url -O $name
            """
}

process mmi {
    label 'M'
    publishDir "progress/ref", mode: params.intermediate_pub_mode

    input:
    path ref

    output:
    tuple val('CCS'), path(ccs_mmi), emit: ccs_mmi
    tuple val('SR'), path(subread_mmi), emit: subread_mmi

    script:
    base = ref.toString().replaceAll('.fa(sta)?$', '')
    ccs_mmi = base + '.ccs.mmi'
    subread_mmi = base + '.subread.mmi'
    """
        pbmm2 index $ref $ccs_mmi --preset CCS --num-threads $task.cpus
        pbmm2 index $ref $subread_mmi --preset SUBREAD --num-threads $task.cpus
        """
}

process fai_dict {
    label 'S'
    publishDir "progress/ref", mode: params.intermediate_pub_mode

    input:
        path ref

    output:
        tuple path(ref), path("${ref}.fai"), path(dict)

    script:
        dict = ref.toString().replaceAll('.fa(sta)?$', '.dict')
        """
        samtools faidx $ref
        samtools dict $ref > $dict
        """
}
