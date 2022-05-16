
process get_pharmvar_vcf {
    label 'S4'
    publishDir "progress/get_pharmvar_vcf",  mode: params.intermediate_pub_mode
    tag {"${pharmvar.ver}:${pharmvar.ref}:${pharmvar.gene}"}

    input:
        tuple val(am), val(pharmvar)
        tuple path(ref), path(fai), path(dict)


    output:
        tuple val(am), path(vcf), path("${vcf}.tbi")

    script:
        url = "https://www.pharmvar.org/get-download-file?name=ALL&refSeq=ALL&fileType=zip&version=${pharmvar.ver}"
        vcf = "${am}.${pharmvar.ver}_${pharmvar.ref}_${pharmvar.gene}.vcf.gz"
        """
        curl "$url" -o bundle.zip
        unzip -l bundle.zip |
            egrep -o pharmvar-.+/${pharmvar.gene}/${pharmvar.ref}/.+vcf\$ |
            xargs -n1 -i bash -c '
                NAME=\$(basename {} .vcf)
                BCF=\$NAME.bcf
                ANN=\$NAME.ann.bcf
                unzip -p bundle.zip {} |
                    sed "s:INFO\t\$:INFO:" |
                    bcftools view -Ob -o \$BCF
                bcftools index \$BCF
                bcftools annotate \$BCF -a \$BCF -m \$NAME -Ob -o \$ANN
                bcftools index \$ANN'
        bcftools merge -m none *.ann.bcf -Ou |
            bcftools norm -f $ref -Oz -o $vcf 
        bcftools index -t $vcf
        rm *.bcf *.csi bundle.zip
        """
}