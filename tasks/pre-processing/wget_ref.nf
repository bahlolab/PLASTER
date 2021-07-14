

process wget_ref {
    label 'S'
    publishDir "output/ref", mode: params.output_pub_mode

    input:
        val(url)

    output:
        path(name)

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