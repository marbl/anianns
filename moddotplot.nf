params.sequence = "$projectDir/*.fa"
params.pythonscript1 = "$projectDir/parsemoddotplot.py"
params.pythonscript2 = "$projectDir/parsediagonals.py"
params.window_size = 10000
params.identity = 90
params.parse_output = "output.bed"

process runModDotPlot {
    input:
    path sequence
    val window_size
    val identity

    output:
    path 'tmp_out/*.bed'

    script:
    """
    moddotplot static -f $sequence -o tmp_out --no-plot -w $window_size --identity $identity
    """
}

process splitDiagonals {
    input:
    path bedfile

    output:
    path "diagonal.bed"

    script:
    """
    awk '\$2 == \$3' $bedfile > diagonal.bed
    """
}

process parseModDotPlot {
    input:
    path bedfile
    path pythonscript

    output:
    path "grouped_identity.bed"

    script:
    """
    python $pythonscript $bedfile grouped_identity.bed
    """
}

process parseDiagonals {
    input:
    path bedfile
    path pythonscript

    output:
    path "satellites.bed"

    script:
    """
    python $pythonscript $bedfile tmpfile_delete satellites.bed
    """
}

workflow {
    runModDotPlot(params.sequence,params.window_size,params.identity)
    parseModDotPlot(runModDotPlot.out, params.pythonscript1)
    splitDiagonals(parseModDotPlot.out)
    parseDiagonals(splitDiagonals.out, params.pythonscript2)
}