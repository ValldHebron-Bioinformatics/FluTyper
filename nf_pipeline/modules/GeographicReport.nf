#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process GeographicReport {
    errorStrategy 'ignore'

    input:
    file(genotyping_file)
    file(metadata_file)

    output:
    file("GeographicReport.html"), emit: report

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import plotly
    import os
    
    """
}
