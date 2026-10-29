param([switch]$Extended)
$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot
function Invoke-Checked {
    param([string]$Program, [string[]]$Parameters)
    & $Program @Parameters
    if ($LASTEXITCODE -ne 0) { throw "Failed: $Program $Parameters" }
}
Invoke-Checked 'Rscript' @('1.prepare_transcriptomic_data.R')
Invoke-Checked 'python' @('2.lasso_nested_cross_validation.py')
Invoke-Checked 'python' @('3.mlp_training_and_deepshap.py')
Invoke-Checked 'python' @('4.candidate_gene_integration.py')
Invoke-Checked 'Rscript' @('5.wgcna_network_analysis.R')
Invoke-Checked 'Rscript' @('6.differential_expression.R')
Invoke-Checked 'Rscript' @('7.functional_enrichment.R')
if ($Extended) {
    Invoke-Checked 'python' @('8.robustness_and_external_associations.py')
    Invoke-Checked 'python' @('9.pca_kernelshap_sensitivity.py')
    foreach ($Target in @('SF2', 'AUC_recomputed')) {
        foreach ($Count in @(3000, 5000, 10000)) {
            if ($Target -eq 'SF2' -and $Count -eq 5000) { continue }
            Invoke-Checked 'python' @('3.mlp_training_and_deepshap.py', '--target', $Target, '--genes', "$Count")
        }
    }
}
