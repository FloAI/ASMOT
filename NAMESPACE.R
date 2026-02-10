# Export Classes and Functions
exportClasses(ASMOT)
export(create_asmot)
export(asmot_calculate)
export(asmot_k_scan)
export(asmot_audit)
export(asmot_classifier)
export(asmot_bias_check)
export(asmot_preprocess)
export(asmot_plot_dashboard)
export(asmot_plot_scan)
export(asmot_plot_sparsity)
export(asmot_plot_classifier)

# Export S3 Methods
export(print.asmot_audit_res)
export(print.asmot_k_scan)
export(print.asmot_k_result)

# Imports (Specific to avoid conflicts)
import(methods)
import(stats)
import(ggplot2)
import(reshape2)
import(randomForest)

# Handle specific function imports to avoid namespace clashes
importFrom(transport, wasserstein1d, transport)
importFrom(T4transport, gw, ugw, sinkhorn)
importFrom(phangorn, cophenetic.phylo)
