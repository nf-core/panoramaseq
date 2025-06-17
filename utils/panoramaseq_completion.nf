// Subworkflow for pipeline completion specific to nf-core-panoramaseq
workflow PANORAMASEQ_COMPLETION {
    take:
    email           // string: email address for completion notification
    email_on_fail   // string: email address for failure notification
    plaintext_email // boolean: send plain-text email instead of HTML
    outdir          // path: output directory for results
    monochrome_logs // boolean: disable ANSI colour codes in log output
    hook_url        // string: hook URL for notifications
    multiqc_report  // string: path to MultiQC report
    fastqc_status   // map: pass/fail status per sample for FASTQC
    mapping_status  // map: pass/fail status per sample for mapping
    umi_status      // map: pass/fail status per sample for UMI counting

    main:
    def pass_fastqc  = [:]
    def pass_mapping = [:]
    def pass_umi     = [:]

    // Collect status for each step
    fastqc_status
        .map{ id, status -> pass_fastqc[id] = status }
    mapping_status
        .map{ id, status -> pass_mapping[id] = status }
    umi_status
        .map{ id, status -> pass_umi[id] = status }

    // Prepare summary parameters (if needed, adapt to your schema)
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    // Completion email and summary
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }
        // Print a summary for the panoramaseq pipeline
        panoramaseqSummary(monochrome_logs, pass_fastqc, pass_mapping, pass_umi)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Panoramaseq pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

// Function to print a summary for the panoramaseq pipeline
// (You can expand this as needed for your reporting)
def panoramaseqSummary(monochrome_logs=true, pass_fastqc=[:], pass_mapping=[:], pass_umi=[:]) {
    def colors = logColours(monochrome_logs)
    def fail_fastqc  = pass_fastqc.count  { key, value -> value == false }
    def fail_mapping = pass_mapping.count { key, value -> value == false }
    def fail_umi     = pass_umi.count     { key, value -> value == false }
    if (workflow.success) {
        log.info colors.green + "Panoramaseq pipeline completed successfully!" + colors.reset
        log.info "FASTQC failures: ${fail_fastqc}"
        log.info "Mapping failures: ${fail_mapping}"
        log.info "UMI counting failures: ${fail_umi}"
    } else {
        log.error colors.red + "Panoramaseq pipeline failed!" + colors.reset
    }
}
