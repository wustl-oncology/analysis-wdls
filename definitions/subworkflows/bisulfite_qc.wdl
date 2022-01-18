version 1.0

import "../tools/bisulfite_qc_conversion.wdl" as bqc
import "../tools/bisulfite_qc_mapping_summary.wdl" as bqms
import "../tools/bisulfite_qc_cpg_retention_distribution.wdl" as bqcrd
import "../tools/bisulfite_qc_coverage_stats.wdl" as bqcs

workflow bisulfiteQc {
  input {
    File vcf
    File bam
    File reference
    File reference_fai
    File QCannotation
  }

  # TODO(john) really these should all be a single task and this needn't be a subworkflow
  call bqc.bisulfiteQcConversion as bisulfiteConversion {
    input:
    vcf=vcf,
    bam=bam,
    reference=reference,
    reference_fai=reference_fai,
    QCannotation=QCannotation
  }

  call bqms.bisulfiteQcMappingSummary as mappingSummary {
    input:
    vcf=vcf,
    bam=bam,
    reference=reference,
    reference_fai=reference_fai,
    QCannotation=QCannotation
  }

  call bqcrd.bisulfiteQcCpgRetentionDistribution as cpgRetentionDistribution {
    input:
    vcf=vcf,
    bam=bam,
    reference=reference,
    reference_fai=reference_fai,
    QCannotation=QCannotation
  }

  call bqcs.bisulfiteQcCoverageStats as coverageStats {
    input:
    vcf=vcf,
    bam=bam,
    reference=reference,
    reference_fai=reference_fai,
    QCannotation=QCannotation
  }

  output {
    Array[File] qc_files = [
    bisulfiteConversion.base_conversion,
    bisulfiteConversion.read_conversion,
    bisulfiteConversion.cph_retention,
    bisulfiteConversion.cpg_retention,
    mappingSummary.strand_table,
    mappingSummary.mapping_quality,
    cpgRetentionDistribution.cpg_retention_dist,
    coverageStats.bga_bed,
    coverageStats.cov_dist,
    coverageStats.bga_bed_dup,
    coverageStats.dup_report,
    coverageStats.cpg_bed,
    coverageStats.cov_dist_cpg,
    coverageStats.cpg_dist
    ]
  }
}
