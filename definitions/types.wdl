version 1.0

struct Sequence {
  File? bam
  File? fastq1
  File? fastq2
}
# assume either bam or fastqs defined
struct SequenceData {
  Sequence sequence
  String? readgroup
}

struct TrimmingOptions {
  File adapters
  Int min_overlap
}

struct LabelledFile {
  File file
  String label
}

# ---- vep_custom_annotation ----
struct Info {
  File file
  Array[File]? secondary_files
  String data_format  # enum, ['bed', 'gff', 'gtf', 'vcf', 'bigwig']
  String name
  Array[String]? vcf_fields
  Boolean? gnomad_filter
  Boolean check_existing
}

struct VepCustomAnnotation {
  Boolean force_report_coordinates
  String method  # enum, ['exact', 'overlap']
  Info annotation
}
