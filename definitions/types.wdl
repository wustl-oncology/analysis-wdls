version 1.0

struct Sequence {
  File? bam
  File? fastq1
  File? fastq2
}
# assume either bam or fastqs defined
struct SequenceData {
  Sequence sequence
  String readgroup
}

struct TrimmingOptions {
  File adapters
  Int min_overlap
}

struct LabelledFile {
  File file
  String label
}
