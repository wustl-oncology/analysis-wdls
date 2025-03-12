version 1.0

task pvacseqAggregatedReportToPreferredTranscriptsList {
  input {
    File pvacseq_aggregated_report
  }

  Int space_needed_gb = 10 + round(size([pvacseq_aggregated_report], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "python:3.11"
    memory: "16GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String out_file = "preferred_transcripts.tsv"
  command <<<
    /usr/bin/python3 -c '
    import csv

    pvacseq_aggregated_report = "~{pvacseq_aggregated_report}"
    preferred_transcript_list = "~{out_file}"

    with open(pvacseq_aggregated_report, "r") as read_fh, open(preferred_transcript_list, "w") as write_fh:
        reader = csv.DictReader(read_fh, delimiter="\t")
        writer = csv.DictWriter(write_fh, delimiter="\t", fieldnames=["CHROM", "POS", "REF", "ALT", "transcript_id"])
        writer.writeheader()
        for line in reader:
            (chrom, start, stop, ref, alt) = line["ID"].split("-")
            if len(ref) == len(alt):
                pos = int(start) + 1
            else:
                pos = start
            out_dict = {
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "transcript_id": line["Best Transcript"]
            }
            writer.writerow(out_dict)
    '
  >>>

  output {
    File preferred_transcripts_tsv = "preferred_transcripts.tsv"
  }
}

workflow wf {
  input {
    File pvacseq_aggregated_report
  }

  call pvacseqAggregatedReportToPreferredTranscriptsList {
    input:
    pvacseq_aggregated_report=pvacseq_aggregated_report
  }
}
