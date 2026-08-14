version 1.0

# Zip a directory so it can cross a WDL task boundary.
#
# WHY THIS EXISTS
# tools/pvacfuse.wdl takes `File input_fusions_zip` and unzips it into
# `agfusion_dir` before handing that directory to `pvacfuse run` -- the zip is a
# transport wrapper, not something pvacfuse ever sees, because WDL 1.0 has no
# portable Directory type.
#
# To let pvacfuse_longread.wdl accept an AGFusion directory the user already has
# on disk, the cheapest correct move is to re-zip it and feed the existing,
# unmodified pvacfuse task. That keeps tools/pvacfuse.wdl untouched (and still
# valid WDL 1.0, since Directory stays out of it), keeps a single pvacfuse call
# so no downstream output types change, and costs a few seconds -- AGFusion
# output is a tree of small text files.
#
# NOTE: uses the Directory type, which is not in the WDL 1.0 spec. This file is
# only reachable from the long-read branch, which already depends on Directory
# via tools/ctat_lr_fusion.wdl.
#
# The default container is the same image tools/agfusion.wdl uses, chosen
# because it demonstrably ships `zip` (agfusion.wdl calls it) and will already
# be present on the node in the normal ctat -> agfusion -> pvacfuse path.
task zipDirectory {
  input {
    Directory dir
    String output_basename = "directory"
    String docker_image = "mgibio/agfusion:1.3.11-ensembl-105"
  }

  runtime {
    preemptible: 1
    maxRetries: 2
    docker: docker_image
    memory: "8GB"
    cpu: 2
    disks: "local-disk 50 HDD"
  }

  command <<<
    set -eou pipefail

    if [ ! -d "~{dir}" ]; then
      echo "ERROR: ~{dir} is not a directory" >&2
      exit 1
    fi

    echo "top-level entries in ~{dir}:"
    ls -1 "~{dir}" | head -20
    echo "(total: $(ls -1 "~{dir}" | wc -l))"

    out="$PWD/~{output_basename}.zip"
    cd "~{dir}"
    # Zip the CONTENTS, not the directory itself. This matches exactly what
    # tools/agfusion.wdl produces (`cd $output_dir && zip -r ../x.zip .`), so
    # pvacfuse sees an identical layout after it unzips.
    zip -r "$out" . 1> /dev/null
  >>>

  output {
    File zipped = output_basename + ".zip"
  }
}

workflow wf {
  input {
    Directory dir
    String? output_basename
  }
  call zipDirectory {
    input: dir=dir, output_basename=select_first([output_basename, "directory"])
  }
}
