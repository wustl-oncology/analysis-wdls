version 1.0

task starFusionDetect {
  input {
    File star_fusion_genome_dir_zip
    File junction_file
    String fusion_output_dir = "STAR-Fusion_outdir"
    String star_path = "/usr/local/bin/STAR"
  }

  Int cores = 10
  Int space_needed_gb = 10 + round(5*size([junction_file, star_fusion_genome_dir_zip], "GB"))
  runtime {
    memory: "64GB"
    cpu: cores
    docker: "trinityctat/starfusion:1.8.0"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

    # https://github.com/STAR-Fusion/STAR-Fusion/issues/175#issuecomment-567913451
  String genome_lib_dir = "`pwd`/" + basename(star_fusion_genome_dir_zip, ".zip")
  command <<<
    mkdir ~{genome_lib_dir} && unzip -qq ~{star_fusion_genome_dir_zip} -d ~{genome_lib_dir}
    /usr/local/src/STAR-Fusion/STAR-Fusion --CPU ~{cores} \
        --genome_lib_dir ~{genome_lib_dir} \
        -J ~{junction_file} --output_dir ~{fusion_output_dir} --STAR_PATH ~{star_path}
  >>>

  output {
    File fusion_predictions = fusion_output_dir + "/star-fusion.fusion_predictions.tsv"
    File fusion_abridged = fusion_output_dir + "/star-fusion.fusion_predictions.abridged.tsv"
  }
}

workflow wf { call starFusionDetect { input: } }
