version 1.0

import "../tools/cram_to_bam.wdl" as cb
import "../tools/index_bam.wdl" as i

workflow cramTobamAndIndex{
    input{
        File cram
        File cram_index
        File reference
        File reference_index
        File reference_dict
    }

    call cb.cramToBam {
        input:
        cram=cram,
        cram_index=cram_index,
        reference=reference,
        reference_fai=reference_fai,
        reference_dict=reference_dict,
    }

    call i.index_bam {
        input:
        bam=cramToBam.bam
    }

    output {
        File indexed_bam = index_bam.indexed_bam
        File indexed_bam_bai = index_bam.indexed_bam_bai
        File indexed_bai = index_bam.indexed_bai

    }
}