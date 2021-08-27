# Differences from CWL

There is no Directory type in CWL. The one instance of this,
`vep_cache_dir`, has been changed to `vep_cache_dir_zip` and expects a
zipped file of the dir, where the zip contains the contents of the dir
not the dir itself.

As far as I can tell there's no equivalent of
"InitialWorkDirRequirement" from CWL in WDL. Where we use this to
create a script to execute, that scripts contents have been moved to
the `command` block of their task. No difference to the caller.

Input files must prefix each argument with the name of the workflow
they're going to run, because a WDL file can contain multiple
workflows or pass inputs over a layer if they aren't propagated
through in the definition. e.g. to call workflow `somaticExome` with
input `foo`, yaml key must be `somaticExome.foo`

WDL does not allow relative imports on the root file, so you'll need a
local copy that strips the leading `../` from each import. In CWL a
similar issue happened but had a workaround of zipping the deps in a
tricky way to allow it. WDL bans it outright.


# TODO:

Various tools need space allocation for their outputs

## tools/vep.wdl

Actually use custom annotations that are passed in


## tools/cnvkit_batch.wdl

Docker version update from 0.9.5 to 0.9.8 causes a regression and we
need a solution because it enables CRAM support. Specifically
permissions to write are denied for commands like `mv`, `touch`, and
generated index files within the script. It _does_ work locally, so
troubleshooting is probably required before filing a report to CNVkit
owners.

## Runtime on cloud is slower than cluster, why?

Timing diagrams for the same run between cloud and cluster show that
the cloud runs about 50% longer than the cluster. Localization
accounts for a portion of this, but task runtimes are also slower on
the cloud.
