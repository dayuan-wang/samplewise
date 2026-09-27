# The annotation checkpoint.
#
# Stage A ends with a DRAFT annotation: marker-gene scores per cluster plus,
# when configured, reference-based labels. It is a proposal. A person edits it
# into the locked file named in config["annotation"]["locked"]; stage B does
# not start without that file. Setting config["annotation"]["accept_draft"]
# to true copies the draft into place, which is only appropriate for test
# fixtures whose annotation is already known.

rule draft_annotation:
    input:
        obj=f"{OUT}/04_cluster/clustered.rds",
        markers=f"{OUT}/04_cluster/markers.tsv",
    output:
        draft=f"{OUT}/05_annotation/draft_annotation.tsv",
        report=f"{OUT}/05_annotation/draft_report.html",
    params:
        annotation=config["annotation"],
        species=config["species"],
    log:
        f"{OUT}/logs/05_draft_annotation.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/05_draft_annotation.R"

def locked_annotation(wc):
    if config["annotation"].get("accept_draft", False):
        return f"{OUT}/05_annotation/draft_annotation.tsv"
    return config["annotation"]["locked"]

rule apply_annotation:
    input:
        obj=f"{OUT}/04_cluster/clustered.rds",
        locked=locked_annotation,
    output:
        obj=f"{OUT}/05_annotation/annotated.rds",
        cells=f"{OUT}/05_annotation/cell_metadata.tsv.gz",
    log:
        f"{OUT}/logs/06_apply_annotation.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/06_apply_annotation.R"
