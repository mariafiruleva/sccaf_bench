SAMPLES = ["SRS2532209", "SRS2974919", "SRS2756401", "SRS2516420", "SRS3984436", "SRS4558714", "pbmc3k"]

rule all:
    input: expand("counts/plots/{sample_id}/before.png", sample_id=SAMPLES),
         expand("counts/plots/{sample_id}/after.png", sample_id=SAMPLES),
         expand("counts/markers/{sample_id}.tsv", sample_id=SAMPLES),
         expand("scale_data/plots/{sample_id}/before.png", sample_id=SAMPLES),
         expand("scale_data/plots/{sample_id}/after.png", sample_id=SAMPLES),
         expand("scale_data/markers/{sample_id}.tsv", sample_id=SAMPLES)

rule cluster_counts:
    input: "counts/data/{sample_id}.h5ad"
    output: plt_before="counts/plots/{sample_id}/before.png",
          plt_after="counts/plots/{sample_id}/after.png",
          markers="counts/markers/{sample_id}.tsv"
    log: "counts/logs/{sample_id}.txt"
    benchmark: "counts/benchmarks/{sample_id}.txt"
    threads: 1
    shell:
         "python clustering.py --h5 {input} --sample_id {wildcards.sample_id} --plt_before {output.plt_before} --plt_after {output.plt_after} --markers {output.markers} --log {log} 2> {log}"

rule cluster_scale_data:
    input: "scale_data/data/{sample_id}.h5ad"
    output: plt_before="scale_data/plots/{sample_id}/before.png",
          plt_after="scale_data/plots/{sample_id}/after.png",
          markers="scale_data/markers/{sample_id}.tsv"
    log: "scale_data/logs/{sample_id}.txt"
    benchmark: "scale_data/benchmarks/{sample_id}.txt"
    threads: 1
    shell:
         "python clustering.py --h5 {input} --sample_id {wildcards.sample_id} --plt_before {output.plt_before} --plt_after {output.plt_after} --markers {output.markers} --log {log} 2> {log}"
