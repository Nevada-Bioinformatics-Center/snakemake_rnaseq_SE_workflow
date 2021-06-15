def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

def get_fastq1(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna().item()

def get_fastq2(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq2"]].dropna().item()

rule uncompress_bz2:
    input:
        get_fastq1,
    output:
        temp("uncompressed/{sample}-{unit}.1.fastq")
    threads: 1
    resources: time_min=120, mem_mb=2000, cpus=1
    shell:
        "bzcat {input} > {output}"
 

rule trimmomatic_se:
    input:
        temp("uncompressed/{sample}-{unit}.1.fastq"),
    output:
        "trimmed/{sample}-{unit}.1.fastq.gz",
    log:
        "logs/trimmomatic_se/{sample}-{unit}.log"
    params:
        # list of trimmers (see manual)
        #trimmer = ["ILLUMINACLIP:/data/gpfs/assoc/inbre/hansvg/common_references/CustomBlacklist.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36"],
	#f"{config['directory']}/{{sample}}_L{{lanenum}}_R1_001.fastq.gz"
        trimmer = [f"ILLUMINACLIP:{config['ref']['adapter']}:2:30:10 SLIDINGWINDOW:4:20 LEADING:30 MINLEN:36"],
        # optional parameters
        extra="",
        compression_level="-9"
    threads: 24
    resources: time_min=820, mem_mb=20000, cpus=24
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    wrapper:
        "v0.69.0/bio/trimmomatic/se"

#rule fastp_pe:
#    input:
#        sample=get_fastq
#    output:
#        trimmed=["trimmed/{sample}-{unit}.1.fastq.gz", "trimmed/{sample}-{unit}.2.fastq.gz"],
#        html="report/pe/{sample}-{unit}.html",
#        json="report/pe/{sample}-{unit}.json"
#    log:
#        "logs/fastp/pe/{sample}-{unit}.log"
#    params:
#        adapters="--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
#        extra=""
#    threads: 4
#    wrapper:
#        "0.72.0/bio/fastp"
