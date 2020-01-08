import pandas as pd

import pathlib
import sys

genetic_map = "data/chr22.OMNI.interpolated_genetic_map.gz"
population = "ldetect/examples/example_data/eurinds.txt"
partitions_file = "ldetect/examples/cov_matrix/scripts/chr22_partitions"

prefix = "/mnt/work/endrebak/ldetect_original/1kg"
vcf_prefix = "/mnt/work/endrebak/ldetect/1kg"

if not pathlib.Path(partitions_file).exists():

    rule all:
        input:
            "ldetect/examples/cov_matrix/scripts/chr22_partitions"


    rule partition_chromosomes:
        input:
            # "ldetect/examples/example_data/interpolated_genetic_map.gz"
            genetic_map
        output:
            partitions_file
        shell:
            "python3 ldetect/examples/P00_00_partition_chromosome.py {input[0]} 379 {output}"

else:

    partitions = pd.read_table(partitions_file, sep=" ", header=None, names="Start End".split())
    p = partitions

    starts = p.Start
    ends = p.End
    prefixes = [prefix] * len(starts)


    include: "rules/helpers.smk"


    rule all:
        input:
            expand("{prefix}/ldetect/examples/example_data/cov_matrix/chr22/chr22.{start}.{end}.gz", zip, start=starts, end=ends, prefix=prefixes),
            expand("{prefix}/ldetect2/examples/example_data/cov_matrix/chr22/chr22.{start}.{end}.gz", zip, start=starts, end=ends, prefix=prefixes)


    rule calculate_covariance_matrix:
        input:
            variants = rules.fetch_vcf.output[0],
            genetic_map = genetic_map,
            population = population
        output:
            "{prefix}/ldetect/examples/example_data/cov_matrix/chr22/chr22.{start}.{end}.gz"
        benchmark:
            "{prefix}/ldetect/examples/example_data/cov_matrix/chr22/chr22.{start}.{end}.txt"
        shell:
            "zcat {input.variants} | python3 ldetect/examples/P00_01_calc_covariance.py {input.genetic_map} {input.population} 11418 1e-7 {output[0]}"

    rule calculate_covariance_matrix2:
        input:
            variants = rules.subset_vcf_on_population.output[0],
            genetic_map = genetic_map,
            population = population
        output:
            "{prefix}/ldetect2/examples/example_data/cov_matrix/chr22/chr22.{start}.{end}.gz"
        benchmark:
            "{prefix}/ldetect2/examples/example_data/cov_matrix/chr22/chr22.{start}.{end}.txt"
        shell:
            """zcat {input.variants} | tr "|" "\\t" |
python scripts/calc_covar.py {input.genetic_map} {input.population} 11418 1e-7 | gzip > {output[0]}"""


    rule matrix_to_vector:
        input:
            matrix = "{prefix}/ldetect/examples/example_data/cov_matrix/chr22/chr22.{start}.{end}.gz",
            all_matrixes = expand("{prefix}/ldetect/examples/example_data/cov_matrix/chr22/chr22.{start}.{end}.gz", zip, start=starts, end=ends, prefix=prefixes)
        output:
            "{prefix}/ldetect/examples/example_data/vector/vector-EUR-chr22-{start}-{end}.txt.gz"
        shell:
            "python3 ldetect/examples/P01_matrix_to_vector_pipeline.py --dataset_path={prefix}/ldetect/examples/example_data/cov_matrix/ --name=chr22 --out_fname={output[0]}"


    rule matrix_to_vector2:
        input:
            matrix = "{prefix}/ldetect2/examples/example_data/cov_matrix/chr22/chr22.{start}.{end}.gz",
            all_matrixes = expand("{prefix}/ldetect2/examples/example_data/cov_matrix/chr22/chr22.{start}.{end}.gz", zip, start=starts, end=ends, prefix=prefixes),
            partitions = partitions_file,
            theta2 = rules.calculate_theta2.output[0]
        output:
            "{prefix}/ldetect2/examples/example_data/vectors/chr22/chr22.{start}.{end}.gz"
        benchmark:
            "{prefix}/ldetect2/examples/example_data/vectors/chr22/chr22.{start}.{end}.txt"
        shell:
            "python scripts/matrix_to_vector_much_ram.py {input.partitions} {input.theta2} {input.covariances} # > {output[0]}"


    rule calculate_minima:
        input:
            vector = "{prefix}/ldetect/examples/example_data/vector/vector-EUR-chr22-{start}-{end}.txt.gz",
            all_matrixes = expand("{prefix}/ldetect/examples/example_data/cov_matrix/chr22/chr22.{start}.{end}.gz", zip, prefix=prefixes, start=starts, end=ends)
        output:
            "{prefix}/ldetect/examples/example_data/minima/minima-EUR-chr22-50-{start}-{end}.pickle"
        shell:
            "python3 ldetect/examples/P02_minima_pipeline.py --input_fname={input.vector} --chr_name=chr22 --dataset_path={prefix}/ldetect/examples/example_data/cov_matrix/ --n_snps_bw_bpoints=50 --out_fname={output[0]}"


    rule calculate_minima2:
        input:
            vector = "{prefix}/ldetect2/examples/example_data/vectors/chr22/chr22.{start}.{end}.gz",
            all_matrixes = expand("{prefix}/ldetect2/examples/example_data/cov_matrix/chr22/chr22.{start}.{end}.gz", zip, prefix=prefixes, start=starts, end=ends)
        output:
            "{prefix}/ldetect2/examples/example_data/minima/minima-EUR-chr22-50-{start}-{end}.pickle"
        shell:
            "python scripts/find_minima.py"



