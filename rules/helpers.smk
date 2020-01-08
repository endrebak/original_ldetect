# f = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"

rule download_vcf:
    output:
        f"{vcf_prefix}/1kg/chr22.vcf.gz"
    shell:
        "axel http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -q -o {output[0]}"


rule download_vcf_index:
    output:
        f"{vcf_prefix}/1kg/chr22.vcf.gz.tbi"
    shell:
        "axel http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz.tbi -q -o {output[0]}"


rule subset_on_population:
    input:
        variants = rules.download_vcf.output[0],
        samples = population,
        index = rules.download_vcf_index.output[0]
    output:
        f"{prefix}/1kg/chr22.vcf.gz"
    shell:
        "bcftools view --threads 48 --force-samples -O z -S {input.samples} {input.variants} > {output[0]}"


rule index_population_vcf:
    input:
        variants = rules.subset_on_population.output[0]
    output:
        f"{prefix}/1kg/chr22.vcf.gz.tbi"
    shell:
        "tabix -f {input[0]}"


rule fetch_vcf:
    output:
        f"{vcf_prefix}/1kg/chr22_{{start}}_{{end}}.vcf.gz"
    shell:
        '''tabix -h http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz 22:{wildcards.start}-{wildcards.end} | grep -v '##' | gzip > {output[0]}'''


rule subset_vcf_on_population:
    input:
        variants = rules.fetch_vcf.output[0],
        samples = population
    output:
        f"{vcf_prefix}/1kg/population/chr22_{{start}}_{{end}}.vcf.gz"
    run:
        p = pd.read_table(input.samples, header=None, squeeze=True).tolist()

        cols = ["POS", "ID"] + p

        df = pd.read_table(input.variants, usecols=cols, index_col=[0, 1])

        df = df.apply(lambda s: s.str.split(":", expand=True)[0])

        df = df.reset_index()

        df.to_csv(output[0], sep="\t", index=False, header=None)



rule calculate_theta2:
    input:
        population
    output:
        "{prefix}/ldetect2/theta2/{chromosome}.txt"
    run:
        inds = pd.read_table(input[0], header=None, squeeze=True).to_list()

        nind_int = len(inds)
        s = 0

        for _i in range(1, 2*nind_int):
            s = s+ 1.0/float(_i)

            s = 1/s

            theta = s/(2.0*float(nind_int)+s)
            thetas2 = (theta/2.0)*(1-theta/2.0)

            with open(output[0], "w+") as o:
                o.write(str(thetas2) + "\n")
