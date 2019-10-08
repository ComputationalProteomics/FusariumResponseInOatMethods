from collections import OrderedDict
import sys

configfile: "config.yaml"

settings = config['runtypes'][config['runtype']]

for setting_name in settings:
    config['settings'][setting_name] = settings[setting_name]

config['paths']['base_out_analysis'] = config['paths']['base_out_analysis'] + '_' + config['runtype']

if config['settings']['input_molecule'] == 'pep':
    molecules = ['pep', 'prot']
elif config['settings']['input_molecule'] == 'prot':
    molecules = ['prot']
else:
    sys.exit("Unknown input molecule: {}, only pep or prot permitted".format(config['settings']['input_molecule']))

if config['settings']['do_norm']:
    norms = config['settings']['select_norm'].split(' ')
else:
    norms = ['nonorm']

if config['settings']['do_batch_split']:
    splits = ['both', 'hours4', 'day1', 'day2', 'day4']
    # splits = ['both', 'batch1', 'batch2']
else:
    splits = ['all']
oats = ['Argamak', 'Belinda']

setupd = config['procdirs']['setup']
normd = config['procdirs']['norm']
rollupd = config['procdirs']['rollup']
statd = config['procdirs']['stat']
annotd = config['procdirs']['annot']
finald = config['procdirs']['final']

blast_headers = OrderedDict()
blast_headers['qseqid'] = 'QueryID'
blast_headers['sseqid'] = 'ProteinID'
blast_headers['pident'] = 'PercentIdentity'
blast_headers['length'] = 'HitLength'
blast_headers['mismatch'] = 'Mismatches'
blast_headers['gapopen'] = 'GapOpens'
blast_headers['qstart'] = 'QueryStart'
blast_headers['qend'] = 'QueryEnd'
blast_headers['sstart'] = 'SubjStart'
blast_headers['send'] = 'SubjEnd'
blast_headers['evalue'] = 'EValue'
blast_headers['bitscore'] = 'Bitscore'
blast_headers['stitle'] = 'Description'

rule all:
    input:
        config['paths']['base_out_setup'] + "/Snakefile",
        config['paths']['base_out_setup'] + "/config.yaml",
        config['paths']['base_out_setup'] + "/comb_oats.blast6e",
        expand(config['paths']['base_out_analysis'] + annotd + '/data_{split}_{norm}_{molecule}_stats_annot.tsv', split=splits, norm=norms, molecule=molecules),
        config['paths']['base_out_analysis'] + finald + '/combined_flat_ses.rds'
        # expand(config['paths']['base_out_analysis'] + annotd + '/data_{split}_{norm}_prot_enrichment.tsv', split=splits, norm=norms)

rule copy_run_files:
    input:
        snakefile="Snakefile",
        configfile="config.yaml"
    output:
        snakefile_proc=config['paths']['base_out_setup'] + "/Snakefile",
        configfile_proc=config['paths']['base_out_setup'] + "/config.yaml",
        snakefile_analysis=config['paths']['base_out_analysis'] + setupd + "/Snakefile",
        configfile_analysis=config['paths']['base_out_analysis'] + setupd + "/config.yaml"
    run:
        shell("""
        cp {input.snakefile} {output.snakefile_proc}
        cp {input.configfile} {output.configfile_proc}
        cp {input.snakefile} {output.snakefile_analysis}
        cp {input.configfile} {output.configfile_analysis}
        """)

rule setup_blast_dbs:
    input:
        arab_fasta = {config['annot']['arabidopsis_proteins']},
        fusarium_fasta = {config['annot']['fusarium_proteins']},
        wheat_fasta = {config['annot']['wheat_proteins']}
    output:
        search_fasta=config['paths']['base_out_setup'] + "/search_protein.fasta",
        search_db=config['paths']['base_out_setup'] + "/search_protein.fasta.phr"
    run:
        if config['settings']['do_annot_subset'] is True:
            print("Doing subset run using " + str(config['settings']['annot_subset_count']) + " database features")
            arab_subset = config['paths']['base_out_setup'] + '/arab_subset.fa'
            fusarium_subset = config['paths']['base_out_setup'] + '/fusarium_subset.fa'
            wheat_subset = config['paths']['base_out_setup'] + '/wheat_subset.fa'
            shell("""
                head -n {config[settings][annot_subset_count]} {input.arab_fasta} > {arab_subset}
                head -n {config[settings][annot_subset_count]} {input.fusarium_fasta} > {fusarium_subset}
                cat {arab_subset} {fusarium_subset} > {output.search_fasta}
                makeblastdb -in {output.search_fasta} -parse_seqids -dbtype prot                
            """)
        else:
            print("Building a full BLAST database")
            shell("""
                cat {input.arab_fasta} {input.fusarium_fasta} > {output.search_fasta}
                makeblastdb -in {output.search_fasta} -parse_seqids -dbtype prot
                """)

rule blasts_runs:
    input:
        query = config['paths']['raw_dir'] + "/180918_databases/Oat_{oat}_w_rev_20150630.fasta",
        blast_db = config['paths']['base_out_setup'] + "/search_protein.fasta",
        blast_db_index = config['paths']['base_out_setup'] + "/search_protein.fasta.phr"
    output:
        blast_results = config['paths']['base_out_setup'] + "/{oat}.blast6e"

    # max_hsps - Maximum high-scoring pairs per sequence
    run:
        out_fields = ' '.join(list(blast_headers.keys()))
        if config['settings']['do_query_subset'] is False:
            print("Performing full BLAST run, MIGHT BE SLOW")
            shell("""
                echo "Number queries"
                grep -c "^>" {input.query}
                
                echo "Database size"
                grep -c "^>" {input.blast_db}
            
                time blastp \
                    -query {input.query} \
                    -db {input.blast_db} \
                    -outfmt "6 {out_fields}" \
                    -max_target_seqs 1 \
                    -max_hsps 1 \
                    -evalue {config[settings][evalue_thres]} \
                    -num_threads {config[settings][blast_threads]} > {output.blast_results}
                """)
        else:
            print("Performing limited BLAST run")
            shell("""
                echo "Subsetting input queries, final count"
                
                seqtk sample {input.query} {config[settings][query_subset_count]} > {config[paths][base_out_setup]}/query_subset_{wildcards.oat}.fasta
                grep -c "^>" {config[paths][base_out_setup]}/query_subset.fasta
                
                echo "Database size"
                grep -c "^>" {input.blast_db}
            
                time blastp \
                    -query query_subset_{wildcards.oat}.fasta \
                    -db {input.blast_db} \
                    -outfmt "6 {out_fields}" \
                    -max_target_seqs 1 \
                    -max_hsps 1 \
                    -evalue {config[settings][evalue_thres]} \
                    -num_threads {config[settings][blast_threads]} > {output.blast_results}
                """)

rule combine_blast_results:
    input:
        all_blasts = expand(config['paths']['base_out_setup'] + "/{oat}.blast6e", oat=oats)
    output:
        combined_blasts = config['paths']['base_out_setup'] + "/comb_oats.blast6e"
    run:
        blast_header_tabs = '\t'.join(list(blast_headers.values()))
        shell("""
            echo -e "{blast_header_tabs}" > {output.combined_blasts}
            cat {input.all_blasts} >> {output.combined_blasts}
            """)

# Copy raw into analysis folder for reference
rule copy_raw:
    input:
        data_matrix=config['settings']['data_path']
    output:
        data_matrix=config['paths']['base_out_analysis'] + setupd + '/data_raw.tsv'
    run:
        shell("""
            cp {input.data_matrix} {output.data_matrix}
        """)

# Clear out entries matching certain patterns
rule cleanup_data:
    input:
        data_matrix=config['paths']['base_out_analysis'] + setupd + '/data_raw.tsv'
    output:
        data_matrix=config['paths']['base_out_analysis'] + setupd + '/data_clean.tsv'
    run:
        shell("""
            Rscript {config[paths][rworkflowbase]}/clean_data_on_column.R \
                --in_fp {input.data_matrix} \
                --out_fp {output.data_matrix} \
                --target_column {config[settings][prot_col]} \
                --target_patterns {config[settings][clean_patterns]} \
                --delim "{config[settings][rdf_name_col_splitter]}" \
                --screen_any_match TRUE
        """)

if config['settings']['do_batch_split']:
    rule split_data_into_batches:
        input:
            data_matrix=config['paths']['base_out_analysis'] + setupd + '/data_clean.tsv',
            design_matrix=config['settings']['design_path']
        output:
            split1_data=config['paths']['base_out_analysis'] + setupd + '/data_clean_batch1.tsv',
            split1_design=config['paths']['base_out_analysis'] + setupd + '/design_batch1.tsv',
            split2_data=config['paths']['base_out_analysis'] + setupd + '/data_clean_batch2.tsv',
            split2_design=config['paths']['base_out_analysis'] + setupd + '/design_batch2.tsv',
            both_data=config['paths']['base_out_analysis'] + setupd + '/data_clean_both.tsv',
            both_design=config['paths']['base_out_analysis'] + setupd + '/design_both.tsv'
        run:
            shell("""
                Rscript {config[paths][projectscripts]}/split_dataset.R \
                    --in_data_fp {input.data_matrix} \
                    --in_design_fp {input.design_matrix} \
                    --sample_col {config[settings][sample_col]} \
                    --split_col {config[settings][batch_col]} \
                    --split1_data_fp {output.split1_data} \
                    --split1_design_fp {output.split1_design} \
                    --split2_data_fp {output.split2_data} \
                    --split2_design_fp {output.split2_design} \
                    --load_tools_script {config[scripts][load_tools]}
                
                echo "Dataset split, now copying original files.."
                cp {input.data_matrix} {output.both_data}
                cp {input.design_matrix} {output.both_design}
                echo "Done!"
            """)
else:
    rule no_split:
        input:
             data_matrix=config['paths']['base_out_analysis'] + setupd + '/data_clean.tsv',
             design_matrix=config['settings']['design_path']
        output:
              out_data=config['paths']['base_out_analysis'] + setupd + '/data_clean_all.tsv',
              out_design=config['paths']['base_out_analysis'] + setupd + '/design_all.tsv'
        run:
            shell("""
                echo "No split specified, copying original files.."
                cp {input.data_matrix} {output.out_data}
                cp {input.design_matrix} {output.out_design}
                echo "Done!"
            """)

rule reparse_ids:
    input: data_matrix=config['paths']['base_out_analysis'] + setupd + '/data_clean_{split}.tsv'
    output: data_matrix=config['paths']['base_out_analysis'] + setupd + '/data_clean_{split}_idsok.tsv'
    run:
        if config['settings']['reparse_ids']:
            print("Reparsing IDs...")
            shell("""
                Rscript {config[paths][rworkflowbase]}/sort_column_fields.R \
                    --in_data_fp {input.data_matrix} \
                    --out_data_fp {output.data_matrix} \
                    --target_column {config[settings][prot_col]} \
                    --delim "{config[settings][rdf_name_col_splitter]}"
            """)
        else:
            print("Skipping reparsing, simply copying")
            shell("""
                cp {input.data_matrix} {output.data_matrix}
            """)

if config['settings']['do_norm']:
    rule normalization:
        input:
            design_matrix=config['paths']['base_out_analysis'] + setupd + '/design_{split}.tsv',
            data_matrix=config['paths']['base_out_analysis'] + setupd + '/data_clean_{split}_idsok.tsv',
        output:
            out_matrix=config['paths']['base_out_analysis'] + normd + '/data_{split}_{norm}_%s.tsv' % config['settings']['input_molecule']
        run:
            shell("""
                Rscript {config[paths][rworkflowbase]}/quick_normalization.R \
                    --ddf_fp {input.design_matrix} \
                    --rdf_fp {input.data_matrix} \
                    --sample_col {config[settings][sample_col]} \
                    --out_fp {output.out_matrix} \
                    --do_log2 {config[settings][do_log2]} \
                    --do_zscale {config[settings][do_zscale]} \
                    --normalization {wildcards.norm} \
                    --na_val {config[settings][na_val]}
            """)
else:
    print("No normalization specified")
    rule skip_normalization:
        input:
            data_matrix=config['paths']['base_out_analysis'] + setupd + '/data_clean_{split}_idsok.tsv'
        output:
            data_matrix=config['paths']['base_out_analysis'] + normd + '/data_{split}_{norm}_%s.tsv' % config['settings']['input_molecule']
        run:
            shell("""
                cp {input.data_matrix} {output.data_matrix}
            """)

# Skip this if already rolled up?
if config['settings']['input_molecule'] == 'pep':
    rule protein_rollup:
        input:
            peptide_fp=config['paths']['base_out_analysis'] + normd + '/data_{split}_{norm}_pep.tsv',
            design_fp=config['paths']['base_out_analysis'] + setupd + '/design_{split}.tsv'
        output:
            protein_fp=config['paths']['base_out_analysis'] + normd + '/data_{split}_{norm}_prot.tsv'
        run:
            shell("""
                Rscript {config[paths][rworkflowbase]}/protein_rollup.R \
                    --ddf_fp {input.design_fp} \
                    --rdf_fp {input.peptide_fp} \
                    --out_fp {output.protein_fp} \
                    --sample_col {config[settings][sample_col]} \
                    --peptide_col {config[settings][pep_col]} \
                    --protein_col {config[settings][prot_col]} \
                    --protein_tools_path {config[scripts][protein_tools]} \
                    --out_protein_name {config[settings][prot_col]}
            """)
else:
    print("Input not peptides, skipping protein rollup")

contrasts = config['contrasts'][config['settings']['contrast_name']]
rule calculate_statistics:
    input:
         design_fp=config['paths']['base_out_analysis'] + setupd + '/design_{split}.tsv',
         data_fp=config['paths']['base_out_analysis'] + normd + '/data_{split}_{norm}_{molecule}.tsv'
    output:
         out_fp=config['paths']['base_out_analysis'] + statd + '/data_{split}_{norm}_{molecule}_stats.tsv'
    run:
        shell("""
            Rscript {config[paths][rworkflowbase]}/limma_contrasts.R \
                --ddf_fp          {input.design_fp} \
                --rdf_fp          {input.data_fp} \
                --out_fp          {output.out_fp} \
                --sample_col      {config[settings][sample_col]} \
                --model           {config[stats][model]} \
                --contrasts       {contrasts[contrasts]} \
                --contrast_names  {contrasts[contrast_names]} \
                --cond_col_name   {config[stats][stat_col]} \
                --omit_absent_contrasts TRUE \
                --do_presence_absence TRUE
        """)

rule annotate_proteins:
    input:
        data_fp=config['paths']['base_out_analysis'] + statd + '/data_{split}_{norm}_{molecule}_stats.tsv',
        annot_fp=config['paths']['base_out_setup'] + '/comb_oats.blast6e'
    output:
        data_fp=config['paths']['base_out_analysis'] + annotd + '/data_{split}_{norm}_{molecule}_stats_annot.tsv'
    run:
        shell("""
            Rscript {config[paths][rworkflowbase]}/annotate_table.R \
                --rdf                    {input.data_fp} \
                --db                     {input.annot_fp} \
                --rdf_name_col           {config[settings][prot_col]} \
                --db_name_col            "QueryID" \
                --out                    {output.data_fp} \
                --rdf_name_col_splitter  "{config[settings][rdf_name_col_splitter]}" \
                --clean_up_descr_type    "TAIR" \
                --clean_up_descr_colname "Description"
        """)

rule assembly_presence:
    input:
        data_fp=config['paths']['base_out_analysis'] + annotd + '/data_{split}_{norm}_{molecule}_stats_annot.tsv',
        design_fp=config['paths']['base_out_analysis'] + setupd + '/design_{split}.tsv'
    output:
        data_fp=config['paths']['base_out_analysis'] + annotd + '/data_{split}_{norm}_{molecule}_stats_annot_annotpres.tsv'
    run:
        shell("""
            Rscript {config[paths][projectscripts]}/assembly_presence.R \
                --in_data_fp         {input.data_fp} \
                --in_design_fp       {input.design_fp} \
                --out_data_fp        {output.data_fp} \
                --sample_col         {config[settings][sample_col]} \
                --annot_col          {config[settings][prot_col]} \
                --load_tools_script  {config[scripts][load_tools]}
        """)

# # http://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/gage.pdf
# # Is this one used?
# rule gene_ontology:
#     input:
#         data_fp=config['paths']['base_out_analysis'] + annotd + '/data_{split}_{norm}_prot_stats_annot_annotpres.tsv',
#         design_fp=config['paths']['base_out_analysis'] + setupd + '/design_{split}.tsv'
#     output:
#         data_fp=config['paths']['base_out_analysis'] + annotd + '/data_{split}_{norm}_prot_enrichment.tsv'
#     run:
#         shell("""
#             Rscript {config[paths][rworkflowbase]}/enrichment.R \
#                 --rdf_fp            {input.data_fp} \
#                 --ddf_fp            {input.design_fp} \
#                 --sample_col        {config[settings][sample_col]} \
#                 --cond_col          {config[settings][group_col]} \
#                 --annot_col         "ProteinID" \
#                 --enrichment_out_fp {output.data_fp} \
#                 --contrasts         {config[stats][contrasts]} \
#                 --contrast_names    {config[stats][contrast_names]} \
#                 --species           "Arabidopsis" \
#                 --set_type          "GO" \
#                 --trim_reg          '\\.\\d$' \
#                 --show_pars         TRUE
#         """)

rule build_ses:
    input:
        design=config['paths']['base_out_analysis'] + setupd + '/design_{split}.tsv',
        datasets=config['paths']['base_out_analysis'] + annotd + '/data_{split}_{norm}_{molecule}_stats_annot_annotpres.tsv'
    output:
        ses_obj=config['paths']['base_out_analysis'] + finald + '/{split}_{norm}_{molecule}_ses_obj.rds'
    run:
        shell("""
            Rscript {config[paths][rworkflowbase]}/combine_to_ses.R \
                --ddf_fp     {input.design} \
                --rdf_fps    {input.datasets} \
                --out_fp     {output.ses_obj} \
                --sample_col {config[settings][sample_col]}
        """)

rule combine_ses:
    input:
        ses_objs=expand(config['paths']['base_out_analysis'] + finald + '/{split}_{norm}_{molecule}_ses_obj.rds', split=splits, norm=norms, molecule=molecules)
    output:
        combined_flat=config['paths']['base_out_analysis'] + finald + '/combined_flat_ses.rds',
        combined_nested=config['paths']['base_out_analysis'] + finald + '/combined_nested_ses.rds'
    run:
        shell("""
            Rscript {config[paths][rworkflowbase]}/concatenate_ses.R \
                --input_ses  {input.ses_objs} \
                --out_nested {output.combined_nested} \
                --out_flat   {output.combined_flat}
        """)





