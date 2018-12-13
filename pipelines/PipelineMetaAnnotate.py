'''
classes and utility functions for pipeline_metaannotate.py

'''

import os
import PipelineMetaAssemblyKit

'''
function to build call to Prodigal
'''

def runProdigal(infile,outfile,params):
    pcall = ["prodigal -i {} -o {}".format(infile,outfile.rstrip("_peptides")+"_positions")]
    pcall.append("-a {}".format(outfile))
    if params["Prodigal_c"] != 'false':
        pcall.append("-c")
    if params["Prodigal_d"] != 'false':
        pcall.append("-d")
    pcall.append("-f {}".format(params["Prodigal_f"]))
    pcall.append("-g {}".format(params["Prodigal_g"]))
    if params["Prodigal_m"] != 'false':
        pcall.append("-m")
    if params["Prodigal_n"] != 'false':
        pcall.append("-n")
    pcall.append("-p {}".format(params["Prodigal_p"]))
    if params["Prodigal_q"] != 'false':
        pcall.append("-q")
    if params["Prodigal_s"] != 'false':
        pcall.append("-s {}".format(outfile.rstrip("_peptides")+"."+params["Prodigal_s"]))
    if params["Prodigal_t"] != 'false':
        pcall.append("-t {}".format(outfile.rstrip("_peptides")+"."+params["Prodigal_t"]))
    return(" ".join(pcall))

'''
function to build call to eggnog-annotation
'''
def runEggmapAnnot(infile,outfile,params):
    return("python2 {} --annotate_hits_table {} --no_file_comments -o {} --cpu {} --data_dir {}".format(params["Eggnogmapper_eggpath"],infile,outfile,params["Eggnogmapper_threads_annot"],params["Eggnogmapper_eggdata"]))


'''
function to build call to eggnog-mapper seed alignment
'''

def runEggmapSeed(infile,outfile,params):
    #eggnog-mapper requires python2
    pcall = ["python2 {}".format(params["Eggnogmapper_eggpath"])]
    #input output
    pcall.append("-i {} -o {}".format(infile, outfile))
    #set database dir and alignment method
    pcall.append("--data_dir {}".format(params["Eggnogmapper_eggdata"]))
    pcall.append("-m {}".format(params["Eggnogmapper_method"]))
    #annotation settings
    pcall.append("--tax_scope {}".format(params["Eggnogmapper_tax_scope"]))
    pcall.append("--target_orthologs {}".format(params["Eggnogmapper_target_orthologs"]))
    pcall.append("--go_evidence {}".format(params["Eggnogmapper_go_evidence"]))
    #HMM specific settings
    if params["Eggnogmapper_method"] == "hmm":
        pcall.append("--database {}".format(params["Eggnogmapper_hmmdb"]))
        pcall.append("--dbtype {}".format(params["Eggnogmapper_dbtype"]))
        pcall.append("--qtype {}".format(params["Eggnogmapper_qtype"]))
        pcall.append("--hmm_maxhits {}".format(params["Eggnogmapper_hmm_maxhits"]))
        pcall.append("--hmm_evalue {}".format(params["Eggnogmapper_hmm_evalue"]))
        pcall.append("--hmm_score {}".format(params["Eggnogmapper_hmm_score"]))
        pcall.append("--hmm_maxseqlen {}".format(params["Eggnogmapper_hmm_maxseqlen"]))
        pcall.append("--hmm_qcov {}".format(params["Eggnogmapper_hmm_qcov"]))
        pcall.append("--hmm_Z {}".format(params["Eggnogmapper_z"]))
    #DIAMOND specific settings
    if params["Eggnogmapper_method"] == "diamond":
        if params["Eggnogmapper_dmnd_db"] != "false":
            pcall.append("--dmnd_db {}".format(params["Eggnogmapper_dmnd_db"]))
        if params["Eggnogmapper_matrix"] != "false":
            pcall.append("--matrix".format(params["Eggnogmapper_matrix"]))
        if params["Eggnogmapper_gapopen"] != "false":
            pcall.append("--gapopen".format(params["Eggnogmapper_gapopen"]))
        if params["Eggnogmapper_gapextend"] != "false":
            pcall.append("--gapextend".format(params["Eggnogmapper_gapextend"]))
    #general settings
    pcall.append("--seed_ortholog_evalue {}".format(params["Eggnogmapper_seed_ortholog_evalue"]))
    pcall.append("--seed_ortholog_score {}".format(params["Eggnogmapper_seed_ortholog_score"]))
    if params["Eggnogmapper_override"] != "false":
            pcall.append("--override")
    if params["Eggnogmapper_no_refine"] != "false":
            pcall.append("--no_refine")
    pcall.append("--no_annot")
    pcall.append("--no_file_comments")
    if params["Eggnogmapper_no_search"] != "false":
            pcall.append("--no_search")
    if params["Eggnogmapper_keep_mapping_files"] != "false":
            pcall.append("--keep_mapping_files")
    if params["Eggnogmapper_translate"] != "false":
            pcall.append("--translate")
    if params["Eggnogmapper_usemem"] != "false":
            pcall.append("--usemem")
    pcall.append("--cpu {}".format(params["Eggnogmapper_threads"]))
    return(" ".join(pcall))


'''
Function to build call to diamond
'''
def runDiamond(infile,outfile,params):
    dcall = ["diamond"]
    #blast mode, input, database, threads and output
    dcall.append("{} -q {} -d {} -p {} -o {}".format(params["Diamond_command"],infile,params["Diamond_db"],params["Diamond_threads"],outfile))
    #general options
    dcall.append("--outfmt {}".format(params["Diamond_outfmt"]))
    if params["Diamond_verbose"] != "false":
        dcall.append("--verbose")
    if params["Diamond_debug"] != "false":
        dcall.append("--debug")
    if params["Diamond_quiet"] != "false":
        dcall.append("--quiet")
    if params["Diamond_strand"] != "false":
        dcall.append("--strand {}".format(params["Diamond_strand"]))
    if params["Diamond_max_target_seqs"] != "false":
        dcall.append("--max-target-seqs {}".format(params["Diamond_max_target_seqs"]))
    if params["Diamond_top"] != "false":
        dcall.append("--top {}".format(params["Diamond_top"]))
    if params["Diamond_range_culling"] != "false":
        dcall.append("--range_culling")
    dcall.append("--evalue {}".format(params["Diamond_evalue"]))
    if params["Diamond_min_score"] != "false":
        dcall.append("--min-score {}".format(params["Diamond_min_score"]))
    if params["Diamond_id"] != "false":
        dcall.append("--id {}".format(params["Diamond_id"]))
    if params["Diamond_query_cover"] != "false":
        dcall.append("--query_cover {}".format(params["Diamond_query_cover"]))
    if params["Diamond_subject_cover"] != "false":
        dcall.append("--subject-cover {}".format(params["Diamond_subject_cover"]))
    if params["Diamond_more_sensitive"] != "false":
        dcall.append("--more-sensitive")
    elif params["Diamond_sensitive"] != "false":
        dcall.append("--sensitive")
    if params["Diamond_block_size"] != "false":
        dcall.append("--block-size {}".format(params["Diamond_block_size"]))
    if params["Diamond_index_chunks"] != "false":
        dcall.append("--index-chunks {}".format(params["Diamond_index_chunks"]))
    if params["Diamond_gapopen"] != "false":
        dcall.append("--gapopen {}".format(params["Diamond_gapopen"]))
    if params["Diamond_gapextend"] != "false":
        dcall.append("--gapextend {}".format(params["Diamond_gapextend"]))
    if params["Diamond_frameshift"] != "false":
        dcall.append("--frameshift {}".format(params["Diamond_frameshift"]))
    if params["Diamond_matrix"] != "false":
        dcall.append("--matrix {}".format(params["Diamond_matrix"]))
    if params["Diamond_comp_based_stats"] != "false":
        dcall.append("--comp-based-stats")
    if params["Diamond_masking"] != "false":
        dcall.append("--masking")
    if params["Diamond_no_self_hits"] != "false":
        dcall.append("--no-self-hits")
    if params["Diamond_taxonmap"] != "false":
        dcall.append("--taxonmap {}".format(params["Diamond_taxonmap"]))
    if params["Diamond_taxonnodes"] != "false":
        dcall.append("--taxonnodes {}".format(params["Diamond_taxonnodes"]))
    if params["Diamond_taxonlist"] != "false":
        dcall.append("--taxonlist {}".format(params["Diamond_taxonlist"]))
    return(" ".join(dcall))
    

'''
function to build call to blast2lca
'''

def runBlast2Lca(infile,outfile,params):
    #input, input format, output
    bcall = ["blast2lca -i {} -f {} -o {}".format(infile,params["Blast2lca_informat"],outfile)]
    #taxonomy lca options
    bcall.append("-sr {}".format(params["Blast2lca_sr"]))
    bcall.append("-oro {}".format(params["Blast2lca_oro"]))
    bcall.append("-tid {}".format(params["Blast2lca_tid"]))
    bcall.append("-ms {}".format(params["Blast2lca_ms"]))
    bcall.append("-me {}".format(params["Blast2lca_me"]))
    bcall.append("-top {}".format(params["Blast2lca_top"]))
    bcall.append("-mid {}".format(params["Blast2lca_mid"]))
    bcall.append("-tn {}".format(params["Blast2lca_tn"]))
    bcall.append("-a2t {}".format(params["Blast2lca_a2t"]))
    #add kegg options if true
    if params["Blast2lca_k"] != "false":
        bcall.append("-k -kr {} -a2kegg {}".format(params["Blast2lca_kr"],params["Blast2lca_a2kegg"]))
        if params["Blast2lca_ktp"] != "false":
            bcall.append("+ktp")
        #set out file for kegg
        kegout = os.path.basename(outfile).replace(".taxonomic.annotations",".megankegg.annotations")
        kegout = "functional_annotations.dir/"+kegout
        bcall.append("-ko {}".format(kegout))
    #other options
    bcall.append("-fwa {}".format(params["Blast2lca_fwa"]))
    bcall.append("-v {}".format(params["Blast2lca_v"]))
    return(" ".join(bcall))
