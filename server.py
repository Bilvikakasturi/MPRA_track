from pathlib import Path
from zipfile import ZipFile
import logging.handlers
import numpy as np
import pyranges as pr
import gzip
from Bio.Seq import Seq
from datetime import datetime
import logging

from biomart import BiomartServer
from flask import Flask, render_template, request, redirect, current_app, send_file, abort, after_this_request
from werkzeug.utils import secure_filename, send_from_directory
import pandas as pd

# snp_select_itr_enh
import os

app = Flask(__name__)
UPLOAD_FOLDER = './uploads'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


@app.route('/')
def index():
    return render_template('upload.html')
    # return render_template('test.html')


@app.route('/uploader', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        f = request.files['file']

        file_name = ''
        if f:
            file_name = f.filename
            f.save(secure_filename(f.filename))
        global utr_enh
        utr_enh = request.form.get("utr_enh")

        true_value = request.form.get("true_value")
        X = request.form.get("X")
        should_roll_over = os.path.isfile('log.txt')

        if should_roll_over:  # if log already exists, roll over!
            with open("log.txt", 'w') as file:
                file.write('readme')
                file.write('\n')
        logging.basicConfig(filename='log.txt', level=logging.INFO,
                            format=f'%(asctime)s %(levelname)s %(name)s %(threadName)s : %(message)s')
        log = logging.getLogger('werkzeug')
        log.setLevel(logging.ERROR)
        logging.info("Process started")
        try:
            snp_select_itr_enh(file_name, utr_enh, true_value, X)
            if utr_enh == 'BothUE':
                snp_select_itr_enh(file_name, "UTR", true_value, X, 1)

        except Exception as e:
            logging.error(e)
            logging.info("The program has stopped abruptly and encountered an error:")
            logging.info(e)
            SudZipWrite()
        else:
            # pass
            logging.info("The program executed sucessfully")
        return render_template('download.html')


@app.errorhandler(400)
def too_large(e):
    return "Invalid file", 400


@app.route('/stream')
def stream():
    my_file = Path("log.txt")
    if my_file.is_file():
        with open("log.txt", "r") as f:
            content = f.read()
            print(content)
    else:
        content = " "
    # as you see, file is loaded only once, no loop here, (loop is on frontend side)
    return app.response_class(content, mimetype='text/plain')


@app.route('/download')
def download_files():
    @after_this_request
    def remove_file(response):
        try:
            os.remove('SUD_files.zip')
        except Exception as error:
            app.logger.error("Error removing or closing downloaded file handle", error)
        return response

    return send_file('SUD_files.zip',
                     mimetype='zip',
                     attachment_filename='SUD_files.zip',
                     as_attachment=True)


def convert_bool_str(list):
    return ''.join(chr(ord('A') + i) if b else ' ' for i, b in enumerate(list))


def snp_select_itr_enh(file_name, utr_or_enhance, true_value, X, flag=0):
    server = BiomartServer("http://uswest.ensembl.org/biomart")
    hsapiens = server.databases['ENSEMBL_MART_SNP'].datasets['hsapiens_snp']
    RSID_values = request.form.get("RSID")
    rsquare = float(request.form.get("r2"))
    data_rsid = RSID_values.split()
    Differ = request.form.get("Differ")
    DifferEnh = request.form.get("DifferEnh")
    X = float(X)

    AF_EUR = "chk1" in request.form
    AF_EAS = "chk2" in request.form
    AF_AMR = "chk3" in request.form
    AF_SAS = "chk4" in request.form
    AF_AFR = "chk5" in request.form
    chk_list = [AF_EUR, AF_EAS, AF_AMR, AF_SAS, AF_AFR]
    print(chk_list, RSID_values, X, data_rsid)
    chk_bool_str = convert_bool_str(chk_list)
    flank = request.form.get("flank")
    flankEnh = request.form.get("flankEnh")
    if flank:
        flank = int(flank)
    p1 = request.form.get("p1")
    p2 = request.form.get("p2")
    p1Enh = request.form.get("p1Enh")
    p2Enh = request.form.get("p2Enh")
    chkENH = False
    chkUTR = False
    if utr_or_enhance == "UTR":
        chkUTR = "chkUTR" in request.form
        if not chkENH:
            flank = 25
            p1 = 'CCGTGTAATTCTAGGAGCTC'
            p2 = 'CGTTCTAGAGTCGGGGCGG'
            Differ = ''
    elif utr_or_enhance == "Enhancer":
        chkENH = "chkENH" in request.form
    elif utr_or_enhance == "BothUE":
        chkENH = "chkENH" in request.form
        chkUTR = "chkUTR" in request.form
    else:
        pass
    if utr_or_enhance == "Enhancer" or utr_or_enhance == "BothUE":
        if not chkENH:
            flank = 66
            p1 = 'ATCTAGAGCATGCACCGG'
            p2 = 'TCGACGAATTCGGCCGGC'
            DifferEnh = ''
        else:
            flank = int(flankEnh)
            p1 = p1Enh
            p2 = p2Enh
            Differ = DifferEnh
    print("values are", p1, p2, flank, utr_or_enhance)
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    if Differ or DifferEnh:
        flag = 0
    if not flag:
        opt_file_path = "CUD_01_snp_lds.txt"
        output_file = open(opt_file_path, "w")
        if file_name:
            data_1 = pd.read_csv(file_name, sep='\n', header=None)
            # if file has 4 columns with rsids column SNP
            # data_1 = pd.read_csv('rsids.txt',header=None, sep=' ')#if file type is only rsids
            data_1.drop_duplicates()
            data_1 = data_1.dropna(axis=1)
            # data_1.columns=['SNP']
            data_1[0] = data_1[0].str.strip()
            unique_snp = data_1[0].unique()
            if data_rsid:
                data_rsid_list = [item.split(',') for item in data_rsid]
                data_rsid_flat = [item for l in data_rsid_list for item in l]
                unique_snp = np.concatenate((unique_snp, data_rsid_flat))
            unique_snp = unique_snp.tolist()
        else:
            unique_snp = data_rsid
        filename = "LD_EUR.tsv.gz"

        with gzip.open(filename, 'rt') as f:
            f.readline()
            for line in f:
                x = line.split()
                if x[0] in unique_snp:
                    output_file.write(line)
                    unique_snp.remove(x[0])
                if not unique_snp:
                    output_file.close()
                    break

        def convert_lst_string(s):
            for ele in s:
                s1 = "".join(c for c in ele if c.isalnum())
            return s1

        # data_2 = pd.read_csv('CUD_01_snp_lds.txt', header=None)
        # data_2.columns = ['Lead_SNP', 'b']
        #
        # merged_df = pd.merge(data_1, data_2, on=['Lead_SNP'], how='inner')
        # merged_df = merged_df[['Lead_SNP', 'b']]
        if os.stat('CUD_01_snp_lds.txt').st_size == 0:
            logging.info("The rsid's provided doesnot return any data please check and provide new rsid values")
            SudZipWrite()
        merged_df = pd.read_csv('CUD_01_snp_lds.txt', header=None, sep='\t')
        merged_df.columns = ['Lead_SNP', 'b']
        final_SNP_list = []

        def extract_snp_greater_pt_eight():
            for key, value in merged_df.iteritems():
                if key == 'b':
                    for first_line in value:
                        f = []
                        snp_split = first_line.split(";")
                        for j in snp_split:
                            val = j.split(',')
                            if float(val[1]) >= rsquare:
                                f.append((val[0]))
                        final_SNP_list.append(f)
            merged_df['final_SNP_list'] = final_SNP_list
            merged_df.to_csv('withcolumn_b.txt', sep=' ', index=None)

        extract_snp_greater_pt_eight()
        merged_df_copy = merged_df.copy()
        merged_df_copy.drop('b', axis=1, inplace=True)
        merged_df_copy.to_csv('drop_b.txt', index=None, header=None, sep=' ')
        Lead_SNP = merged_df_copy['Lead_SNP']
        headers = ['Lead_SNP']

        def extract_bedfiles_fromSNP():
            global x
            with open('my_rsids.txt', 'w') as output_file:
                with open("drop_b.txt", 'r') as data_file:
                    for line in data_file:
                        data = line.split()
                        output_file.writelines("\n\n")
                        for rs in data:
                            data1 = rs.split(" ")
                            string_data = convert_lst_string(data1)
                            # s1 = "".join(c for c in data1 if c.isalnum())
                            output_file.writelines(string_data)
                            output_file.writelines("\n")
                data_file.close()
            output_file.close()

        extract_bedfiles_fromSNP()
        big_bed_tools_output = './bigBedNamedItems -nameFile http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153.bb my_rsids.txt ' \
                               'output_my_rsids_location.txt '
        os.system(big_bed_tools_output)
        file_path = 'output_my_rsids_location.txt'

        if os.stat(file_path).st_size == 0:
            logging.info("Error in the SNP files input not given correctly, locations of the SNP's not found")
            SudZipWrite()
            return
        else:
            logging.info("BedNameditems location files created succesfully")
        output_rsids = pd.read_csv('output_my_rsids_location.txt', sep='\t', header=None)
        first_4_columns = output_rsids.iloc[:, :4]
        first_4_columns.columns = ['CHR', 'min_range', 'max_range', 'SNP']
        onlyrsids = pd.read_csv('my_rsids.txt', sep='\t', header=None)
        onlyrsids = onlyrsids.drop_duplicates().reset_index(drop=True)
        onlyrsids.columns = ['SNP']
        unique_vals = onlyrsids[~onlyrsids.SNP.isin(first_4_columns.SNP)].append(
            first_4_columns[~first_4_columns.SNP.isin(onlyrsids.SNP)],
            ignore_index=True)['SNP'].to_list()
        biomart_op = biomart_Synonumvariant(hsapiens, unique_vals, 'snp_synonym_filter')
        # to check how many snp's returned
        synonym_filter = biomart_op['SNP'].to_list()
        list_snp_selection = list(set(synonym_filter).symmetric_difference(set(unique_vals)))
        if list_snp_selection:
            biomart_op1 = biomart_Synonumvariant(hsapiens, list_snp_selection, 'snp_filter')
            biomart_op = biomart_op.append(biomart_op1, ignore_index=True)
        biomart_snp = biomart_op.groupby('SNP').agg(
            {'CHR': 'min', 'min_range': 'min', 'max_range': 'max'}).reset_index()
        biomart_snp['CHR'] = 'chr' + biomart_op['CHR'].astype(str)
        biomart_snp = biomart_snp[['CHR', 'min_range', 'max_range', 'SNP']]

        new_4_rsids = first_4_columns.append(biomart_snp, ignore_index=True)
        new_4_rsids = new_4_rsids.astype({'min_range': 'int', 'max_range': 'int'})
        i = (new_4_rsids[new_4_rsids.CHR.str.contains(r'[_]')]).index
        new_4_rsids = new_4_rsids.drop(labels=i, axis=0)
        opt_rsids = pd.read_csv('drop_b.txt', sep=' ', header=None)
        opt_rsids.columns = ['lead_SNP', 'SNP']
        new_df = opt_rsids
        new_df['new_column'] = new_df['lead_SNP'] + ',' + new_df['SNP']
        # min_range_list = []
        failed_list = []
        every_line_list = []
        # every_char_list = []
        chr_range_lis = []
        new_chr_list = []
        chr_count_list = []
        every_line_list_max = []
        # chr_count_list_key = []
        f = []

        def extr_Chr_min_max(Differ):
            global f, i

            for j in value:
                minimum_range_list = []
                maximum_range_list = []
                f = []
                val = j.split(',')
                try:
                    main_min_range = (new_4_rsids['min_range'].loc[new_4_rsids['SNP'] == val[0]].values.tolist()[0])
                    main_max_range = (new_4_rsids['max_range'].loc[new_4_rsids['SNP'] == val[0]].values.tolist()[0])
                except Exception as x:

                    main_min_range = 'empty'
                    main_max_range = 'empty'
                # print(main_min_range)
                # print(main_max_range)
                assign_min_max_chr(f, maximum_range_list, minimum_range_list, val)

                chr_range_lis.append(f)
                dict_of_counts = {item: f.count(item) for item in f}
                chr_count_list.append(dict_of_counts)
                res = []
                while ("" in f):
                    f.remove("")
                for i in f:
                    if i not in res:
                        res.append(i)
                str_chr = ''.join(res)
                new_chr_list.append(str_chr)
                every_line = min(x for x in minimum_range_list if x is not None)

                every_line_max = max(y for y in maximum_range_list if y is not None)
                if Differ and main_min_range != 'empty' and main_max_range != 'empty':
                    Differ = int(Differ)

                    if (every_line_max - main_max_range > Differ):
                        every_line_max = main_max_range + Differ
                    if (main_min_range - every_line > Differ):
                        every_line = main_min_range - Differ
                every_line_list.append(every_line)
                every_line_list_max.append(every_line_max)

        def assign_min_max_chr(f, maximum_range_list, minimum_range_list, val):
            for x in val:
                minimum_range = None
                maximum_range = 0
                chr_range = ''
                s1 = "".join(c for c in x if c.isalnum())
                try:
                    if s1:
                        chr_range = new_4_rsids['CHR'].loc[new_4_rsids['SNP'] == s1].values.tolist()[0]
                    else:
                        chr_range = ''
                    minimum_range = (new_4_rsids['min_range'].loc[new_4_rsids['SNP'] == s1].values.tolist()[0])
                    maximum_range = (new_4_rsids['max_range'].loc[new_4_rsids['SNP'] == s1].values.tolist()[0])
                    chr_range = new_4_rsids['CHR'].loc[new_4_rsids['SNP'] == s1].values.tolist()[0]

                except Exception as x:
                    if s1:
                        failed_list.append(s1)
                minimum_range_list.append(minimum_range)
                maximum_range_list.append(maximum_range)

                f.append(chr_range)
        print("The list is",failed_list)

        failed_lit = pd.DataFrame(failed_list)
        failed_lit.to_csv('testingfailed.txt')
        failed_lit = failed_lit.dropna(axis=1)
        print("The dataframe of failed list is",failed_lit)
        if not failed_lit[failed_lit.isin([0])].empty:
            if flag == 1:
                failed_lit.to_csv('failedSNPListUTR.txt')
            else:
                failed_lit.to_csv('failedSNPList.txt')
        for key, value in new_df.iteritems():
            if key == 'new_column':
                extr_Chr_min_max(Differ)
        opt_rsids['CHR'] = new_chr_list
        opt_rsids['min'] = every_line_list
        opt_rsids['max'] = every_line_list_max
        opt_rsids['count_char'] = chr_count_list
        opt_rsids.to_csv('final_output.csv', index=None, sep=',')
        opt_rsids = opt_rsids.reindex(columns=['CHR', 'min', 'max', 'lead_SNP', 'SNP', 'new_column', 'count_char'])
        first_4_columns_leadSNP = opt_rsids.iloc[:, :4]
        new_4_rsids_frst = first_4_columns_leadSNP
        new_4_rsids_frst.to_csv('SNP_LD.bed', header=None, index=None, sep='\t')

    def sum_true(row):
        row1 = sum(row)
        return row1

    if utr_or_enhance == 'UTR':
        name_file = '3UTR'
        path2 = '/Users/bilvikakasturi/PycharmProjects/pythonProject1/MPRA/hg38.ensGene.3UTR.gtf.bed'
        gr3 = pr.read_gtf(path2)
    else:
        name_file = 'enh'
        path2 = '/Users/bilvikakasturi/PycharmProjects/pythonProject1/MPRA/GRCh38-ELS.bed'
        gr3 = pr.read_bed(path2)
    f_gtf = 0
    f_bed = 0

    if chkUTR:
        f_gtf = request.files['filegtfbed']
    if chkENH:
        f_bed = request.files['filebed']

    if f_gtf:
        file_name_gtf = f_gtf.filename
        f_gtf.save(secure_filename(f_gtf.filename))
        path2 = file_name_gtf
        gr3 = pr.read_gtf(path2)
    if f_bed:
        file_name_bed = f_bed.filename
        f_bed.save(secure_filename(f_bed.filename))
        path2 = file_name_bed
        gr3 = pr.read_bed(path2)

    # segment_all_vcf = 'bedtools intersect -header -a SNP_LD.bed -b '+path2+' > intersect.bed'
    # os.system(segment_all_vcf)
    # gr3 = pr.read_bed(path2)
    def intersect_gz():
        path = '/Users/bilvikakasturi/PycharmProjects/pythonProject1/MPRA/SNP_LD.bed'
        gr = pr.read_bed(path)
        grinn = gr.intersect(gr3)
        # gr4.to_csv('intersect.bed', sep='\t')
        return grinn

    grintersect = intersect_gz()
    grintersect.to_csv('intersect.bed', sep='\t')
    intersect_path = 'intersect.bed'
    if os.stat(intersect_path).st_size == 0 or (not grintersect):
        logging.info(
            "Intersection of the files with GRCh38-ELS.bed for enhancer or hg38.ensGene.3UTR.gtf.bed for UTR is not "
            "successful")
        SudZipWrite()
        return
    else:
        logging.info(
            "Intersection of the files with GRCh38-ELS.bed for enhancer or hg38.ensGene.3UTR.gtf.bed for UTR bed "
            "files creation is sucessfull")
    df = pd.read_csv('intersect.bed', sep='\t')
    unique_chromosomes = df.Chromosome.unique()

    def extract_indv_chrfiles():
        for i in unique_chromosomes:
            chr_file = 'CUD_ldreg_' + i + '.bed'
            df[df.Chromosome == i].to_csv(chr_file, index=None, header=None, sep='\t')

    def bed_to_vcf():
        for c in unique_chromosomes:
            vcf_files = 'bcftools view -G -H -R CUD_ldreg_' + c + \
                        '.bed ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working' \
                        '/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_' + \
                        c + '.recalibrated_variants.annotated.vcf.gz  > CUD_1000g_' + c + '.vcf'
            os.system('conda activate bcftools_env')
            os.system(vcf_files)

    def concat_header_allvcf():
        list_vcf = []
        header = 'bcftools view -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage' \
                 '/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr1' \
                 '.recalibrated_variants.annotated.vcf.gz  > CUD_1000g_header.vcf '
        os.system(header)
        list_vcf.append('CUD_1000g_header.vcf')
        for i in unique_chromosomes:
            vcf_file = 'CUD_1000g_' + i + '.vcf'
            list_vcf.append(vcf_file)
        with open("CUD_1000g_all.vcf", "w") as new_file:
            for name in list_vcf:
                with open(name) as file:
                    for line in file:
                        new_file.write(line)
        vcf_files = 'CUD_1000g_all.vcf'
        if os.stat(vcf_files).st_size == 0:
            logging.info("VCF files not created error in bed files")
            SudZipWrite()
            return
        else:
            logging.info("VCF files created succesfully ")

    def extract_desired_frmallvcf():
        segment_all_vcf = 'bcftools view -v snps CUD_1000g_all.vcf |gzip > CUD_1000g_snps.vcf.gz'
        os.system(segment_all_vcf)
        extract_desired = "bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF_EUR\t%AF_EAS\t%AF_AMR\t%AF_SAS\t%AF_AFR" \
                          "\n' CUD_1000g_snps.vcf.gz |gzip > CUD_06_1000g_snps_combo_maf.txt.gz "
        os.system(extract_desired)

    extract_indv_chrfiles()
    after_bed = datetime.now()
    after_time = after_bed.strftime("%H:%M:%S")
    logging.info("After Bed file time =")
    logging.info(after_time)
    bed_to_vcf()
    concat_header_allvcf()
    after_vcf = datetime.now()
    after_time = after_vcf.strftime("%H:%M:%S")
    logging.info("After VCF files extraction time =")
    logging.info(after_time)
    extract_desired_frmallvcf()
    x_df = pd.read_csv('CUD_06_1000g_snps_combo_maf.txt.gz', compression='gzip', header=None, sep='\t')
    # true_value = 2
    x_df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'AF_EUR', 'AF_EAS', 'AF_AMR', 'AF_SAS', 'AF_AFR']
    x_df['AF_EUR'] = x_df['AF_EUR'].replace(['.'], 0)
    x_df = x_df[x_df.AF_EUR != '.']
    x_df.AF_EUR = x_df.AF_EUR.astype(float)

    x_df['true_EUR'] = (x_df['AF_EUR'].between(X, 1 - X))
    x_df['true_EAS'] = (x_df['AF_EAS'].between(X, 1 - X))
    x_df['true_AMR'] = (x_df['AF_AMR'].between(X, 1 - X))
    x_df['true_SAS'] = (x_df['AF_SAS'].between(X, 1 - X))
    x_df['true_AFR'] = (x_df['AF_AFR'].between(X, 1 - X))
    x_df["true"] = x_df[["true_EUR", "true_EAS", "true_AMR", "true_SAS", "true_AFR"]].values.tolist()
    x_df['true_count'] = x_df['true'].apply(lambda row: sum_true(row))
    x_df['bool_str'] = x_df['true'].apply(convert_bool_str)
    if true_value == 'Select':
        cud_u3_max = x_df[x_df['bool_str'] == chk_bool_str]
    else:
        cud_u3_max = x_df.loc[x_df['true_count'] >= int(true_value)]
    # cud_u3_max = x_df.loc[(x_df['AF_EUR'].between(X, 1-X)) | (x_df['AF_EAS'].between(X, 1-X)) | (
    #     x_df['AF_AMR'].between(X, 1-X)) | (x_df['AF_SAS'].between(X, 1-X)) | (x_df['AF_AFR'].between(X, 1-X))]
    cud_u3_max.to_csv('cud_u3_maf_uniq.txt.gz', index=None, sep='\t')  # verify if this file is needed
    if cud_u3_max.empty:
        logging.info('values which lie between the AF is null so no files generated!!!')
        SudZipWrite()
        return
    after_truevalue = datetime.now()
    after_time = after_truevalue.strftime("%H:%M:%S")
    print("After after_truevalue =", after_time)
    selected_columns = cud_u3_max[['AF_EUR', 'AF_EAS', 'AF_AMR', 'AF_SAS', 'AF_AFR']]
    SUD_dataframe = pd.DataFrame()
    SUD_dataframe['CHROM'] = cud_u3_max['CHROM']
    SUD_dataframe['POS1'] = cud_u3_max['POS'] - 1 - flank
    SUD_dataframe['POS2'] = cud_u3_max['POS'] + flank
    # SUD_dataframe = selected_columns.copy()
    cols = ['CHROM', 'POS', 'REF', 'ALT']
    SUD_dataframe['combined'] = cud_u3_max[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    SUD_dataframe['combined'] = 'SUD_pass_' + SUD_dataframe['combined']
    SUD_dataframe.to_csv('SUD_pass_oligoloc.bed', header=None, index=None, sep='\t')
    SUD_dataframe = SUD_dataframe.join(selected_columns)
    # SUD_dataframe.to_csv('SUD_pass_oligoloc.bed', header=None, index=None, sep='\t')
    # generating .fa file using bedtools
    bedtools_execution = "bedtools getfasta -fi hg38.fa -bed SUD_pass_oligoloc.bed -name -fo test1.fa"
    os.system(bedtools_execution)
    # SUD_dataframe.to_csv('SUD_pass_ol.bed', index=None, sep='\t')
    sud_oligo_df = pd.read_csv('SUD_pass_oligoloc.bed', sep='\t')
    # sud_oligo_df = pd.read_csv('SUD_pass_ol.bed', sep='\t')
    sud_oligo_df.columns = ['CHR', 'POS1', 'POS2', 'name']
    extract_position = sud_oligo_df['name'].str.split('_', expand=True)
    sud_oligo_df['POS1'] = extract_position[3]
    sud_oligo_df['POS1'] = sud_oligo_df['POS1'].astype(int)
    sud_oligo_df['POS2'] = sud_oligo_df['POS1'] + 1
    sud_oligo_df.to_csv('SUD_pass_snploc.bed', index=None, sep='\t')
    bedtools_rsid = "bedtools intersect -a SUD_pass_snploc.bed -b SNP_LD.bed -wo > SUD_pass_rsid.txt"
    os.system(bedtools_rsid)
    new_sud_pass_3UTR = []

    def three_UTR_region():
        global sud_pass_df

        bedtools_3UTR = "bedtools intersect -a SUD_pass_snploc.bed -b hg38.ensGene.3UTR.gtf.gz -wo > SUD_pass_3UTR.txt"

        os.system(bedtools_3UTR)

        sud_pass_df = pd.read_csv('SUD_pass_3UTR.txt', header=None, sep='\t')
        SUD_rsid = pd.read_csv('SUD_pass_rsid.txt', header=None, sep='\t')
        rsid_name = SUD_rsid.iloc[:, [3, 7]]
        rsid_name.to_csv('Name_rsid.text', index=None, header=None, sep='\t')
        sud_pass_df.columns = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', 'gene', '14']
        sud_pass_df['name'] = '>' + sud_pass_df['4'] + ':' + sud_pass_df['11']

    def enhancer_region():
        global new_sud_pass_enhancer
        new_sud_pass_enhancer = sud_oligo_df[['name']]
        new_sud_pass_enhancer['plus'] = sud_oligo_df['name'] + ':+'
        new_sud_pass_enhancer['minus'] = sud_oligo_df['name'] + ':-'
        new_sud_pass_enhancer = pd.melt(new_sud_pass_enhancer, id_vars=['name'], value_vars=['plus', 'minus'])
        new_sud_pass_enhancer = new_sud_pass_enhancer.drop(['name', 'variable'], axis=1)
        new_sud_pass_enhancer.columns = ['name']
        new_sud_pass_enhancer['name'] = '>' + new_sud_pass_enhancer['name']

    if utr_or_enhance == 'UTR':
        three_UTR_region()
        new_sud_pass_3UTR = sud_pass_df[['name']]
        # generate new file with sudpassdf[name]and geneid
        sign_geneid = sud_pass_df[['name']].copy()
        sign_geneid['geneid'] = (sud_pass_df.gene.str.split(pat='"', expand=True))[1]
        sign_geneid = sign_geneid.drop_duplicates().reset_index(drop=True)
        sign_geneid.to_csv('chr_sign_geneid.txt', index=None, sep='\t')
    # if enhancer
    else:
        enhancer_region()
        new_sud_pass_3UTR = new_sud_pass_enhancer[['name']]
    drop_duplicates = new_sud_pass_3UTR.drop_duplicates()
    drop_duplicates = drop_duplicates.reset_index(drop=True)
    # drop_duplicates.to_csv('sud_u3_strandinfo.txt', index=None, sep='\t')
    drop_duplicates[['to_verify', 'SIGN']] = drop_duplicates['name'].str.split('\:', expand=True)
    # convert test1.fa to csv file
    list1 = []
    list2 = []
    with open('test1.fa', mode='r', newline='') as f:
        lines = f.readlines()
        for number, line in enumerate(lines):
            if number % 2:
                list2.append(line)
            if not number % 2:
                list1.append(line)
    data_tuples = list(zip(list1, list2))
    name_sequence = pd.DataFrame(data_tuples, columns=['name', 'ref_sequence'])
    name_sequence = name_sequence.replace('\n', '', regex=True)
    Alt_alpha_sequence = (name_sequence.name.str.split(pat='::', expand=True))[0]
    Alt_alpha = Alt_alpha_sequence.str.split(pat='_', expand=True)[5]
    name_sequence['combined'] = Alt_alpha_sequence
    SUD_dataframe['combined'] = '>' + SUD_dataframe['combined']
    name_sequence = pd.merge(name_sequence, SUD_dataframe, on=['combined'], how='inner')
    # making the sequence column strings to upper case
    name_sequence['ref_sequence'] = name_sequence['ref_sequence'].str.upper()
    name_sequence['alt_sequence'] = name_sequence['ref_sequence'].str[:flank] + Alt_alpha + name_sequence[
                                                                                                'ref_sequence'].str[
                                                                                            flank + 1:]
    # if utr_or_enhance!='UTR':
    #     name_sequence=pd.concat([name_sequence]*2, ignore_index=True)
    name_sequence['to_verify'] = name_sequence['name'].str.split('\::').str[0]
    new_dataframe_seq = drop_duplicates.merge(name_sequence, on=['to_verify'], suffixes=['_df1', '_df2'])
    new_dataframe_seq = new_dataframe_seq.drop(['to_verify', 'name_df2'], axis=1)

    def compliment_DNA(row):
        row1 = Seq(row).reverse_complement()
        return row1

    split_name_df = new_dataframe_seq['name_df1'].str.split(pat=":", expand=True)
    y1 = split_name_df[0].str.split(pat="_", expand=True)
    new_dataframe_seq['sequence_alternative'] = new_dataframe_seq['ref_sequence'].apply(
        lambda row: compliment_DNA(row)).astype('string')
    new_dataframe_seq['alt_sequence_alternative'] = new_dataframe_seq['alt_sequence'].apply(
        lambda row: compliment_DNA(row)).astype('string')
    new_dataframe_seq.loc[new_dataframe_seq['SIGN'] != '+', 'ref_sequence'] = new_dataframe_seq['sequence_alternative']
    new_dataframe_seq.loc[new_dataframe_seq['SIGN'] == ' ', 'ref_sequence'] = 'data mismatch'
    new_dataframe_seq.loc[new_dataframe_seq['SIGN'] != '+', 'alt_sequence'] = new_dataframe_seq[
        'alt_sequence_alternative']
    new_dataframe_seq.loc[new_dataframe_seq['SIGN'] == ' ', 'alt_sequence'] = 'data mismatch'
    new_dataframe_seq["gREF"] = y1[4]
    new_dataframe_seq["gALT"] = y1[5]
    new_dataframe_seq["oALT"] = new_dataframe_seq["gALT"]
    new_dataframe_seq["oREF"] = new_dataframe_seq["gREF"]
    new_dataframe_seq['oREF_ALT'] = new_dataframe_seq['gREF'].apply(lambda row: compliment_DNA(row)).astype('string')
    new_dataframe_seq['oALT_ALT'] = new_dataframe_seq['gALT'].apply(lambda row: compliment_DNA(row)).astype('string')
    new_dataframe_seq.loc[new_dataframe_seq['SIGN'] != '+', 'oREF'] = new_dataframe_seq['oREF_ALT']
    new_dataframe_seq.loc[new_dataframe_seq['SIGN'] != '+', 'oALT'] = new_dataframe_seq['oALT_ALT']
    new_dataframe_seq = new_dataframe_seq.drop(
        ['sequence_alternative', 'alt_sequence_alternative', 'SIGN', 'oREF_ALT', 'oALT_ALT'], axis=1)
    # y2 = y1.apply(lambda y: '%s_%s_%s_%s_%s_%s' % (y[0], y[1], y[2], y[3], y[4], y[5]), axis=1)
    # new_dataframe_seq['Name'] = y2 + '_' + split_name_df[1]
    new_dataframe_seq['Name'] = split_name_df[0] + '_' + split_name_df[1]
    new_dataframe_seq = new_dataframe_seq.drop(['name_df1'], axis=1)
    new_dataframe_seq = new_dataframe_seq[
        ['Name', 'gREF', 'gALT', 'ref_sequence', 'alt_sequence', 'oREF', 'oALT', 'AF_EUR', 'AF_EAS', 'AF_AMR', 'AF_SAS',
         'AF_AFR']]
    new_dataframe_seq.to_csv('SUD_' + name_file + '_oligos.txt', index=None, sep='\t')
    oligos_AF = new_dataframe_seq
    new_dataframe_seq = new_dataframe_seq.drop(
        ['AF_EUR', 'AF_EAS', 'AF_AMR', 'AF_SAS', 'AF_AFR', 'gREF', 'gALT', 'oREF', 'oALT'], axis=1)
    sud_u3_oligoinfo = pd.melt(new_dataframe_seq, id_vars=["Name"],
                               var_name="SEQ_Behavior", value_name="SEQUENCE")
    sud_u3_oligoinfo['final_Name'] = sud_u3_oligoinfo['Name'].astype(str) + '_' + sud_u3_oligoinfo[
        'SEQ_Behavior'].astype(
        str)
    sud_u3_oligoinfo = sud_u3_oligoinfo.drop(['Name', 'SEQ_Behavior'], axis=1)
    sud_u3_oligoinfo = sud_u3_oligoinfo[['final_Name', 'SEQUENCE']]
    sud_u3_oligoinfo.to_csv('sud_' + name_file + '_oligoinfo.txt', index=None, sep='\t')
    oligos_AF['REF_SEQUENCE'] = p1 + oligos_AF['ref_sequence'].astype(str) + p2
    oligos_AF['ALT_SEQUENCE'] = p1 + oligos_AF['alt_sequence'].astype(str) + p2
    new_dataframe_seq = oligos_AF
    new_dataframe_seq = new_dataframe_seq.drop(['ref_sequence', 'alt_sequence'], axis=1)
    new_dataframe_seq.to_csv('SUD_' + name_file + '_oligos_p1_p2.txt', index=None, sep='\t')
    new_dataframe_seq = new_dataframe_seq.drop(
        ['AF_EUR', 'AF_EAS', 'AF_AMR', 'AF_SAS', 'AF_AFR', 'gREF', 'gALT', 'oREF', 'oALT'], axis=1)

    sud_u3_oligoinfo = pd.melt(new_dataframe_seq, id_vars=["Name"],
                               var_name="SEQ_Behavior", value_name="SEQUENCE")
    sud_u3_oligoinfo['final_Name'] = sud_u3_oligoinfo['Name'].astype(str) + '_' + sud_u3_oligoinfo[
        'SEQ_Behavior'].astype(
        str)
    sud_u3_oligoinfo = sud_u3_oligoinfo.drop(['Name', 'SEQ_Behavior'], axis=1)
    sud_u3_oligoinfo = sud_u3_oligoinfo[['final_Name', 'SEQUENCE']]
    sud_u3_oligoinfo.to_csv('sud_' + name_file + '_oligoinfo_p1_p2.txt', index=None, sep='\t')
    # logging.debug("Logging test...")
    # logging.info("The program is working as expected")
    # logging.warning("The program may not function properly")
    # logging.error("The program encountered an error")
    # logging.critical("The program crashed")
    if os.path.exists("failedSNPList.txt"):
        logging.info("There are few rs SNIP's which did not return any values, Please find them in failedSNPList.txt "
                     "file in ZIP File.")
        failed_exists = True
    else:
        failed_exists = False

    if flag == 1:
        with ZipFile("SUD_files.zip", "a") as newzip:
            newzip.write('sud_' + name_file + '_oligoinfo_p1_p2.txt')
            newzip.write('SUD_' + name_file + '_oligos_p1_p2.txt')
            newzip.write('SUD_' + name_file + '_oligos.txt')
            newzip.write('sud_' + name_file + '_oligoinfo.txt')
            if failed_exists:
                newzip.write('failedSNPListUTR.txt')
            newzip.write('log.txt')
            newzip.write('chr_sign_geneid.txt')
        return
    with ZipFile("SUD_files.zip", "w") as newzip:
        newzip.write('sud_' + name_file + '_oligoinfo_p1_p2.txt')
        newzip.write('SUD_' + name_file + '_oligos_p1_p2.txt')
        newzip.write('SUD_' + name_file + '_oligos.txt')
        newzip.write('sud_' + name_file + '_oligoinfo.txt')
        if failed_exists:
            newzip.write('failedSNPList.txt')
        if not utr_or_enhance == 'BothUE':
            newzip.write('log.txt')
        if name_file == '3UTR':
            newzip.write('chr_sign_geneid.txt')
    after = datetime.now()
    after_time = after.strftime("%H:%M:%S")
    print("After Time =", after_time)
    time_1 = datetime.strptime(current_time, "%H:%M:%S")
    time_2 = datetime.strptime(after_time, "%H:%M:%S")
    time_interval = time_2 - time_1
    print(time_interval)


def biomart_Synonumvariant(hsapiens, unique_vals, passed_filter):
    global x
    chunks = [unique_vals[x:x + 500] for x in range(0, len(unique_vals), 500)]
    tlist = []
    for x in chunks:
        response = hsapiens.search({
            'filters': {passed_filter: x},
            'attributes': ['chr_name', 'chrom_start', 'chrom_end', 'synonym_name']
        })
        # 'synonym_name',
        llist = []
        for line in response.iter_lines():
            x = (line.decode())
            llist.append(x.split())
        tlist = tlist + llist
    biomart_op = pd.DataFrame(tlist,
                              columns=['CHR', 'min_range', 'max_range', 'SNP'])
    return biomart_op


def SudZipWrite():
    with ZipFile("SUD_files.zip", "w") as newzip:
        newzip.write('log.txt')


if __name__ == '__main__':
    # app.config['TEMPLATES_AUTO_RELOAD'] = True
    app.run(debug=True, port=5000)
