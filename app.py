import io
import os
from json import loads
from os import listdir
from celery import Celery
import numpy as np
import pandas as pd
from numpy import int64, float64

from BioinformaticsLab.ComputationalBiology.data_analysis import all_features_calculator as fc
from BioinformaticsLab.ComputationalBiology.data_analysis.main_parser_features_calc import features_on_each_gene
from BioinformaticsLab.ComputationalBiology.file_utils.genbank_parser import read_genbank_file
from BioinformaticsLab.ComputationalBiology.data_analysis.gene_features_calculator import compute_gc_content
import boto3
from flask import Flask, send_from_directory, jsonify, request, json
from flask_cors import CORS, cross_origin
from random import randrange
from Bio import Entrez, SeqIO

import time
from werkzeug.utils import secure_filename

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from flask import Response

# b = back.Back()
app = Flask(__name__)
cors = CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'
# absoulte path for dir
cwd = os.getcwd()

aws_storage_bucket_name = "bio-upload-files"
ACCESS_ID = "AKIARDVFNR2ZY5JRSV7Z"
ACCESS_KEY = "HgkvAFPGms/KfGVUW/YIyA4cq+TopM1uaxVj2ocx"
s3_client = boto3.resource('s3', aws_access_key_id=ACCESS_ID, aws_secret_access_key= ACCESS_KEY)

path_to_input_file = os.path.join(cwd, "BioinformaticsLab", "data", "data_inputs", "GenBank")
path_to_pickle_files = os.path.join(cwd, "BioinformaticsLab", "data", "data_outputs")
# path = f"{cwd}/TranscriptionalFactors/data_inputs"
features_dict = {
    'GC CONTENT': 'GC_CONTENT',
    'DNA LENGTH': 'DNA_LENGTH',
    'GENE NAME': 'GENE_NAME',
    'GENE ID': '',
    'TYPE': 'TYPE',
    'PRODUCT TYPE': 'PRODUCT_TYPE',
    'STRAND': 'STRAND',
    'PRODUCT DESCRIPTION': 'PRODUCT_DESCRIPTION',
    'HYDROPHOBIC AA': 'HYDROPHOBIC_AA',
    'HYDROPHILIC AA': 'HYDROPHILIC_AA',
    'POLAR AA': 'POLAR_AA',
    'AROMATIC AA': 'AROMATIC_AA',
    'POSITIVE AA': 'POSITIVE_AA',
    'NEGATIVE AA': 'NEGATIVE_AA',
    'NONPOLAR AA': 'NONPOLAR_AA',
    'AA LENGTH': 'AA_LENGTH',
    'H1': 'H1',
    'H2': 'H2',
    'H3': 'H3',
    'V': 'V',
    'P1': 'P1',
    'P2': 'P2',
    'SASA': 'SASA',
    'NCI': 'NCI',
    'MASS': 'MASS',
    'PKA COOH': 'PKA_COOH',
    'PKA NH': 'PKA_NH',
    'PI': 'PI'

}
features_description = {
    'GC CONTENT': 'GC_CONTENT',
    'DNA LENGTH': 'DNA_LENGTH',
    'GENE NAME': 'GENE_NAME',
    'GENE ID': 'GENE ID',
    'TYPE': 'TYPE',
    'PRODUCT TYPE': 'PRODUCT_TYPE',
    'STRAND': 'STRAND',
    'PRODUCT DESCRIPTION': 'PRODUCT_DESCRIPTION',
    'HYDROPHOBIC AA': 'HYDROPHOBIC_AA',
    'HYDROPHILIC AA': 'HYDROPHILIC_AA',
    'POLAR AA': 'POLAR_AA',
    'AROMATIC AA': 'AROMATIC_AA',
    'POSITIVE AA': 'POSITIVE_AA',
    'NEGATIVE AA': 'NEGATIVE_AA',
    'NONPOLAR AA': 'NONPOLAR_AA',
    'AA LENGTH': 'AA_LENGTH',
    'H1': 'H1',
    'H2': 'H2',
    'H3': 'H3',
    'V': 'V',
    'P1': 'P1',
    'P2': 'P2',
    'SASA': 'SASA',
    'NCI': 'NCI',
    'MASS': 'MASS',
    'PKA COOH': 'PKA_COOH',
    'PKA NH': 'PKA_NH',
    'PI': 'PI'
}


genome_feature_list = ['GC_CONTENT', 'DNA_LENGTH']
gene_feature_list = ['GC_CONTENT', 'DNA_LENGTH']
general_feature_list = ['GENE_ID', 'GENE_NAME', 'TYPE', 'PRODUCT_TYPE', 'STRAND', 'PRODUCT_DESCRIPTION']
# protein_feature_list = ['HYDROPHOBIC_AA', 'HYDROPHILIC_AA', 'POLAR_AA', 'AROMATIC_AA', 'POSITIVE_AA', 'NEGATIVE_AA',
#                         'NONPOLAR_AA', 'AA_LENGTH', 'H1', 'H2', 'H3', 'V', 'P1', 'P2', 'SASA', 'NCI', 'MASS', 'PKA_COOH',
#                         'PKA_NH', 'PI']
protein_feature_list = ['HYDROPHOBIC_AA', 'HYDROPHILIC_AA', 'POLAR_AA', 'AROMATIC_AA', 'POSITIVE_AA', 'NEGATIVE_AA',
                        'NONPOLAR_AA', 'AA_LENGTH']
numeric_feature_list =['GC_CONTENT', 'DNA_LENGTH', 'HYDROPHOBIC_AA', 'HYDROPHILIC_AA', 'POLAR_AA', 'AROMATIC_AA',
                       'POSITIVE_AA', 'NEGATIVE_AA', 'NONPOLAR_AA', 'AA_LENGTH', 'H1', 'H2', 'H3', 'V', 'P1', 'P2',
                       'SASA', 'NCI', 'MASS', 'PKA_COOH', 'PKA_NH', 'PI']


def match_client_feature_to_df(feature_list_by_user):
    feature_list_by_user_match_server = list()
    for feature_name in feature_list_by_user:
        print(feature_name)
        if feature_name in features_dict:
            feature_list_by_user_match_server.append(features_dict[feature_name])
    return feature_list_by_user_match_server


@cross_origin()
# @app.route('/')
# def get_current_time():
#     return {'hello': b.printHello()}


@app.route('/api/features', methods=['GET'])
def feature():
    file_list_names = request.args.getlist('fileList[]')
    feature_list_by_user = request.args.getlist('featureList[]')
    print(feature_list_by_user)
    arr_match_feature_server = match_client_feature_to_df(feature_list_by_user)
    print("after match", arr_match_feature_server)
    path_to_pickle_files = {}
    to_return = dict()

    for file in file_list_names:

        to_return[file] = dict()
        fileName = os.path.splitext(file)[0]
        path_to_pickle_files[file] = './BioinformaticsLab/data/data_outputs/features_' + fileName + '.pickle'
        print(fileName)
        if len(check_existing_files(file)) != 2:
            features_on_each_gene(fileName)

        for fileName in path_to_pickle_files:
            data_frame_file = pd.read_pickle(path_to_pickle_files[fileName])
            full_data_features = dict()
            genome_dict = dict()
            gene_bank_file ='./BioinformaticsLab/data/data_inputs/GenBank/' + fileName
            with open(gene_bank_file, "r") as input_handle:
                gen = SeqIO.parse(input_handle, "genbank")
                record_gb = next(gen)
            genome_dict['Description'] = record_gb.description
            genome_dict['Publish date'] = record_gb.annotations['date']
            count_type= dict(data_frame_file['PRODUCT_TYPE'].value_counts())
            for type in count_type:
                genome_dict['Number of '+type] = int(count_type[type])
            for feature_name in arr_match_feature_server:
                if feature_name in genome_feature_list:
                    if feature_name == 'GC_CONTENT':
                        genome_dict[feature_name] = data_frame_file[[feature_name]].mean().to_json().split(':')[1][:-1]
                    elif feature_name == 'DNA_LENGTH':
                        genome_dict[feature_name] = data_frame_file[[feature_name]].sum().to_json().split(':')[1][:-1]
            array_genes_features = intersection(arr_match_feature_server, gene_feature_list)
            if len(array_genes_features) != 0:
                array_genes_features.append('GENE_NAME')
                object_genes = create_object_from_data_frame(data_frame_file, array_genes_features,'Gene')
                full_data_features['Gene'] = object_genes
            array_protein_features = intersection(arr_match_feature_server, protein_feature_list)
            if len(array_protein_features) != 0:
                array_protein_features.append('GENE_NAME')
                object_protein = create_object_from_data_frame(data_frame_file,array_protein_features, 'Protein')
                full_data_features['Protein'] = object_protein

            array_general_features = intersection(arr_match_feature_server, general_feature_list)
            if len(array_general_features) != 0:
                if 'GENE_NAME' not in array_general_features:
                    array_general_features.append('GENE_NAME')
                object_general = create_object_from_data_frame(data_frame_file, array_general_features,'General')
                full_data_features['General'] =object_general
            full_data_features['Genome'] = genome_dict
            to_return[fileName] = full_data_features
    return to_return


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def create_object_from_data_frame(data_frame, arr_features_name, type_of_gene):
    if type_of_gene == 'Protein':
        data_frame = data_frame[data_frame["PRODUCT_TYPE"] == 'CDS']

    object_genes = list()
    for index in range(data_frame.shape[0]):
        object_per_gene = {"key":index}
        for feature in arr_features_name:
            # print(feature,type(data_frame[feature].iloc[index,]))
            if type(data_frame[feature].iloc[index,]) != int64 and type(data_frame[feature].iloc[index, ]) != float64:
                object_per_gene[list(features_dict.keys())[list(features_dict.values()).index(feature)]] = data_frame[feature].iloc[index, ]

            elif not np.isnan(data_frame[feature].iloc[index, ]):
                object_per_gene[list(features_dict.keys())[list(features_dict.values()).index(feature)]] = float(
                    data_frame[feature].iloc[index, ])
        object_genes.append(object_per_gene)
    return object_genes


@app.route('/api/gc-content', methods=['GET'])
def gcContent():
    file_name = request.args.getlist('fileName')[0]
    file_name = os.path.splitext(file_name)[0]
    path_to_pickle_files = './BioinformaticsLab/data/data_outputs/features_' + file_name + '.pickle'
    data_frame_file = pd.read_pickle(path_to_pickle_files)
    gc_content = list(data_frame_file[['GC_CONTENT']].values)
    gc_content = [number[0] for number in gc_content]

    return json.dumps(gc_content)


# @app.route('/api/featuresHist', methods=['GET'])
# def numeric_feature_to_hist():
#     file_name = request.args.getlist('fileName')[0]
#     print(file_name)
#     feature_list_by_user = request.args.getlist('featureList[]')
#     print(feature_list_by_user)
#     arr_match_feature_server = match_client_feature_to_df(feature_list_by_user)
#     print("after match", arr_match_feature_server)
#     to_return = {}
#     path_to_pickle_files = './BioinformaticsLab/data/data_outputs/features_' + file_name[:-3] + '.pickle'
#     data_frame_file = pd.read_pickle(path_to_pickle_files)
#     for feature_to_compute in arr_match_feature_server:
#         if feature_to_compute in numeric_feature_list:
#             raw_data = list(data_frame_file[[feature_to_compute]].values)
#             raw_data = [0 if np.isnan(float64(number[0])) else float64(number[0]) for number in raw_data] #TODO: check NaN
#             to_return[feature_to_compute] = raw_data
#     # to_return={'GC_CONTENT':to_return['GC_CONTENT'],'POLAR_AA':to_return['POLAR_AA']}
#     return to_return

# need to be check with postman
@app.route('/api/featuresHist', methods=['GET'])
def numeric_feature_to_hist():
    file_list_names = request.args.getlist('fileList[]')
    feature_list_by_user = request.args.getlist('featureList[]')
    arr_match_feature_server = match_client_feature_to_df(feature_list_by_user)
    print("after match", arr_match_feature_server)
    to_return = {}
    for file_name in file_list_names:
        numeric_of_files = {}
        path_to_pickle_files = './BioinformaticsLab/data/data_outputs/features_' + file_name[:-3] + '.pickle'
        fileName = os.path.splitext(file_name)[0]
        if len(check_existing_files(fileName)) != 2:
            features_on_each_gene(fileName)
        data_frame_file = pd.read_pickle(path_to_pickle_files)
        for feature_to_compute in arr_match_feature_server:
            if feature_to_compute in numeric_feature_list:
                raw_data = list(data_frame_file[[feature_to_compute]].values)
                raw_data = [float64(number[0]) for number in raw_data if not np.isnan(float64(number[0]))] #TODO: check NaN
                # raw_data = [0 if np.isnan(float64(number[0])) else float64(number[0]) for number in raw_data] #TODO: check NaN
                numeric_of_files[feature_to_compute] = raw_data
        to_return[file_name] = numeric_of_files

    return to_return



# file from server
# @app.route('/api/chooseFileFromServer', methods=['POST'])
# def file_from_server():
#     post_data = request.get_json()
#     fileName = post_data.get("fileName")
#     #valid_file = valid(fileName)
#     print("route", fileName)
#     to_return = {
#         'faild_list': [],
#         'status': status_to_client("Uploaded")}
#     return to_return


# download file from s3 server and save it to data_inputs
@app.route('/api/uploadBucketFile', methods=['POST'])
def file_bucket_download():
    post_data = request.get_json()
    fileNameKey = post_data.get("fileName")
    valid_file = valid(fileNameKey)
    if valid_file:
        full_path = os.path.join(path_to_input_file, fileNameKey)
        req = s3_client.meta.client.download_file(aws_storage_bucket_name, fileNameKey, full_path)
        to_return = {
            'faild_list': [],
            'status': status_to_client("Uploaded")}
    else:
        to_return = {
            'faild_list': [],
            'status': status_to_client("")}

    return to_return


@app.route('/api/uploadFile', methods=['POST'])
def file_input_manger():
    to_return = {
        'faild_list': [],
        'status': ''}
    if len(request.files) != 0:
        file = request.files['file']
        status = get_file(file)
        to_return = {
            'faild_list': [],
            'status': status}
    else:
        # return a list of falid acc numbers if one is ok it will save it and others not
        print(request.get_json())
        faild_list, status = download_gb_file_by_id(request.get_json())
        to_return = {
            'faild_list': faild_list,
            'status': status}
    return to_return


def get_file(file):
    # file = request.files['file']
    filename = file.filename
    valid_file = valid(filename)
    if valid_file:

        full_path = os.path.join(path_to_input_file, filename)
        file.save(full_path)
        return status_to_client("Uploaded")
    return status_to_client("")


def valid(filename):
    allowed_file_type = ['gb', 'GB']
    type_of_file = filename.split(".")[1]
    if type_of_file in allowed_file_type:
        return True
    return False


def check_existing_files(filename):
    # method check if the file is already uploaded to the web and check if he got 2 pickle file.
    # fx return a list with the result
    current_input_files = os.listdir(path_to_input_file)
    current_pickle_files = get_current_pickle_files()
    # ---------------remove extension-------------------
    file_name_without_ex = os.path.splitext(filename)[0]
    # current_input_files_without_ex = [os.path.splitext(name)[0] for name in current_input_files]
    get_current_pickle_files_without_ex = [os.path.splitext(name)[0] for name in current_pickle_files]
    # ---------------------------------------------------
    list_of_present = []
    # if file_name_without_ex in current_input_files_without_ex:
    #     list_of_present.append(file_name_without_ex)
    counterOfPickeleFiles = 0
    for name in get_current_pickle_files_without_ex:

        if file_name_without_ex in name:
            counterOfPickeleFiles += 1
            list_of_present.append(name)
            if counterOfPickeleFiles == 2:
                break

    return list_of_present


# @app.route('/api/accession')
def download_gb_file_by_id(acc_id_dict):
    str_acc_id_dict = list(acc_id_dict.values())[0]
    acc_id_list = str_acc_id_dict.split(',')
    Entrez.email = "jcecomputationalbiology@gmail.com"  # Provide an email address
    faild_list_acc = []
    for acc_id in acc_id_list:
        # print("acc", acc_id)
        try:
            with Entrez.efetch(db="nucleotide", id=acc_id, rettype="gbwithparts", retmode="text") as in_handle:
                file_name = f"{acc_id}.gb"  # if going to changed we need to see what to do
                # data = in_handle.read()
                with open((os.path.join(path_to_input_file, file_name)), "w+") as out_handle:
                    out_handle.write(in_handle.read())  # was data
        except:
            faild_list_acc.append(acc_id)

    if not faild_list_acc:
        return [], status_to_client("Uploaded")
    else:
        return faild_list_acc, status_to_client("")


@app.route('/api/chart_data')
def get_chart_data():
    array = list(map(lambda x: {'x': x, 'y': randrange(20)}, range(10)))
    print(jsonify(array))
    return jsonify(array)


def status_to_client(status):
    # need to add more status ?
    # todo : maybe chenged numbers and transfr for words
    dict_status = {
        "Uploaded": {"Uploaded": 200},
        "Failed": {"Failed": 500},
        "Exist": {"Exist": 900},
        "FeatureReceived": {"FeatureReceived": 200}
    }
    return dict_status[status]
    # if status == "Success":
    #     to_return = {"status": 200}
    #     return to_return
    # return {"status": 500}


@app.route('/api/existinglistFiles', methods=['GET'])
def get_existing_files():
    arr = os.listdir(path_to_input_file)
    return jsonify(arr)


def get_current_pickle_files():
    array_of_pickle_files = os.listdir(path_to_pickle_files)
    return array_of_pickle_files


@app.route('/api/listFeatures', methods=['GET'])
def get_features_list():
    dict_to_return = dict()
    dict_to_return['General_Features'] = []
    dict_to_return['Gene_Features'] = []
    dict_to_return['Protein_Features'] = []
    for feature in [*fc.general_features_map]:
        dict_to_return['General_Features'].append(str(feature).partition(".")[-1].replace("_", " "))
    for feature in [*fc.gene_features_map]:
        dict_to_return['Gene_Features'].append(str(feature).partition(".")[-1].replace("_", " "))
    for feature in [*fc.protein_features_map]:
        dict_to_return['Protein_Features'].append(str(feature).partition(".")[-1].replace("_", " "))
    dict_to_return['Genome_Features'] = ['GC CONTENT', 'DNA LENGTH']
    return jsonify(dict_to_return)


@app.route('/api/getFeaturesDescription', methods=['GET'])
def get_features_description():
    return jsonify(features_description)


# need file name for getting the data
@app.route('/api/test', methods=['GET'])
def get_number_of_null_gene_name():
    file_name = request.args.getlist('fileList[]')
    path_to_pickle_file = './BioinformaticsLab/data/data_outputs/features_' + file_name[0][:-3] + '.pickle'
    data_frame_file = pd.read_pickle(path_to_pickle_file)
    df_new = data_frame_file[['GENE_NAME', 'PRODUCT_TYPE']]
    counter_list = {}
    unique_names = df_new['PRODUCT_TYPE'].unique()

    unique_names = list(filter(None, unique_names)) # removed None vaule type

    for index, row in df_new.iterrows():
        if row['GENE_NAME'] == '':
            if row['PRODUCT_TYPE'] in unique_names:
                if row['PRODUCT_TYPE'] in counter_list:
                    counter_list[row['PRODUCT_TYPE']] = counter_list[row['PRODUCT_TYPE']] + 1
                else:
                    counter_list[row['PRODUCT_TYPE']] = 1
    return counter_list
    # cds_counter = 0
    # tRNA_counter = 0
    # rRNA_counter = 0
    # tmRNA_counter = 0
    # for index, row in df_new.iterrows():
    #     if row['GENE_NAME'] == '':
    #         if row['PRODUCT_TYPE'] == 'CDS':
    #             cds_counter = cds_counter + 1
    #         elif row['PRODUCT_TYPE'] == 'tRNA':
    #             tRNA_counter = tRNA_counter + 1
    #         elif row['PRODUCT_TYPE'] == 'rRNA':
    #             rRNA_counter = rRNA_counter + 1
    #         elif row['PRODUCT_TYPE'] == 'tmRNA':
    #             tmRNA_counter = tmRNA_counter + 1
    #
    # to_return ={'Missing CDS': str(cds_counter),
    #             'Missing tRNA': str(tRNA_counter),
    #             'Missing rRNA': str(rRNA_counter),
    #             'Missing tmRNA': str(tmRNA_counter)}
    # return to_return

