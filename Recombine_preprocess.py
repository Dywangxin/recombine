import os
import subprocess
from collections import defaultdict,OrderedDict
import joblib
import warnings
warnings.filterwarnings("ignore")

if os.sep == '/':
    libpath = './recombine_lib/'
else:
    libpath = '.\\recombine_lib\\'

#########################################################Quality Control#########################################################
def quality_control(rawfastafilepath,newfastafilepath,library_strategy,min_seq_length,min_cluster_size):
    print('Quality control start',library_strategy)
    count_dict=defaultdict(int)
    file_obj=open(rawfastafilepath,'r')
    linenum=0
    length_stat_list=[]
    label_list=[]
    for line in file_obj.readlines():
        if linenum%2==0:
            label_id=(line.strip()).split('__')[-1]
            seqid=line.strip()
            label_list.append(label_id)
        else:
            seqlength=len(line.strip())
            length_stat_list.append(seqlength)
            if seqlength>=32000:
                print('Quality WARNING:',seqid,seqlength,'bp long.')
            if seqlength>=min_seq_length:
                count_dict[label_id]=count_dict[label_id]+1
        linenum=linenum+1
    file_obj.close()
    label_list_set=set(label_list)
    raw_dataoutline=(length_stat_list,label_list_set)

    file_obj=open(rawfastafilepath,'r')
    file_obj2=open(newfastafilepath,'w')
    file_obj3=open(newfastafilepath[:-14]+'discarded.fasta','w')
    length_stat_list_filtered=[]
    label_list_filtered=[]
    linenum=0
    for line in file_obj.readlines():
        if linenum%2==0:
            label_id=(line.strip()).split('__')[-1]
            idline=line
        else:
            seqlength=len(line.strip())
            if count_dict[label_id]>=min_cluster_size and seqlength>=min_seq_length:
                file_obj2.write(idline)
                file_obj2.write(line)
                seqlength=len(line.strip())
                length_stat_list_filtered.append(seqlength)
                label_list_filtered.append(label_id)
            else:
                file_obj3.write(idline)
                file_obj3.write(line)
        linenum=linenum+1
    label_list_set_filtered=set(label_list_filtered)
    quality_dataoutline=(length_stat_list_filtered,label_list_set_filtered)
    file_obj.close()
    file_obj2.close()
    file_obj3.close()
    print('Quality control complete',library_strategy)
    return count_dict,raw_dataoutline,quality_dataoutline


#########################################################Fragment Cutting#########################################################
def fragments_from_seq(sequence,fraglength=18):
    sequence=sequence.upper()
    length=len(sequence)
    fragment_list=[sequence[i:i+fraglength] for i in range(0,length-fraglength+1)]
    return fragment_list


def fragment_cutting(library_strategy, filtered_path, fragment_dict_path, rebundantseq_path, all_label_list_path, all_label_dict_path, count_dict, label_cuts_number=-1, fraglength=36):
    print(library_strategy,'Fragment',str(fraglength),'bp generation start.')
    label_all_fragments_set=set([])
    file_obj=open(filtered_path,'r')
    file_obj2=open(rebundantseq_path,'w')
    linenum=0
    test_seq_num=0
    all_label_list=[]
    length_stat_list=[]
    for line in file_obj.readlines():
        if linenum%2==0:
            label_id=(line.strip()).split('__')[-1]
            seq_id=(line.strip()).split('.')[0].split(' ')[0][1:]
            idline=line.strip()
            if label_id not in all_label_list:
                all_label_list.append(label_id)
        else:
            if label_cuts_number_dict[label_id]<label_cuts_number or label_cuts_number==-1:
                fragment_list=fragments_from_seq(line.strip().upper(),fraglength)
                fragment_set=set(fragment_list)

                new_fragment=fragment_set-label_fragments_set_dict[label_id]
                label_fragments_set_dict[label_id]=label_fragments_set_dict[label_id].union(new_fragment)

                for index,frag in enumerate(new_fragment):
                    fragment_label_set_dict[frag]=fragment_label_set_dict[frag].union([label_id])
                label_cuts_number_dict[label_id]+=1
                test_seq_num=0

                seqlength=len(line.strip())
                length_stat_list.append(seqlength)
                if seqlength>=50000:
                    print('Lib WARNING:',seq_id,seqlength,'bp long.')

            if label_cuts_number_dict[label_id]==label_cuts_number or label_cuts_number==-1:
                label_all_fragments_set=label_all_fragments_set.union(label_fragments_set_dict[label_id])

            if label_cuts_number_dict[label_id]==label_cuts_number and label_cuts_number!=-1: #and test_seq_num<test_Seq_generate_num
                # file_obj2.write(idline+'\n')
                # file_obj2.write(line.strip()+'\n')
                test_seq_num+=1
        linenum=linenum+1
    file_obj.close()
    file_obj2.close()

    joblib.dump(fragment_label_set_dict,fragment_dict_path)
    print(library_strategy,'Fragment',str(fraglength),'bp generation complete.')

    library_dataoutline=(length_stat_list,set(all_label_list))

    all_label_list.append('Unknown_fragments')
    all_label_dict=dict()
    for index,label in enumerate(all_label_list):
        all_label_dict[label]=index
    joblib.dump(all_label_list,all_label_list_path)
    joblib.dump(all_label_dict,all_label_dict_path)

    # file_obj4=open('Recombine_all_label_list'+'.csv','w')
    # for index,label in enumerate(all_label_list):
    #     file_obj4.write(str(index)+','+label+','+str(count_dict[label])+'\n')
    # file_obj4.close()

    return library_dataoutline


#########################################################BLAST DB Building#########################################################
def build_one_blast_db(fasta_path, blastdb_path):
    cmd = f"makeblastdb -in {fasta_path} -dbtype nucl -out {blastdb_path}"
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        raise Exception("Error occurred while building BLAST database "+blastdb_path)


def extract_fasta_for_blast(whole_fasta_path,merged_label):
    file_in=open(whole_fasta_path,'r')
    outpath=libpath+'tmp_fasta_for_blast.fasta'
    file_out=open(outpath,'w')
    linenum=0
    for line in file_in.readlines():
        if linenum%2==0:
            line1=line.strip()
            label_id=(line.strip()).split('__')[-1]
        else:
            line2=line.strip()
            if merged_label in label_id:
                file_out.write(line1+'\n')
                file_out.write(line2+'\n')
        linenum=linenum+1
    file_in.close()
    file_out.close()
    return outpath

#BLAST lib
def build_blast_db_caller(whole_fasta_path1,whole_fasta_path2):
    all_label_list=joblib.load(libpath+'lib_all_label_list')
    merged_all_label_list=[]
    for label in all_label_list:
        merged_l=label.split('-Cluster')[0]
        if merged_l not in merged_all_label_list:
            merged_all_label_list.append(merged_l)

    for merged_l in merged_all_label_list:
        if merged_l=='Unknown_fragments':
            continue
        tmp_fasta_path=extract_fasta_for_blast(whole_fasta_path1,merged_l)
        blastdb_path=libpath+'parentseq_BLASTDB__'+merged_l #merged_label->blastdb_path
        build_one_blast_db(tmp_fasta_path, blastdb_path)

    blastdb_path=libpath+'parentseq_BLASTDB__'+'All' #
    build_one_blast_db(whole_fasta_path1, blastdb_path)

    blastdb_path=libpath+'parentseq_BLASTDB__'+'Raw' #
    build_one_blast_db(whole_fasta_path2, blastdb_path)

    os.remove(tmp_fasta_path)
    print('Building BLASTDB complete')


#########################################################Preprocessing Main#########################################################
#[WARNING]For all .fasta files used in this project, the sequence must be in one line rather than in several lines, otherwise ERROR.
rawfilepath=libpath+'COV-reselect-Cluster7.fasta'


label_fragments_set_dict=defaultdict(set)
label_cuts_number_dict=defaultdict(int)
fragment_label_set_dict=defaultdict(set)

#parameters
standard_fraglength=36         #standard fragments 36 bp
second_level_fralength=9        #secondary fragments 9 bp
global_label_cuts_number=100    #recommended 100; -1 includes all filtered seq


def fragment_library_main(library_strategy=''):
    #parameters
    min_seq_length=3500 #recommended
    min_cluster_size=10 #recommended
    filtered_path=rawfilepath[:-6]+'filtered.fasta'
    print('Fragment library strategy: min_seq_length >=',min_seq_length,'bp; min_cluster_size >=',min_cluster_size)

    #quality control
    count_dict,raw_dataoutline,quality_dataoutline=quality_control(rawfilepath,filtered_path,library_strategy,min_seq_length=min_seq_length,min_cluster_size=min_cluster_size)

    #1st fragment cutting
    label_fragments_set_dict.clear()
    label_cuts_number_dict.clear()
    fragment_label_set_dict.clear()
    rebundantseq_path=libpath+'lib_rebundant_seq.fasta'
    all_label_list_path=libpath+'lib_all_label_list'
    all_label_dict_path=libpath+'lib_all_label_dict'
    fragment_dict_path=libpath+'lib_dict1'
    library_dataoutline=fragment_cutting(library_strategy, filtered_path, fragment_dict_path, rebundantseq_path, all_label_list_path, all_label_dict_path, count_dict, label_cuts_number=global_label_cuts_number, fraglength=standard_fraglength)

    #2nd fragment cutting
    label_fragments_set_dict.clear()
    label_cuts_number_dict.clear()
    fragment_label_set_dict.clear()
    fragment_dict_path2=libpath+'lib_dict2'
    fragment_cutting(library_strategy, filtered_path, fragment_dict_path2, rebundantseq_path, all_label_list_path, all_label_dict_path, count_dict, label_cuts_number=global_label_cuts_number, fraglength=second_level_fralength)

    #BLAST lib
    build_blast_db_caller(whole_fasta_path1=filtered_path,whole_fasta_path2=rawfilepath)


if __name__ == '__main__':
    fragment_library_main()

