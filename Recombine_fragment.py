import joblib
import warnings
import os
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Application import ApplicationError
from Bio import pairwise2
import multiprocessing
import pandas as pd

warnings.filterwarnings("ignore")

if os.sep == '/':
    libpath = './recombine_lib/'
else:
    libpath = '.\\recombine_lib\\'

uncertain_symbol='#Uncertain#'

#########################################################Fragment Vote Util#########################################################
#prior knowledge
def BLAST_Prior(sequence):
    tmp_blast_fasta=open(libpath+'tmp_blast_seq_prior.fasta','w')
    tmp_blast_fasta.write('>tmp_Seq_for_blast\n'+sequence+'\n')
    tmp_blast_fasta.close()
    blastdb_path=libpath+'parentseq_BLASTDB__All' #blastdb_path
    tmp_blast_result_path=libpath+'tmp_blast_identity_prior.txt'
    num_threads = multiprocessing.cpu_count()
    blastp_cline = NcbiblastnCommandline(query=libpath+'tmp_blast_seq_prior.fasta', db=blastdb_path, outfmt=6, out=tmp_blast_result_path, evalue=1e-5, word_size=7, num_threads=num_threads)
    try:
        stdout, stderr = blastp_cline()
    except ApplicationError as err:
        print(err.cmd)
        print(err.returncode)
        print(err.stderr)
        print(err.stdout)
        raise
    with open(tmp_blast_result_path, "r") as file:
        first_line = file.readline().strip()
        elements = first_line.split("\t")
        if len(elements)<2:
            return -1,0,0
        third_element = elements[2]
        identity = float(third_element)
        parent_idline=elements[1].strip()
        wholeseq_flag=(parent_idline.strip()).split('__')[-1]
        high_identity = 0 if identity<80 else 1  #Uncertain#：primitive BLAST coverage<80，uncertain result
        return wholeseq_flag,high_identity,identity

#cut
def fragments_from_seq(sequence,fraglength):
    sequence=sequence.upper()
    length=len(sequence)
    fragment_list=[sequence[i:i+fraglength] for i in range(0,length-fraglength+1)] #
    return fragment_list

def find_segments(seqlabels):
    segments = []
    start, end = None, None
    for i, x in enumerate(seqlabels):
        if x == -1:
            if start is None:
                start,end = i,i
            else:
                end = i
        elif start is not None:
            segments.append((start, end))
            start, end = None, None
    if start is not None and end is not None:
        segments.append((start, end))
    return segments


#########################################################Fragment Vote Main#########################################################
#Fragment Vote
def fragment_vote(sequence,fraglength_for_context,fragment_label_set_dict1,fragment_label_set_dict2,all_label_dict,standard_lib=1,fraglength=36,wholeseq_flag=-1):
    #cut
    fragment_list=fragments_from_seq(sequence,fraglength=fraglength)
    #set lib
    fragment_label_set_dict=fragment_label_set_dict1 if standard_lib==1 else fragment_label_set_dict2

    #vote
    vote_pool=[]
    for frag in fragment_list:
        vote_pool.extend(list(fragment_label_set_dict[frag]))

    if len(vote_pool) == 0 and standard_lib==1:
        unknown_seqlabels=fragment_vote(sequence,fraglength_for_context,fragment_label_set_dict1,fragment_label_set_dict2,all_label_dict,standard_lib=2,fraglength=second_level_fralength,wholeseq_flag=wholeseq_flag)
        return unknown_seqlabels

    elif len(vote_pool) == 0:
        # print('Vote for Unknown_fragments',len(sequence),'bp')
        return [all_label_dict['Unknown_fragments']]*len(sequence)

    max_label_tmp=max(set(vote_pool),key=vote_pool.count)
    maxlabel_count=vote_pool.count(max_label_tmp)
    max_label=-1
    if wholeseq_flag!=-1 and vote_pool.count(wholeseq_flag)==maxlabel_count:
        max_label=wholeseq_flag
    else:
        max_label=max(set(vote_pool),key=vote_pool.count)

    seqlabels=[-1]*(len(sequence))
    color_index_list=[index for index,frag in enumerate(fragment_list) if max_label in fragment_label_set_dict[frag]]
    for index in color_index_list:
        seqlabels[index:index+fraglength]=[all_label_dict[max_label]]*(fraglength)
    # print(seqlabels)

    if wholeseq_flag==-1:
        wholeseq_flag=max_label

    unknown_fragments=find_segments(seqlabels)
    for unknown in unknown_fragments:
        start,end=unknown[0],unknown[1]
        if end-start+1<=(fraglength_for_context):
            if start!=0:
                seqlabels[start:end+1]=[seqlabels[start-1]]*(end-start+1)
            elif start==0 and end!=(len(seqlabels)-1):
                seqlabels[start:end+1]=[seqlabels[end+1]]*(end-start+1)
            elif start==0 and end==(len(seqlabels)-1):
                pass
            continue
        else:
            unknown_seqlabels=fragment_vote(sequence[start:end+1],fraglength_for_context,fragment_label_set_dict1,fragment_label_set_dict2,all_label_dict,standard_lib=1,fraglength=fraglength,wholeseq_flag=wholeseq_flag)
            seqlabels[start:end+1]=unknown_seqlabels
    return seqlabels

#after vote
def smooth_tiny_label(seqlabels,all_label_list,fraglength_for_context):
    label_list,merged_label_list=[],[]
    max_labelid=max(set(seqlabels),key=seqlabels.count)
    maxlabel=all_label_list[max_labelid]

    start=0
    for index, label in enumerate(seqlabels):
        if index==(len(seqlabels)-1):
            label_list.append((start,index,seqlabels[start]))
        elif label==seqlabels[start]:
            continue
        else:
            label_list.append((start,index-1,seqlabels[start]))
            start=index
    if label_list[0][1]==(len(seqlabels)-1):
        pass
    else:
        for tup in label_list:
            if all_label_list[tup[2]]!=maxlabel and   (tup[1]-tup[0]+1)<=fraglength_for_context:
                start,end=tup[0],tup[1]
                if start!=0:
                    seqlabels[start:end+1]=[seqlabels[start-1]]*(end-start+1)
                elif start==0 and end!=(len(seqlabels)-1):
                    seqlabels[start:end+1]=[seqlabels[end+1]]*(end-start+1)
    return seqlabels

#BLAST for check fragment
def BLAST_for_check_fragment(tmp_blast_fasta_path,filepath_id,merged_label):
    if merged_label=='Unknown_fragments':
        # merged_label='All'
        return 0
    blastdb_path=libpath+'parentseq_BLASTDB__'+merged_label
    tmp_blast_result_path=libpath+'tmp_blast_results_check.txt'
    num_threads = multiprocessing.cpu_count()
    blastp_cline = NcbiblastnCommandline(query=tmp_blast_fasta_path, db=blastdb_path, outfmt=6, out=tmp_blast_result_path, evalue=1e-5, word_size=7, num_threads=num_threads)
    try:
        stdout, stderr = blastp_cline()
    except ApplicationError as err:
        print(err.cmd)
        print(err.returncode)
        print(err.stderr)
        print(err.stdout)
        raise
    with open(tmp_blast_result_path, "r") as file:
        first_line = file.readline().strip()
        elements = first_line.split("\t")
        if len(elements)<2:
            return 0
        third_element = elements[2]
        identity = float(third_element)
    return identity

#seq_id_short_for_object
def check_fragments_by_BLAST(seq_id_short,seq,merged_seqlabels,merged_all_label_list,fraglength_for_context):
    merged_label_list=[]
    start=0
    for index, label in enumerate(merged_seqlabels):
        if index==(len(merged_seqlabels)-1):
            merged_label_list.append((start,index,merged_seqlabels[start]))
        elif label==merged_seqlabels[start]:
            continue
        else:
            merged_label_list.append((start,index-1,merged_seqlabels[start]))
            start=index
    if len(merged_label_list)==1:
        return merged_seqlabels

    for start, end, label_index in merged_label_list:
        fraglength=end-start+1
        if fraglength>fraglength_for_context*2:
            continue
        merged_label=merged_all_label_list[label_index]
        subseq=seq[start:end+1]
        filepath_id=start
        tmp_blast_fasta_path=libpath+'tmp_blast_seq_checkBLAST.fasta'
        tmp_blast_fasta=open(tmp_blast_fasta_path,'w')
        tmp_blast_fasta.write('>tmp_Seq_for_blast\n'+subseq+'\n')
        tmp_blast_fasta.close()

        identity=BLAST_for_check_fragment(tmp_blast_fasta_path,filepath_id,merged_label)

        if identity < 95:
            if start!=0:
                merged_seqlabels[start:end+1]=[merged_seqlabels[start-1]]*(end-start+1)
            elif start==0 and end!=(len(merged_seqlabels)-1):
                merged_seqlabels[start:end+1]=[merged_seqlabels[end+1]]*(end-start+1)

    return merged_seqlabels


#########################################################Fragment Output#########################################################
def merge_seqlabels(seqlabels,all_label_list,merged_all_label_dict):
    merged_seqlabels=[0]*len(seqlabels)
    for index,l in enumerate(seqlabels):
        merged_seqlabels[index]=merged_all_label_dict[all_label_list[l].split('-Cluster')[0]]
    return merged_seqlabels

#
def output_recombine(seq_id,seqlabels,merged_seqlabels,fragment_final_list,all_label_list,merged_all_label_list,high_identity):
    label_list,merged_label_list=[],[]
    max_labelid=max(set(seqlabels),key=seqlabels.count)
    maxlabel=all_label_list[max_labelid]
    max_percentage=int(seqlabels.count(max_labelid)/len(seqlabels)*100)

    start=0
    for index, label in enumerate(seqlabels):
        if index==(len(seqlabels)-1):
            label_list.append((start,index,seqlabels[start]))
        elif label==seqlabels[start]:
            continue
        else:
            label_list.append((start,index-1,seqlabels[start]))
            start=index

    max_labelid_merged=max(set(merged_seqlabels),key=merged_seqlabels.count)
    maxlabel_merged=merged_all_label_list[max_labelid_merged]
    max_percentage_merged=int(merged_seqlabels.count(max_labelid_merged)/len(merged_seqlabels)*100)

    start=0
    for index, label in enumerate(merged_seqlabels):
        if index==(len(merged_seqlabels)-1):
            merged_label_list.append((start,index,merged_seqlabels[start]))
        elif label==merged_seqlabels[start]:
            continue
        else:
            merged_label_list.append((start,index-1,merged_seqlabels[start]))
            start=index

    if merged_label_list[0][1]==(len(merged_seqlabels)-1):
        lines_merged=seq_id+',No,'+merged_all_label_list[merged_label_list[0][2]]+'='+str(max_percentage_merged)+'%'+'\n'
    else:
        lines_merged=seq_id+',Yes,'+maxlabel_merged+'='+str(max_percentage)+'%'
        for tup in merged_label_list:
            if merged_all_label_list[tup[2]]!=maxlabel_merged:
                lines_merged=lines_merged+(','+str(tup[0]+1)+'-'+str(tup[1]+1)+'bp='+str(tup[1]-tup[0])+'bp '+ merged_all_label_list[tup[2]] )

        lines_merged=lines_merged+'\n'

    if len(fragment_final_list)==0:
        pass
    else:
        dict_seq_fraglength=dict()
        for tup in fragment_final_list:
            fragment_loci,parent_seqid,parent_loci=tup[0],tup[1][0],tup[2]
            fragment_length=fragment_loci[1]-fragment_loci[0]
            if parent_seqid in dict_seq_fraglength:
                dict_seq_fraglength[parent_seqid] += fragment_length
            else:
                dict_seq_fraglength[parent_seqid] = fragment_length
        max_parent_seq = max(dict_seq_fraglength, key=dict_seq_fraglength.get)
        max_parent_seq_percentage = int(dict_seq_fraglength[max_parent_seq]/len(seqlabels)*100)

    if len(fragment_final_list)==0: #no recombination
        if high_identity==1:
            lines_final=seq_id+',No,'+merged_all_label_list[merged_label_list[0][2]]+'|~='+str(max_percentage_merged)+'%'+'\n'
        else: #uncertain
            lines_final=seq_id+','+uncertain_symbol+'No,'+merged_all_label_list[merged_label_list[0][2]]+'|~='+str(max_percentage_merged)+'%'+'\n'
    else: #recombination
        if high_identity==1:
            lines_final=seq_id+',Yes,'+maxlabel_merged+'|'+max_parent_seq+'='+str(max_parent_seq_percentage)+'%'
        else: #uncertain
            lines_final=seq_id+','+uncertain_symbol+'Yes,'+maxlabel_merged+'|'+max_parent_seq+'='+str(max_parent_seq_percentage)+'%'
        for tup in fragment_final_list:
            #fragment_final_list tuple=(fragment_loci=(fragment_start,fragment_end), parent_seq=(parent id, parent merged label), parent_loci=(parent_start,parent_end) )
            fragment_loci,parent_seq,parent_loci=tup[0],tup[1],tup[2]
            if parent_loci[0]=='None':
                # print(fragment_loci,parent_seq,parent_loci)
                lines_final=lines_final+(','+str(fragment_loci[0]+1)+'-'+str(fragment_loci[1]+1)+'bp<-'+parent_seq[1]+'|'+parent_seq[0]+'|'+'None'+'-'+'None'+'bp')
            else: #high_identity=1
                lines_final=lines_final+(','+str(fragment_loci[0]+1)+'-'+str(fragment_loci[1]+1)+'bp<-'+parent_seq[1]+'|'+parent_seq[0]+'|'+str(parent_loci[0]+1)+'-'+str(parent_loci[1]+1)+'bp')
        lines_final=lines_final+'\n'

    return lines_merged,lines_final


#########################################################Parents Detection & Positioning#########################################################
#positioning
def find_similar_subseq(subseq, parentseq):
    alignments = pairwise2.align.localms(parentseq, subseq, 1, -1, -1, -0.1) #pairwise2.align.localms > pairwise2.align.localxx
    best_alignment = max(alignments, key=lambda x: x.score)
    seqB=best_alignment.seqB
    start_pos=seqB.find(seqB.strip('-')) #
    return start_pos, start_pos+len(subseq)-1

#1st BLAST->2nd BLAST->3rd BLAST
def BLAST_find_parent(seq_id,merged_label):
    if merged_label=='Unknown_fragments':
        merged_label='All'
    blastdb_path=libpath+'parentseq_BLASTDB__'+merged_label #merged_label->blastdb_path
    tmp_blast_result_path=libpath+'tmp_blast_results.xml'
    num_threads = multiprocessing.cpu_count()
    blastp_cline = NcbiblastnCommandline(query=libpath+'tmp_blast_seq.fasta', db=blastdb_path, outfmt=5, out=tmp_blast_result_path, evalue=1e-5, word_size=7, num_threads=num_threads)
    # cmd = f"blastn -query {libpath+'tmp_blast_seq.fasta'} -db {blastdb_path} -out {tmp_blast_result_path} -outfmt 5"
    # blastp_cline()
    try:
        stdout, stderr = blastp_cline()
    except ApplicationError as err:
        print(err.cmd)
        print(err.returncode)
        print(err.stderr)
        print(err.stdout)
        raise
    result_handle = open(tmp_blast_result_path)
    blast_records = NCBIXML.parse(result_handle)
    top_hit = None
    top_score = 0
    for record in blast_records:
        for alignment in record.alignments:
            current_hit_idline=alignment.title.split(' ')[1]
            for hsp in alignment.hsps:
                if (hsp.score > top_score) and (seq_id not in current_hit_idline):
                # if hsp.score > top_score:
                    top_score = hsp.score
                    top_hit = alignment.title
    if top_hit==None:
        if merged_label=='Raw':   #3rd BLAST fail->fail
            return None
        elif merged_label=='All': #2nd BLAST fail->3rd BLAST (Raw)
            top_hit = BLAST_find_parent(seq_id,'Raw')
        else:                     #1st BLAST fail->2nd BLAST (All)
            top_hit = BLAST_find_parent(seq_id,'All')
    return top_hit

#1st BLAST[spatio-temporal condition]
def BLAST_find_parent_spacetime(seq_id,merged_label,seq_region,seq_year,seqs_info_dataframe):
    if merged_label=='Unknown_fragments':
        merged_label='All'
    #BLAST
    blastdb_path=libpath+'parentseq_BLASTDB__'+merged_label #merged_label->blastdb_path
    tmp_blast_result_path=libpath+'tmp_blast_results.xml'
    num_threads = multiprocessing.cpu_count()
    blastp_cline = NcbiblastnCommandline(query=libpath+'tmp_blast_seq.fasta', db=blastdb_path, outfmt=5, out=tmp_blast_result_path, evalue=1e-5, word_size=7, num_threads=num_threads)
    # cmd = f"blastn -query {libpath+'tmp_blast_seq.fasta'} -db {blastdb_path} -out {tmp_blast_result_path} -outfmt 5"
    # blastp_cline()
    try:
        stdout, stderr = blastp_cline()
    except ApplicationError as err:
        print(err.cmd)
        print(err.returncode)
        print(err.stderr)
        print(err.stdout)
        raise
    result_handle = open(tmp_blast_result_path)
    blast_records = NCBIXML.parse(result_handle)
    top_hit = None
    top_score = 0
    # [spatio-temporal condition]
    for record in blast_records:
        for alignment in record.alignments:
            current_hit_idline=alignment.title.split(' ')[1]
            #spatio-temporal information
            current_hit_seq_id=(current_hit_idline.strip()).split('.')[0].split(' ')[0].split('__')[0]
            current_hit_region = seqs_info_dataframe.loc[current_hit_seq_id, 'Region']
            current_hit_year = seqs_info_dataframe.loc[current_hit_seq_id, 'Year']
            # [spatio-temporal condition] same location & within 10 years
            if seq_region=='Unknown' or seq_year=='Unknown':
                break
            if current_hit_region=='Unknown' or current_hit_year=='Unknown' or seq_region=='Unknown' or seq_year=='Unknown' or current_hit_region!=seq_region or not ( abs(int(seq_year)-int(current_hit_year))<=10 ):
                continue
            for hsp in alignment.hsps:
                if (hsp.score > top_score) and (seq_id not in current_hit_idline):
                    # if hsp.score > top_score:
                    top_score = hsp.score
                    top_hit = alignment.title
    # [temporal condition]
    if top_hit == None:
        for record in blast_records:
            for alignment in record.alignments:
                current_hit_idline=alignment.title.split(' ')[1]
                #解析时空信息
                current_hit_seq_id=(current_hit_idline.strip()).split('.')[0].split(' ')[0].split('__')[0]
                current_hit_year = seqs_info_dataframe.loc[current_hit_seq_id, 'Year']
                # [temporal condition] within 10 years
                if seq_year=='Unknown':
                    break
                if current_hit_year=='Unknown' or seq_year=='Unknown' or not ( ( abs(int(seq_year)-int(current_hit_year))<=10 ) ):
                    continue
                for hsp in alignment.hsps:
                    if (hsp.score > top_score) and (seq_id not in current_hit_idline):
                        # if hsp.score > top_score:
                        top_score = hsp.score
                        top_hit = alignment.title
    # [spatio condition] same location
    if top_hit == None:
        for record in blast_records:
            for alignment in record.alignments:
                current_hit_idline=alignment.title.split(' ')[1]
                current_hit_seq_id=(current_hit_idline.strip()).split('.')[0].split(' ')[0].split('__')[0]
                current_hit_region = seqs_info_dataframe.loc[current_hit_seq_id, 'Region']
                # print(current_hit_seq_id,current_hit_region,current_hit_year)
                # [spatio condition] same location
                if seq_region=='Unknown':
                    break
                if current_hit_region=='Unknown' or seq_region=='Unknown' or current_hit_region!=seq_region:
                    continue
                for hsp in alignment.hsps:
                    if (hsp.score > top_score) and (seq_id not in current_hit_idline):
                        # if hsp.score > top_score:
                        top_score = hsp.score
                        top_hit = alignment.title
    #no 2nd/3rd BLAST with [spatio-temporal condition]
    if top_hit == None:
        top_hit = BLAST_find_parent(seq_id,merged_label)
    return top_hit


# parent_seq=merged_label_to_parent_seq_dict[merged_all_label_list[label_index]]
def multi_Aligning(parent_seq,seq,start,end,parent_merged_label):
    if parent_seq==None:
        fragment_final=((start,end),('NO-HIT','NO-HIT'),('None','None'))
        return fragment_final
    parent_start_index, parent_end_index=find_similar_subseq(parentseq=parent_seq[2], subseq=seq[start:end+1])
    #[NOTE]parent_merged_label may have been changed due to the 3rd BLAST
    parent_real_merged_label=((parent_seq[1]).split('__')[-1]).split('-Cluster')[0]
    fragment_final=((start,end),(parent_seq[0],parent_real_merged_label),(parent_start_index,parent_end_index))
    return fragment_final

#parent seq & fragment positioning
def find_parent_seqs(seq_id_long,seq_id_short_for_object,seq,merged_seqlabels,merged_all_label_list):
    merged_label_list=[]
    start=0
    for index, label in enumerate(merged_seqlabels):
        if index==(len(merged_seqlabels)-1):
            merged_label_list.append((start,index,merged_seqlabels[start]))
        elif label==merged_seqlabels[start]:
            continue
        else:
            merged_label_list.append((start,index-1,merged_seqlabels[start]))
            start=index
    if len(merged_label_list)==1: #no recombination
        return []

    #mergelabel_fragloc_dict:label->loc_list
    mergelabel_fragloc_dict = {}
    for start, end, label_index in merged_label_list:
        merged_label=merged_all_label_list[label_index]
        if merged_label in mergelabel_fragloc_dict:
            mergelabel_fragloc_dict[merged_label].append((start, end))
        else:
            mergelabel_fragloc_dict[merged_label] = [(start, end)]

    merged_labels_for_this_seq=mergelabel_fragloc_dict.keys()
    merged_label_to_parent_seq_dict={} #dict: merged label-> parent seq tuple
    for merged_label in merged_labels_for_this_seq:
        fraglocs=mergelabel_fragloc_dict[merged_label]
        subseq_fraglist=[seq[tup[0]:tup[1]+1] for tup in fraglocs]
        subseq=''.join(subseq_fraglist)
        tmp_blast_fasta=open(libpath+'tmp_blast_seq.fasta','w')
        tmp_blast_fasta.write('>tmp_Seq_for_blast\n'+subseq+'\n')
        tmp_blast_fasta.close()

        #parent search
        if spacetime_BLAST_flag==True:#[spatio-temporal condition]
            # df = pd.read_csv('Recombine_spacetime_DICT.csv')
            # location_region_dict = df.set_index('Location')['Region'].to_dict()
            seqs_info_dataframe=pd.read_csv('Recombine_spacetime_info_raw.csv', index_col='ID')
            #[NOTE]The .csv file (Recombine_spacetime_info_raw.csv) includes the spatio-temporal information of filtered sequences for local library construction.
            seq_region = seqs_info_dataframe.loc[seq_id_short_for_object, 'Region']
            seq_year = seqs_info_dataframe.loc[seq_id_short_for_object, 'Year']
            top_hit = BLAST_find_parent_spacetime(seq_id_short_for_object,merged_label,seq_region,seq_year,seqs_info_dataframe)
        if spacetime_BLAST_flag==False: #[no spatio-temporal condition]
            top_hit = BLAST_find_parent(seq_id_short_for_object,merged_label)

        #3rd BLAST fail
        if top_hit==None:
            merged_label_to_parent_seq_dict[merged_label]=None
            continue
        parent_idline=top_hit.split(' ')[1]

        #parent fasta
        file_searchseqs=open(rawfilepath[:-6]+'.fasta','r')
        linenum=0
        for line in file_searchseqs.readlines():
            if linenum%2==0:
                id_line=line.strip()
                seq_id=(line.strip()).split('.')[0].split(' ')[0].split('__')[0][1:]
            else:
                if seq_id in parent_idline:
                    seq_line=line.strip().upper()
                    break
            linenum+=1
        file_searchseqs.close()
        parent_seq=(seq_id,parent_idline,seq_line)
        merged_label_to_parent_seq_dict[merged_label]=parent_seq

    #fragment positioning
    #multiprocess
    parent_seq_args=[merged_label_to_parent_seq_dict[merged_all_label_list[label_index]] for start,end,label_index in merged_label_list]
    seq_args=[seq for s in merged_label_list]
    start_args=[start for start,end,label_index in merged_label_list]
    end_args=[end for start,end,label_index in merged_label_list]
    merged_label_args=[merged_all_label_list[label_index] for start,end,label_index in merged_label_list]
    args = zip(parent_seq_args, seq_args, start_args, end_args, merged_label_args)
    # print('Multi-Core:',multiprocessing.cpu_count())
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    fragment_final_list_unsorted = pool.starmap(multi_Aligning, args) #multi_Aligning(parent_seq,seq,start,end)
    pool.close()
    pool.join()
    # sort frags by loci
    fragment_final_list = sorted(fragment_final_list_unsorted, key=lambda x: x[0][0])
    return fragment_final_list


#########################################################Fragment Main#########################################################
def predict_and_evaluation(running_mode, inpath, outpath_basis):
    print('Running Mode',running_mode,'start')
    #load
    all_label_list=joblib.load(libpath+'lib_all_label_list')
    all_label_dict=joblib.load(libpath+'lib_all_label_dict')
    fragment_label_set_dict1=joblib.load(libpath+'lib_dict1')
    fragment_label_set_dict2=joblib.load(libpath+'lib_dict2')
    merged_all_label_list=[]
    for label in all_label_list:
        merged_l=label.split('-Cluster')[0]
        if merged_l not in merged_all_label_list:
            merged_all_label_list.append(merged_l)
    merged_all_label_dict=dict()
    for index,label in enumerate(merged_all_label_list):
        merged_all_label_dict[label]=index

    file_obj1=open(inpath,'r')
    file_obj4=open(outpath_basis+'_FinalFrag.csv','w')  ##standard fragment result, showing genetic origins, for interspecies recombination only
    file_obj4.write('Seq,Recombination,Main,GeneticFragments\n')
    linenum=0
    for line in file_obj1.readlines():
        if linenum%2==0:
            seq_id=line.strip()[1:]
            seq_id_short=(line.strip()).split('.')[0].split(' ')[0].split('__')[0][1:]
            label_id=(line.strip()).split('__')[-1]
        else:
            seq=line.strip().upper()

            if len(seq)<15000: #Too short sequence shall not be analyzed
                # print('Too short seq:',seq_id_short)
                linenum+=1
                continue

            #BLAST
            wholeseq_flag,high_identity,identity=BLAST_Prior(sequence=seq)
            print('Doing',seq_id_short)

            #vote
            seqlabels=fragment_vote(seq,fraglength_for_context,fragment_label_set_dict1,fragment_label_set_dict2,all_label_dict,standard_lib=1,fraglength=standard_fraglength,wholeseq_flag=wholeseq_flag)

            #tiny fragments
            seqlabels=smooth_tiny_label(seqlabels,all_label_list,fraglength_for_context)

            #recombination detection
            merged_seqlabels=merge_seqlabels(seqlabels,all_label_list,merged_all_label_dict)

            if running_mode==1:
                fragment_final_list=[]
            elif running_mode==2:
                merged_seqlabels=check_fragments_by_BLAST(seq_id_short,seq,merged_seqlabels,merged_all_label_list,fraglength_for_context)
                fragment_final_list=find_parent_seqs(seq_id,seq_id_short,seq,merged_seqlabels,merged_all_label_list)
            else : #running_mode==3
                fragment_final_list=find_parent_seqs(seq_id,seq_id_short,seq,merged_seqlabels,merged_all_label_list)

            lines2,lines4=output_recombine(seq_id,seqlabels,merged_seqlabels,fragment_final_list,all_label_list,merged_all_label_list,high_identity)
            if running_mode==1:
                file_obj4.write(lines2)
            else:
                file_obj4.write(lines4)

        linenum+=1
    file_obj1.close()
    file_obj4.close()
    print('Recombine detection complete')


#########################################################Caller Main#########################################################
rawfilepath=libpath+'COV-reselect-Cluster7.fasta' #fasta for lib
standard_fraglength=36
second_level_fralength=9
fraglength_for_context=24

inpath=libpath+'COV-reselect-Cluster7.fasta' #input file for recombination detection

#Fast Detection:                    RunningMode=1,  vote only
#Regular Detection[recommended]:    RunningMode=2,  vote + fragment_check + parental&fragments
#Precise Detection:                 RunningMode=3,  vote + parental&fragments
running_mode=2
# spacetime_BLAST_flag=False ##Enable spatio-temporal conditions
spacetime_BLAST_flag=False ##Disable spatio-temporal conditions
if __name__ == '__main__':
    predict_and_evaluation(running_mode, inpath, outpath_basis='Recombine')
