import os
import re
import multiprocessing
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Application import ApplicationError
import matplotlib.pyplot as plt

if os.sep == '/':
    picpath = './recombine_pic/'
    libpath = './recombine_lib/'
else:
    picpath = '.\\recombine_pic\\'
    libpath = '.\\recombine_lib\\'

uncertain_symbol='#Uncertain#'

#label->viral group
renameLabel={
    'Alphacoronavirus':                         'AlphaCoV',

    'Porcine_epidemic_diarrhea_virus-swine':    'PEDV-swine',
    'Swine_coronavirus-swine':                  'SwineCoV',
    'Coronavirus_HKU15-swine':                  'HKU15-swine',

    'Human_coronavirus_229E-human':             '229E-human',
    'Human_coronavirus_HKU1-human':             'HKU1-human',
    'Human_coronavirus_NL63-human':             'NL63-human',
    'Human_coronavirus_OC43-human':             'OC43-human',

    'Avian_coronavirus-avian':                  'AvianCoV',

    'SARS-rat':                                 'SARS-rat',
    'SARS-human':                               'SARS-human',
    'SARS-bat':                                 'SARS-bat',
    'SARS-CoV-2':                               'SARS-CoV-2',
    'SARS-CoV-2_Omicron':                       'SARS-CoV-2',
    'SARS-CoV-2_Alpha':                         'SARS-CoV-2',

    'Others_coronavirus-bat':                   'OtherCoV',
    'Others_coronavirus-avian':                 'OtherCoV',
    'Others_coronavirus-dolphin':               'OtherCoV',
    'Human_coronavirus_Others-human':           'OtherCoV',

    'Murine_coronavirus-mouse':                 'MurineCoV',
    'Murine_coronavirus-human':                 'MurineCoV-human',
    'Mouse_coronavirus-mouse':                  'MurineCoV',
    'Mouse_coronavirus-rat':                    'MurineCoV',
    'Rodent_coronavirus-mouse':                 'MurineCoV',

    'MERS-Camel':                               'MERS-camel',
    'MERS-human':                               'MERS-human',
    'MERS-bat':                                 'MERS-bat',

    'Bat_coronavirus-bat':                      'BatCoV',
    'Bat_coronavirus-swine':                    'BatCoV', #only 4 sequences
    'Pangolin_coronavirus-pangolin':            'PangolinCoV',

    'Sarbecovirus-bat':                         'SarbeCoV-bat',
    'NO-HIT':                                   'NO-HIT',
}

#########################################################Recombine_FinalFrag Visualize#########################################################
def plot_color_bars_enhanced(data_row, loci_mode=1 , input_fasta=libpath+'COV-reselect-Cluster7.fasta' ):
    colors = {}
    segments = []
    current_position = 0
    if 'No' in data_row[1]: #
        return
    raw_pic_title=data_row[0].strip()
    pic_title=(raw_pic_title.split('__')[0]+'__'+raw_pic_title.split('__')[1]).replace('_complete_genome','')
    if uncertain_symbol in data_row[1]: #Uncertain result
        pic_title='#Uncertain# '+pic_title
    seq_id_short=(raw_pic_title.strip()).split('.')[0].split(' ')[0].split('__')[0] #seq id
    print('Drawing',seq_id_short)

    rawfilepath=input_fasta
    file_in=open(rawfilepath,'r')
    linenum=0
    for line in file_in.readlines():
        if linenum%2==0:
            seq_id=(line.strip()).split('.')[0].split(' ')[0].split('__')[0][1:]
        else:
            if seq_id_short==seq_id:
                genome_seq=line.strip()
                break
        linenum=linenum+1
    file_in.close()
    whole_length=len(genome_seq)

    #recombine
    for item in data_row[3:]:
        if item=='':
            break
        segment_info = item.split('<-')[0]
        # print(item,segment_info.split('-'))
        segment_length = int(segment_info.split('-')[1][:-2]) - int(segment_info.split('-')[0]) + 1
        label = item.split('<-')[1].split('|')[0]
        source = renameLabel[item.split('<-')[1].split('|')[0]]+'|'+item.split('<-')[1].split('|')[1] #label重命名
        #将NO-HIT|NO-HIT转换为NA|NA
        if source=='NO-HIT|NO-HIT':
            source='NA|NA'
        if source not in colors:
            colors[source] = len(colors)

        #BLAST score
        identity=-1
        if segment_length/whole_length>0.025:
            subseq=genome_seq[int(segment_info.split('-')[0]):(int(segment_info.split('-')[1][:-2])-1)]
            identity=BLAST_identity_score(subseq,label)
            # print(label,identity,len(subseq))
        segments.append((current_position, current_position + segment_length, colors[source], identity))
        current_position += segment_length
    # whole_length=current_position

    #genome
    fig, ax = plt.subplots(figsize=(15, 2))
    for segment in segments:
        ax.axhspan(0.5, 0.75, xmin=segment[0]/whole_length, xmax=segment[1]/whole_length, facecolor=plt.cm.tab10.colors[segment[2]])
        identity=segment[3]
        #BLAST score
        if identity>0:
            text = ax.text((segment[0] + segment[1])/2, 0.82, 'I='+str(identity), weight='bold', ha='center', va='top', color='darkslategray', fontsize=8)#color=colors_orf[i]
            # text.set_path_effects([Stroke(linewidth=0.5, foreground='black'), Normal()])

    #title
    oneline_title_length=125
    if len(pic_title) <= oneline_title_length:
        plt.text(0.00, 0.83+0.15, pic_title, transform=plt.gca().transAxes, verticalalignment='top', horizontalalignment='left',fontsize=12.5)
    else:
        substring = pic_title[:oneline_title_length]
        match = re.search(r"[^A-Za-z0-9]", substring[::-1])
        if match:
            split_index = oneline_title_length - match.start()
            part1 = pic_title[:split_index]
            part2 = pic_title[split_index:]
            if len(part2)<5:
                plt.text(0.00, 0.83+0.15, pic_title, transform=plt.gca().transAxes, verticalalignment='top', horizontalalignment='left',fontsize=12.5)
            else:
                plt.text(0.00, 0.83+0.15+0.15, part1, transform=plt.gca().transAxes, verticalalignment='top', horizontalalignment='left',fontsize=12.5)
                plt.text(0.00, 0.83+0.15, part2, transform=plt.gca().transAxes, verticalalignment='top', horizontalalignment='left',fontsize=12.5)
        else :
            print('No char to split:',pic_title)
            plt.text(0.00, 0.83+0.15, pic_title, transform=plt.gca().transAxes, verticalalignment='top', horizontalalignment='left',fontsize=12.5)

    # legend
    legend_labels = list(colors.keys())
    if "NA|NA" in legend_labels:
        legend_labels.remove("NA|NA")
        legend_labels.append("NA|NA")
    legend_handles = [plt.Rectangle((0,0),1,1, color=plt.cm.tab10.colors[colors[label]]) for label in legend_labels]
    if len(legend_labels)<=6:
        ax.legend(legend_handles, legend_labels, loc='lower center', bbox_to_anchor=(0.5, 0), ncol=3,  frameon=False) #ncol=len(legend_labels)
    elif len(legend_labels)<=9:
        ax.legend(legend_handles, legend_labels, loc='lower center', bbox_to_anchor=(0.5, -0.175), ncol=3,  frameon=False) #ncol=len(legend_labels)
    else:
        ax.legend(legend_handles, legend_labels, loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=3,  frameon=False) #ncol=len(legend_labels)
    #
    ax.set_xlim(0, whole_length)
    ax.set_ylim(0, 1)
    ax.axis('off')
    #
    if loci_mode == 1:
        for i in range(0, whole_length, 2000):
            if i ==0:
                ax.text(i, 0.4, f"{1}bp", ha='center')
            else:
                ax.text(i, 0.4, f"{i}bp", ha='center')
            ax.plot([i, i], [0.475, 0.525], color='black', linewidth=0.75)
    elif loci_mode == 2:
        color_intervals = [int(whole_length * i / len(colors)) for i in range(len(colors) + 1)]
        for interval in color_intervals:
            ax.text(interval, 0.4, f"{interval}bp", ha='center')
    seq_id='EPI_ISL'+data_row[0].strip().split('EPI_ISL')[1].split('|')[0] if ('EPI_ISL' in data_row[0].strip()) else data_row[0].strip().split('.')[0].split(' ')[0].split('__')[0]
    plt.savefig(picpath+seq_id+'.png', dpi=300, format='png')
    # plt.savefig(picpath+seq_id+'.pdf', dpi=300, format='pdf')
    # plt.show()
    plt.close()


#identity
def BLAST_identity_score(sequence,merged_label):
    if merged_label=='NO-HIT':
        return 0
    tmp_blast_fasta=open(libpath+'tmp_blast_seq_pic.fasta','w')
    tmp_blast_fasta.write('>tmp_Seq_for_blast\n'+sequence+'\n')
    tmp_blast_fasta.close()
    blastdb_path=libpath+'parentseq_BLASTDB__'+merged_label #merged_label->blastdb_path
    tmp_blast_result_path=libpath+'tmp_blast_identity.txt'
    num_threads = multiprocessing.cpu_count()
    blastp_cline = NcbiblastnCommandline(query=libpath+'tmp_blast_seq_pic.fasta', db=blastdb_path, outfmt=6, out=tmp_blast_result_path, evalue=1e-5, word_size=7, num_threads=num_threads)
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

    with open(tmp_blast_result_path, "r") as file:
        first_line = file.readline().strip()
        elements = first_line.split("\t")
        if len(elements)<2:
            return 0
        third_element = elements[2]
        identity = round(float(third_element),2)
    return identity



def draw_each_seq(filepath='Recombine_FinalFrag.csv'):
    infile=open(filepath,'r')
    lines=infile.readlines()
    for line in lines[1:]:
        data_row=line.strip().split(',')
        plot_color_bars_enhanced(data_row, loci_mode=1, input_fasta=libpath+'COV-reselect-Cluster7.fasta')


if __name__ == '__main__':
    draw_each_seq()

