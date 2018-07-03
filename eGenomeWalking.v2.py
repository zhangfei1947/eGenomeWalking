#coding=utf-8
#by zhangfei@novogene.com

import os,argparse,sys,re,glob

parser = argparse.ArgumentParser(description="e-GenomeWalking")
parser.add_argument('--cdsfasta',help='Initiation_only_unigene.blast.cds.fasta or Termination_only_unigene.blast.cds.fasta',required=True)
parser.add_argument('--cleandata',help='clean data,split by ",", *_1.clean.fq or *_1.clean.fq.gz is OK',required=True)
parser.add_argument('--terminal',help='which Terminal to prolong',choices=['5','3'],default='3')
argv=vars(parser.parse_args())
cdsfasta=argv['cdsfasta'].strip()
cleandata_list=argv['cleandata'].strip().split(',')
terminal=argv['terminal'].strip()

root = os.getcwd()

def prepare_database():
        allfq = root+'/DB.fq'
        allfa = root+'/DB.fa'
        for each in cleandata_list:
                if each.endswith('.gz'):
                        os.system('zcat {each} >> {allfq}'.format(each=each,allfq=allfq))
                elif each.endswith('.fq'):
                        os.system('cat {each} >> {allfq}'.format(each=each,allfq=allfq))
                else:
                        print (each+" uncorrect format! ")
        assert os.path.isfile(allfq)
        fq2fa(allfq,allfa)
        assert os.path.isfile(allfa)
        os.system('formatdb -i {allfa} -p F'.format(allfa=allfa))

def fq2fa(fq,fa):
        out = open(fa,'w')
        i = 1
        j = 1
        with open(fq) as f:
                for line in f:
                        if i%4 == 1:
                                id = '>'+'0'*(10-len(str(j)))+str(j)+'\n'
                                out.write(id)
                                j += 1
                        elif i%4 == 2:
                                out.write(line)
                        i += 1
        out.close()

def trans(fasta):
	trans_out = ''
	trans_dic = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
	for i in range(len(fasta))[::-1]:
		trans_out += trans_dic[fasta[i]]
	return trans_out
	
def get_query():
	query_txt_all = ''
	f = open(cdsfasta,'r').read()
	tmp_list = f.strip().strip('>').split('>')
	for each in tmp_list:
		each_list = each.strip().split('\n')
		tmp_head = each_list[0]
		unigene_id = tmp_head.split(';')[0]
		tmp_fa = ''.join(each_list[1:])
		query_txt = '>'+unigene_id
		if terminal == '3':
			query_txt += ';3\n'
			query_txt += tmp_fa[-100:]+'\n'
		elif terminal == '5':
			query_txt += ';5\n'
			query_txt += trans(tmp_fa[:100])+'\n'
		query_txt_all += query_txt
		dic_query[unigene_id] = query_txt
		query_out = open(root+'/query.fasta','w')
		query_out.write(query_txt_all)
		query_out.close()
		f = ''
		
def prepare_dic():
	with open(allfa) as f:
                for line in f:
                        if line[0] == '>':
                                tmp_id = line[1:].strip()
                                dic_fa[tmp_id] = line
                        else:
                                dic_fa[tmp_id] += line

def blast_parse():
  	os.mkdir('blast_parse_out')
	os.system('/PUBLIC/software/public/Alignment/ncbi-blast-2.2.28+/bin/blastn -query query.fasta -out egw.out.xml -db DB.fa -outfmt 5 -evalue 1e-8 -num_threads 4')
	assert os.path.isfile(root+'/egw.out.xml')
	f = open(root+'/egw.out.xml','r').read()
	tmp_list = re.findall('<Iteration>[\s\S]+?</Iteration>',f)
	for each in tmp_list:
		m = re.search('<Iteration_query-def>(.+?)</Iteration_query-def>',each)
		query_id = m.group(1)
		query_id = query_id.split(';')[0]
		tmp_list = re.findall('<Hit_def>(\d+?)</Hit_def>',each)
		if len(tmp_list) > 0:
			bp_txt = dic_query[query_id]
			for each in tmp_list:
				bp_txt += dic_fa[each]
		open('blast_parse_out/'+query_id+'.blastout','w').write(bp_txt)

def mafft_align(blastout):
	mafftout = blastout.replace('blastout','mafftout')
	os.system('/BJPROJ/RNA/zhangfei/test/MAFFT/mafft-7.305-with-extensions/scripts/mafft --inputorder --adjustdirection --anysymbol --auto  {blastout} > {mafftout}'.format(blastout=blastout,mafftout=mafftout))

def prolong(mafftout):
#------- prepare
        pl_dic = {}
        f = open(mafftout,'r').readlines()
        i = 1
        for line in f:
                if line[0] == '>':
                        tmp_id = line[1:].strip()
                        pl_dic[tmp_id] = ''
                        if i == 1:
                                unigene_id = tmp_id
                        i += 1
                else:
                        pl_dic[tmp_id] += line.strip()
        uni_fa_long = pl_dic[unigene_id]
        uni_fa = uni_fa_long.replace('-','')
        del pl_dic[unigene_id]
#------- remove bad reads


#-------- get ** bp k-mer ( from tail **bp begin )
        count = 0
        for i in range(len(uni_fa_long))[::-1]:
                if uni_fa_long[i] != '-':
                        count += 1
                if count == philip:
                        sta = i
                        break

        dbj_list = []
        dbj_list.append(uni_fa_long[sta:].replace('-',''))

        for each_id in pl_dic.keys():
                pl_dic[each_id] = pl_dic[each_id][sta+1:].replace('-','')
                if len(pl_dic[each_id]) < philip:
                        del pl_dic[each_id]

        dic_kmer = {}
        while len(pl_dic.keys()) >= 3:
                for each_id in pl_dic.keys():
                        if pl_dic[each_id][:philip] in dic_kmer:
                                dic_kmer[pl_dic[each_id][:philip]] += 1
                        else:
                                dic_kmer[pl_dic[each_id][:philip]] = 1
                        pl_dic[each_id] = pl_dic[each_id][1:]
                        if len(pl_dic[each_id]) < philip:
                                del pl_dic[each_id]

        for each in dic_kmer.keys():
                if dic_kmer[each] < 3:
                        del dic_kmer[each]
#-------- prolong kmer (construct de Brujin graph, every prolong must >= 3 kmer )
	dic_kmer_use_count = {}
        flag = 'yes'
        while flag == 'yes':
                flag = 'no'
                for dbj in dbj_list:
                        new_dbj = dbj_list
                        for kmer in dic_kmer.keys():
                                if dbj_match(dbj,kmer) == 'TRUE':
                                        flag = 'yes'
                                        if dbj in new_dbj:
                                                new_dbj.remove(dbj)
                                        new_dbj.append(dbj+kmer[-1])
					if kmer not in dic_kmer_use_count:
						dic_kmer_use_count[kmer] = 0
					dic_kmer_use_count[kmer] += 1
					if dic_kmer_use_count[kmer] > 10:
						return unigene_id,uni_fa,uni_fa,''
                dbj_list = new_dbj

        max_len = 0
        max_dbj = ''
        for each in dbj_list:
                if len(each) > max_len:
                        max_len = len(each)
                        max_dbj = each
        return unigene_id,uni_fa,uni_fa[:-philip]+max_dbj,max_dbj[philip:]


def dbj_match(dbj,kmer):
        for i in range(philip-1):
                if kmer[i] != dbj[i-philip+1]:
                        return 'FALSE'
        return 'TRUE'


def update_cds(cdsfasta,egwout):
        upd_dic = {}
	longest = []
        f = open(egwout,'r').readlines()
        for line in f:
                if line[0] == '>':
                        tmp = line[1:].strip().split(';')
                        id = tmp[0]
                        term = tmp[1]
                        i = 1
                else:
                        i += 1
                        if i == 4:
                        	prolong = line.strip()
	                        if len(prolong) > 0:
					upd_dic[id] = [term,'']
	                                upd_dic[id][1] = prolong
				else:
					longest.append(id)

	longest_fa = ''
        cds_dic = {}
        f = open(cdsfasta,'r').readlines()
        for line in f:
                if line[0] == '>':
                        id = line[1:].strip().split(';')[0]
                        if id in upd_dic:
                                cds_dic[id] = [line,'']
			elif id in longest:
				longest_fa += '>'+id+'\n'
				id = 'longest'
			else:
				id = ''
                else:
			if id == '':
				pass
			elif id == 'longest':
				longest_fa += line.strip()
			else:
	                        cds_dic[id][1] += line.strip()
        update_out_part = ''
	update_out_comp = ''
	update_out = ''
        for id in upd_dic.keys():
                if upd_dic[id][0] == '3':
                        cds_dic[id][1] = cds_dic[id][1] + upd_dic[id][1]
                if upd_dic[id][0] == '5':
                        cds_dic[id][1] = trans(upd_dic[id][1]) + cds_dic[id][1]
		up_fa = cds_dic[id][1]
		up_fa_new = ''
                for i in range(len(up_fa)/60+1):
                        up_fa_new += up_fa[60*(i):60*(i+1)]+'\n'
		if orf_complete(up_fa,upd_dic[id][0]) == 'complete':
			update_out_comp += cds_dic[id][0]
			update_out_comp += up_fa_new
		else:
			update_out_part += cds_dic[id][0]
			update_out_part += up_fa_new	
		update_out += cds_dic[id][0]
		update_out += up_fa_new
        open('CDS_prolong.complete.fasta','w').write(update_out_comp)
	open('CDS_prolong.part.fasta','w').write(update_out_part)
	open('CDS_prolong.fasta','w').write(update_out)
	if longest_fa != '':
		longest_fa_new = ''
		for line in longest_fa.split('\n'):
			if line[0] == '>':
				longest_fa_new += line+'\n'
			else:
				for i in range(len(line)/60+1):
					longest_fa_new += line[i*60:(i+1)*60]+'\n'
		open('CDS_prolong_longest.fasta','a').write(longest_fa_new+'\n')

def orf_complete(fa,part):
        qs = ['ATG']
        zz = ['TAA','TAG','TGA']
        if part == '3':
                i = 0
                while i <= len(fa)-3:
                        if fa[i:i+3] in zz:
                                return 'complete'
                        i += 3
                return 'no'
        elif part == '5':
                i = len(fa)
                while i >= 3:
                        if fa[i-3:i] in qs:
                                return 'complete'
                        i -= 3
                return 'no'

########################################################################

prepare_database()
print 'ok-1'

dic_query = {}
get_query()
print 'ok-2'

allfq = root+'/DB.fq'
allfa = root+'/DB.fa'
dic_fa = {}
prepare_dic()

blast_parse()		
print 'ok-3'

all_blastout_list = glob.glob('blast_parse_out/*.blastout')
for blastout in all_blastout_list:
	mafft_align(blastout)
print 'ok-4'

out_txt = ''
philip = 20
all_mafftout_list = glob.glob('blast_parse_out/*.mafftout')	
for mafftout in all_mafftout_list:
	output_id,output_fa,output_fa_long,fa_prolong = prolong(mafftout)
	out_txt += '>'+output_id+'\n'+output_fa+'\n'+output_fa_long+'\n'+fa_prolong+'\n'
egwout = root+'/eGenomeWalking_output.xls'
open(egwout,'w').write(out_txt)

egwout = root+'/eGenomeWalking_output.xls'
update_cds(cdsfasta,egwout)

#----------- prolong to max
i = 1
while os.path.getsize('CDS_prolong.fasta') > 0:
	os.mkdir('temp_'+str(i))
	os.system('mv query.fasta egw.out.xml blast_parse_out eGenomeWalking_output.xls CDS_prolong.complete.fasta CDS_prolong.part.fasta temp_'+str(i))
	i += 1

	cdsfasta = 'CDS_prolong.fasta'
	get_query()

	blast_parse()

	all_blastout_list = glob.glob('blast_parse_out/*.blastout')
	for blastout in all_blastout_list:
		mafft_align(blastout)

	out_txt = ''
	philip = 20
	all_mafftout_list = glob.glob('blast_parse_out/*.mafftout')

	for mafftout in all_mafftout_list:
		output_id,output_fa,output_fa_long,fa_prolong = prolong(mafftout)
		out_txt += '>'+output_id+'\n'+output_fa+'\n'+output_fa_long+'\n'+fa_prolong+'\n'
	egwout = root+'/eGenomeWalking_output.xls'
	open(egwout,'w').write(out_txt)
	egwout = root+'/eGenomeWalking_output.xls'
	update_cds(cdsfasta,egwout)

