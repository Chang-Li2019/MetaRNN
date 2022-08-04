import os, re, h5py, sys, allel
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
import numpy as np
import tensorflow as tf

import pandas as pd
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

# input file preparation for annovar annotation
in_vcf = allel.read_vcf(str(sys.argv[2]))
try:
    genome_build_version = str(sys.argv[1])
    in_vcf = allel.read_vcf(str(sys.argv[2]))
except:
    print('Use format: ./MetaRNN genome_build(hg19 or hg38) input_vcf\nExample:\n./MetaRNN hg38 test.vcf')

else:

    annovar_in = open(str(sys.argv[2])+'.temp','w')
    nssnv_dict = {}
    CHROMS = {}
    old_cor_dict = {} 
    # Convert vcf coordinates to ANNOVAR coordinates for indels;
    # Save SNVs to a dictionary for annotation
    for i in range(len(in_vcf['variants/CHROM'])):
        pos = in_vcf['variants/POS'][i]
        alts = [a for a in in_vcf['variants/ALT'][i] if a !='']
        chrom = str(in_vcf['variants/CHROM'][i]).upper().replace('CHR','')
        for alt in alts:
            if len(alt)==1 and len(in_vcf['variants/REF'][i]) <50 and len(in_vcf['variants/REF'][i])>1:
                l_pos = str(pos +1)
                r_pos = str(int(l_pos) +len(in_vcf['variants/REF'][i])- 2)
                old_cor = '_'.join([chrom, str(pos), in_vcf['variants/REF'][i], alt])
                alt = '-'
                ref = in_vcf['variants/REF'][i][1:]
                old_cor_dict['_'.join([chrom,l_pos, r_pos, ref, alt])] = old_cor
                outline = '\t'.join([chrom,l_pos, r_pos, ref, alt])+'\n'
                annovar_in.write(outline)
            elif len(in_vcf['variants/REF'][i]) == 1 and len(alt) <50 and len(alt)>1:
                l_pos = str(pos +1)
                r_pos = str(pos +1)
                old_cor = '_'.join([chrom, str(pos), in_vcf['variants/REF'][i], alt])
                ref = '-'
                alt = alt[1:]
                old_cor_dict['_'.join([chrom,l_pos, r_pos, ref, alt])] = old_cor
                outline = '\t'.join([chrom,
                    l_pos, r_pos, ref, alt])+'\n'
                annovar_in.write(outline)
            elif len(in_vcf['variants/REF'][i])==1 and len(alt)==1:
                key = '_'.join([chrom,str(pos),str(in_vcf['variants/REF'][i]),str(alt)])
                CHROMS[in_vcf['variants/CHROM'][i]] = 0
                nssnv_dict[key] = 0
            elif len(in_vcf['variants/REF'][i]) < 50 and len(alt) <50:
                l_pos = str(pos) 
                r_pos = str(int(l_pos) + len(in_vcf['variants/REF'][i]) -1)
                old_cor = '_'.join([chrom, str(pos), in_vcf['variants/REF'][i], alt])
                ref = in_vcf['variants/REF'][i]
                old_cor_dict['_'.join([chrom,l_pos, r_pos, ref, alt])] = old_cor
                outline = '\t'.join([chrom,
                    l_pos, r_pos, ref, alt])+'\n'
                annovar_in.write(outline)
    annovar_in.close()

    # Processing INDELs
    if "l_pos" in locals():
        # get annovar annotation and load MetaRNN-indel model
        input_file = str(sys.argv[2])+'.temp'
        anno_exe = "./annotate_variation.pl -geneanno -dbtype ensGene --transcript_function  -buildver {} {} humandb/ --separate".format(genome_build_version,input_file)
        os.system(anno_exe)
        Annotation = open('{}.exonic_variant_function'.format(input_file),'r')
        model = tf.keras.models.load_model('RNN_DNN_indel_adam.h5', compile=False)
        model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy','AUC'])
        print('#######################################')
        print('Annotating INDELs ...')
        Annotated = open('{}.indel.temp'.format(str(sys.argv[2])),'w')
        Annotated.write('\t'.join(['#CHROM','POS','REF','ALT','TRANSCRIPT','MetaRNN-indel_score'])+'\n')
        # Save INDEL variants' IDs
        indel_var = open(input_file,'r')
        indel_var_dict = {}
        for idx,line in enumerate(indel_var):
            idx += 1
            line = line.strip('\n').split('\t')
            line[0] = line[0].upper().replace('CHR','')
            value = '_'.join(line)
            indel_var_dict[idx] = value
        indel_var.close()
        # Annotate INDELS line by line 
        for line in Annotation:
            line = line.strip('\n').split('\t')
            if 'nonframeshift' in line[1]:
                idx = line[0].replace('line','')
                info = [x for x in line[2].split(',') if x != '']
                for enst in info:
                    enst = enst.split(':')
                    if len(enst)>4:
                        trans_name = re.sub('\..+','',enst[1])
                        pos = re.sub('[A-Za-z]+','',enst[4])
                        pos = re.split('\D',enst[4])
                        pos = [int(x) for x in pos if x!='']
                        l_pos = max(0,min(pos))
                        r_pos = max(pos)
                        key = indel_var_dict[int(idx)].split('_')
                        anno_file = h5py.File('./Transcript_Anno/Chr{}.hdf5'.format(key[0]),'r')
                        if trans_name in anno_file.keys():
                            anno = anno_file[trans_name][()]
                            anno = anno[(l_pos<=anno[:,-1]) & (anno[:,-1]<=r_pos),:]
                            input_instance = np.pad(anno,((0,58-anno.shape[0]),(0,0)))[:,:-1]
                            score = model.predict(np.expand_dims(input_instance,axis=0)).flatten()[0]
                            new_cor = '_'.join(key)
                            key = old_cor_dict[new_cor].split('_')
                            Annotated.write('\t'.join(key)+'\t'+trans_name+'\t'+str(score)+'\n')
        Annotated.close()
        enst_dict = {}
        score_dict = {}
        
        # reformat annotated INDELs into each variant per row.
        with open('{}.indel.temp'.format(str(sys.argv[2])),'r') as f:
            for line in f:
                line = line.strip('\n').split('\t')
                key = '_'.join(line[:4])
                if key not in enst_dict:
                    enst_dict[key] = [line[4]]
                    score_dict[key] = [line[5]]
                else:
                    enst_dict[key].append(line[4])
                    score_dict[key].append(line[5])
        with open('{}.indel.annotated'.format(str(sys.argv[2])),'w') as f:
            for key,val in enst_dict.items():
                coords = key.split('_')
                scores = ';'.join(score_dict[key])
                ensts = ';'.join(val)
                outstring = '\t'.join(coords)+'\t'+ensts+'\t'+scores+'\n'
                f.write(outstring)
        print('Finished INDEL annotation! Annotated file written to {}.indel.annotated'.format(str(sys.argv[2])))
        os.system(f'rm {str(sys.argv[2])}.indel.temp')

    else:
        print('#######################################')
        print('No indels found!')
    
    
    
    if nssnv_dict=={}:
        print('No SNVs found! Analysis finished!')
    else:
        print('Annotating nsSNVs ...')
        for loop_num,i in enumerate(CHROMS.keys()):
            i = i.upper().replace('CHR','')
            db = pd.read_pickle(f'./SNV_Anno/Chr{i}.pkl','gzip')
            if loop_num ==0:
                if genome_build_version=='hg38':
                    key = db.iloc[:,0].astype(str)+'_'+db.iloc[:,1].astype(str)+'_'+db.iloc[:,2].astype(str)+'_'+db.iloc[:,3]
                    idx1 = key.isin(nssnv_dict)
                    # remove SNVs that are not nonsynonymous
                    idx2 = -(db.iloc[:,4]=='X')   
                    idx3 = -(db.iloc[:,5]=='X')
                    idx4 = ((db.iloc[:,4]+'_'+db.iloc[:,5])!='._.')
                    idx = idx1 & idx2 & idx3 & idx4
                    db.loc[idx,['#chr', 'pos(1-based)', 'ref', 'alt','Ensembl_transcriptid', 'MetaRNN_score']].to_csv(
                        f'{sys.argv[2]}.nsSNV.annotated',sep='\t',header=['#CHROM','POS','REF','ALT','TRANSCRIPT','MetaRNN_score'],index=False)

                elif genome_build_version=='hg19':
                    key = db.iloc[:,6].astype(str)+'_'+db.iloc[:,7].astype(str)+'_'+db.iloc[:,2].astype(str)+'_'+db.iloc[:,3]
                    idx1 = key.isin(nssnv_dict)
                    idx2 = -(db.iloc[:,4]=='X')
                    idx3 = -(db.iloc[:,5]=='X')
                    idx4 = ((db.iloc[:,4]+'_'+db.iloc[:,5])!='._.')
                    idx = idx1 & idx2 & idx3 & idx4
                    db.loc[idx,['hg19_chr', 'hg19_pos(1-based)', 'ref', 'alt','Ensembl_transcriptid', 'MetaRNN_score']].to_csv(
                        f'{sys.argv[2]}.nsSNV.annotated',sep='\t',header=['#CHROM','POS','REF','ALT','TRANSCRIPT','MetaRNN_score'],index=False)
            else:
                if genome_build_version=='hg38':
                    key = db.iloc[:,0].astype(str)+'_'+db.iloc[:,1].astype(str)+'_'+db.iloc[:,2].astype(str)+'_'+db.iloc[:,3]
                    idx1 = key.isin(nssnv_dict)
                    idx2 = -(db.iloc[:,4]=='X')
                    idx3 = -(db.iloc[:,5]=='X')
                    idx4 = ((db.iloc[:,4]+'_'+db.iloc[:,5])!='._.')
                    idx = idx1 & idx2 & idx3 & idx4
                    db.loc[idx,['#chr', 'pos(1-based)', 'ref', 'alt','Ensembl_transcriptid', 'MetaRNN_score']].to_csv(
                        f'{sys.argv[2]}.nsSNV.annotated',sep='\t',mode = 'a',header=None,index=False)

                elif genome_build_version=='hg19':
                    key = db.iloc[:,6].astype(str)+'_'+db.iloc[:,7].astype(str)+'_'+db.iloc[:,2].astype(str)+'_'+db.iloc[:,3]
                    idx1 = key.isin(nssnv_dict)
                    idx2 = -(db.iloc[:,4]=='X')
                    idx3 = -(db.iloc[:,5]=='X')
                    idx4 = ((db.iloc[:,4]+'_'+db.iloc[:,5])!='._.')
                    idx = idx1 & idx2 & idx3 & idx4
                    db.loc[idx,['hg19_chr', 'hg19_pos(1-based)', 'ref', 'alt','Ensembl_transcriptid', 'MetaRNN_score']].to_csv(
                        f'{sys.argv[2]}.nsSNV.annotated',sep='\t',mode = 'a',header=None,index=False)
        print('Finished nsSNV annotation! Annotated file written to {}.nsSNV.annotated'.format(str(sys.argv[2])))
    
    print('#######################################')