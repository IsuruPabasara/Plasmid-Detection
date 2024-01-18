clc;
clear all;

//find best positions
function [positions,all_scores]=find_positions(filename,seg_of_c)
    n=10;
    score_mat=read(filename,-1,seg_of_c);
    selected_scores=gsort(unique(score_mat))(1:n);
    all_scores=[]
    positions=[]
    for q=1:n
        positions=[positions,find(score_mat==selected_scores(q))]
        for k=1:length(find(score_mat==selected_scores(q)))
            all_scores=[all_scores,selected_scores(q)]
        end
    end
endfunction

//find distance matrix
function [distance_mat]=calc_distance_matrix(start_positions,end_positions,seg_of_c)
    no_of_start_pos=length(start_positions);
    no_of_end_pos=length(end_positions);
    distance_mat=zeros(no_of_end_pos,no_of_start_pos)
    for q=1:no_of_end_pos
        for p=1:no_of_start_pos
            if start_positions(p)<end_positions(q)
                distance_mat(q,p)=end_positions(q)-start_positions(p)
            else
                distance_mat(q,p)=seg_of_c-(start_positions(p)-end_positions(q))
            end        
        end
    end
endfunction

//Remove EOL, overhead and segment the FASTA considering circular nature
function [gen_seg_M,gen_array]=get_segments(file_name,m)
      fd = mopen(file_name,'rb');
      file_details = fileinfo(file_name);
      file_len = file_details(1);
      mseek(0,fd);
      header = mgetl(fd,1);
      g_start = length(header)+1;      
      mseek(g_start,fd);      
      raw_code = mget(file_len-g_start,'c',fd);      
      mclose(fd);
      code_i=raw_code(~members(raw_code,[10]));
      disp(length(code_i))
      gen_array=code_i;      
      n=length(code_i);
      col=floor(n/m);
      gen_seg_M=matrix(code_i(1:m*col),m,col)';
endfunction

function comp_mat = build_comp_mat(x, y, score_m, score_mm)
   comp_mat = (~(repmat(x', 1, length(y)) - repmat(y, length(x), 1)))*(score_m - score_mm) + score_mm;
   comp_mat = flipdim(comp_mat,1);
endfunction

//Calculate an alignment score based on local alignment
function matrix_result = scorematrix_local(seq_x,seq_y,score_m,score_mm,gap_penalty)
    len_x = length(seq_x);
    len_y = length(seq_y);
    basic_mat = zeros(len_x+1,len_y+1);
    comp_mat = build_comp_mat(seq_x, seq_y, score_m, score_mm);
    for c = -len_x+1:len_y-1
        gap_p = (c < 1)*1;
        gap_q = (c > len_y - len_x)*1;
        
        score_match = diag(basic_mat,c-1)(2 - gap_p:$-gap_q)+diag(comp_mat,c);
        score_gap_x = diag(basic_mat,c)(1:$-1)+gap_penalty;
        score_gap_y = diag(basic_mat,c)(2:$)+gap_penalty;
        max_match = max(score_match, score_gap_x, score_gap_y, 0);
        
        l = length(max_match);
        index = 1+max(0,-c) + (2+max(0,c) - 1)*(len_x+1):len_x+2:max(0,-c)+l + (1+max(0,c)+l - 1)*(len_x+1);
        basic_mat(index) = max_match;
    end
    matrix_result = max(basic_mat);
endfunction

chrom_file_name='D:\Work\Aca Sem 7\Genomic Signal Processing\Assignmnet 2\E_C_Genome.fasta';
plas_file_name='D:\Work\Aca Sem 7\Genomic Signal Processing\Assignmnet 2\E_C_Plasmid_Short.fasta';
     
n_segc=1000; //elements in a segment
n_segp=200;

[chrom_segments,chrom_array]=get_segments(chrom_file_name,n_segc);    //segmenting the genome
[plas_segments,plas_array]=get_segments(plas_file_name,n_segp);    //segmenting the plasmid

seg_of_c=floor(length(chrom_segments)/n_segc); //no of segments of chromosome
seg_of_p=floor(length(plas_segments)/n_segp);  //no of segments of plasmid

//comparing first and last segments of plasmid with the chromosome segments
score_mat_start=[]  
score_mat_end=[]
for k=1:seg_of_c
    if modulo(k,10)==0
        disp(k)
    end
    score = scorematrix_local(chrom_segments(k,:),plas_segments(1,:),2,-2,-3);
    score_mat_start=[score_mat_start,score];
    score = scorematrix_local(chrom_segments(k,:),plas_segments(seg_of_p,:),2,-2,-3);
    score_mat_end=[score_mat_end,score];
end
filename='D:\Work\Aca Sem 7\Genomic Signal Processing\Assignmnet 2\Score_Mat_Start.csv';
deletefile(filename);
write(filename,score_mat_start);
filename='D:\Work\Aca Sem 7\Genomic Signal Processing\Assignmnet 2\Score_Mat_End.csv';
deletefile(filename);
write(filename,score_mat_end);

//Finding the alignment positions of best end alignments
filename='D:\Work\Aca Sem 7\Genomic Signal Processing\Assignmnet 2\Score_Mat_Start.csv';
[start_positions_unsorted,all_start_scores]=find_positions(filename,seg_of_c);
start_positions=gsort(start_positions_unsorted,'g','i');
scatter(start_positions_unsorted,all_start_scores);
filename='D:\Work\Aca Sem 7\Genomic Signal Processing\Assignmnet 2\Score_Mat_End.csv';
[end_positions_unsorted,all_end_scores]=find_positions(filename,seg_of_c);
end_positions=gsort(end_positions_unsorted,'g','i');
scatter(end_positions_unsorted,all_end_scores);

//Finding reigions of interest
filename='D:\Work\Aca Sem 7\Genomic Signal Processing\Assignmnet 2\Distance_Mat.csv';
deletefile(filename);
distance_mat=calc_distance_matrix(start_positions,end_positions,seg_of_c)
write(filename,distance_mat);
Matplot(flipdim(distance_mat,1));
gcf().color_map=jetcolormap(200);

//Creating FASTA files of reigion of interest
start_position_index=[9,13,20,31,36]   //manually chosen values
end_position_index=[23,25,31,40,47]
no_of_reg=length(start_position_index)
for i=1:no_of_reg
    filename='D:\Work\Aca Sem 7\Genomic Signal Processing\Assignmnet 2\Reigion_Of_Interest_'+string(i)+'.FASTA'
    fd = mopen(filename,'wb');    
    reg_start_index=((start_positions(start_position_index(i))-1)*n_segc)+1;
    reg_end_index = end_positions(end_position_index(i))*n_segc;
    disp("csada")
    disp(reg_start_index,reg_end_index)    
    header=">Reigion"+string(i)
    mputl(header,fd)
    save_reg=chrom_array(reg_start_index:reg_end_index)
    k=1
    for j=reg_start_index:reg_end_index
        mput(chrom_array(j),'c')
        if modulo(k,70)==0
            mputl('',fd)
        end
        k=k+1
    end
    mputl('',fd)
    mclose(fd);
end
