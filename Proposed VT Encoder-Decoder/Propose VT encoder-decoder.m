%=== This Simulation performs our proposed VT encoder-decoder scheme ===%
clc
clear all
close all
%% Parameters
data_num=2*10000000; % Number of different data-bit strings 
k1=4; % 2^k1: Length of each VT codeword section
disp('k1:');
disp(k1);
k2=4;  % 2^k2: Number of VT codeword sections
disp('k2:');
disp(k2);
p=0.5; % Bernoulli probability of each element in data-bit string d
data_len=(2^k2)*(2^k1-k1-1); % Length of each data bit sub-string
n=(2^k2)*(2^k1);  %Length of DNA string
disp('DNA Length:');
disp(n);
cut_rule=log2(n); % Breaking pattern
confined_iters=1000000000;  % Limit for decoding iterations
error_new=0; % Number of cases that error happens
algo_num_1=0; % decoding iterations

%% Main-procedure
for loop=1:data_num
    error_happen=0;
    algo_num=0;
    for DN=1:2
        %=== Generate data-bit string ===%
        C=[];
        data_bits_total=[];
        for i=1:2^k2
            ni=2^k1;
            data_bits_length=2^k1-k1-1;
            data_bits=(rand(1,data_bits_length)<p);
            data_bits_total=[data_bits_total data_bits];
            
            %===== Decoding Procedure =====%
            
            %= Initial VT codeword =%
            ci=zeros(1,ni);
            j=0;
            jj=1;
            for pp=1:ni
                if pp==2^j
                   ci(1,pp)=0;
                   j=j+1;
                else
                   ci(1,pp)=data_bits(1,jj);
                   jj=jj+1;
                end
            end      
            
            VT_sum=0;
            for pp=1:ni
                VT_sum=VT_sum+(pp*ci(1,pp));
            end
            r=mod(VT_sum,(ni+1));  % Residue of initial VT code
            
            %=== Modify remainder value r to achieve desired residue ===%
            if r~=(i-1)
                if (ni+1-r+(i-1))<2^(k1+1)
                modify_remainder=de2bi((ni+1-r+(i-1)),(k1+1));
                else
                  modify_remainder=de2bi((ni-1),(k1+1));
                end
                ci_modified=ci;
                j=0;
                for pp=1:ni
                    if pp==2^j
                       ci_modified(1,pp)=modify_remainder(1,(j+1));
                       j=j+1;      
                    end
                end 
                ci=ci_modified;
            end    
            C=[C ci]; % VT codeword sub-string
        end
        
        %=== Generate breaking points of DNA string ===%
        N=floor(cut_rule)-1;
        if DN==1
            ran_vec=geornd(cut_rule/n,1,2*(N+1));
            if size(find(ran_vec==0),2)~=0
                ran_vec=ran_vec+ones(1,2*(N+1));
            end
            ran_vec; 

            pre_rran=0;
            rran=1;
            stop_ran=0;
            rran_vec=0;
            while stop_ran==0 & rran<2*(N+1)+1
                if ran_vec(1,rran)+pre_rran<n
                   rran_vec(1,rran)=ran_vec(1,rran)+pre_rran;
                   pre_rran=rran_vec(1,rran);
                   rran=rran+1;
                else
                   stop_ran=1;
                end
            end
            ran_vec;
            ran_vec=rran_vec;
        end
         
            if ran_vec==0
            else
                %=== Generate set of DNA fragments ===%
               section_number=length(ran_vec)+1;
               half_algo=algo_num;
               C_seperated=zeros(section_number,n);
               po_0=0;
                for i=1:(section_number-1)
                    po_1=ran_vec(1,i);
                    C_seperated(i,1)=po_1-po_0;
                    C_seperated(i,2:(po_1-po_0+1))=C(1,(po_0+1):po_1);
                    po_0=po_1;
                end
                po_1=n;
                C_seperated(section_number,1)=po_1-po_0;
                C_seperated(section_number,2:(po_1-po_0+1))=C(1,(po_0+1):po_1);
                
                
                %===== Decoding Procedure =====%
                
                %=== Define controlling parameters ===%
                % section_number+1= Variable V5
                % section_number+2= Variale V3
                % section_number+3= Flag F2
                % section_number+4= Flag F1
                % num_mu_case= Variable V4
                
                loc=zeros(1,section_number+2);                  
                loc(:,section_number+3)=-1;  
                loc(:,section_number+4)=1;   
                
                %=== Algorithms 2 and 3 ===%
                while sum(loc(:,section_number+4))~=0 & algo_num<confined_iters+1

                      u=size(loc,1);
                      row=u;
                      for i=1:u
                          if loc(i,section_number+3)==-1 && loc(i,section_number+4)==1
                              pos=loc(i,section_number+1);
                              ni_current=loc(i,section_number+2);

                              temp=1:section_number;
                              if pos~=0
                                  temp2=loc(i,1:pos);
                                  for j=1:section_number
                                      for jj=1:pos
                                          if temp(1,j)==temp2(1,jj)
                                              temp(j)=0;
                                          end
                                      end
                                  end
                              end
                              num_mu_case=1;
                              for j=1:section_number
                                  if temp(1,j)~=0
                                      algo_num=algo_num+1;
                                      
                                      first_sec=[];
                                      if pos~=0
                                          for jj=1:pos
                                              i_num=loc(i,jj);
                                              num=C_seperated(i_num,1);
                                              first_sec=[first_sec C_seperated(i_num,2:(num+1))];
                                          end
                                      end
                                      
                                      num2=C_seperated(j,1);
                                      second_sec=C_seperated(j,2:(num2+1));
                                      test_sec=[first_sec second_sec];

                                      ni_number=floor(length(test_sec)/ni)-ni_current;
                                      
                                      %=== Algorithm 2 ===%
                                      if ni_number~=0
                                          test_sec=test_sec(1,ni_current*ni+1:(ni_current+ni_number)*ni);
                                          case_correct=1;
                                          for ni_c=1:ni_number
                                              test_sec1=test_sec((ni_c-1)*ni+1:ni_c*ni);
                                              VT_sum=0;
                                              for pp=1:ni
                                                  VT_sum=VT_sum+(pp*test_sec1(1,pp));
                                              end
                                              r=mod(VT_sum,(ni+1));

                                              if r==(ni_current+ni_c-1) && case_correct==1
                                                  case_correct=1;
                                              else
                                                  case_correct=0;
                                              end
                                          end
                                          
                                          if case_correct==1
                                              %== If new row has common
                                             %component with other rows ==%
                                              if num_mu_case==1
                                                  loc(i,pos+1)=j;
                                                  loc(i,section_number+1)=loc(i,section_number+1)+1;
                                                  loc(i,section_number+2)=loc(i,section_number+2)+ni_number;
                                                  ni_number_temp=ni_number;
                                                  num_mu_case=num_mu_case+1;
                                              else
                                                  row=row+1;
                                                  loc(row,1:section_number+4)=0;
                                                  if pos~=0
                                                  loc(row,1:pos)=loc(i,1:pos);
                                                  end
                                                  loc(row,pos+1)=j;
                                                  loc(row,section_number+1)=loc(i,section_number+1);
                                                  loc(row,section_number+2)=loc(i,section_number+2)+ni_number-ni_number_temp;
                                                  loc(row,section_number+3)=-1;
                                                  loc(row,section_number+4)=1;
                                              end
                                          else
                                              if num_mu_case==1
                                                  loc(i,section_number+3)=-2;
                                                  loc(i,section_number+4)=0;
                                              end
                                          end
                                      end
                                      if ni_number==0
                                          %== If new row has common
                                          %component with other rows ==%
                                          if num_mu_case==1
                                             loc(i,pos+1)=j;
                                             loc(i,section_number+1)=loc(i,section_number+1)+1;
                                             ni_number_temp=ni_number;
                                             num_mu_case=num_mu_case+1;
                                          else                                            
                                             row=row+1;
                                             loc(row,1:section_number+4)=0;
                                             if pos~=0
                                             loc(row,1:pos)=loc(i,1:pos);
                                             end
                                             loc(row,pos+1)=j;
                                             loc(row,section_number+1)=pos+1;
                                             loc(row,section_number+2)=loc(i,section_number+2)+ni_number-ni_number_temp;
                                             loc(row,section_number+3)=-1;
                                             loc(row,section_number+4)=1;                                             
                                          end
                                      end
                                  end
                              end
                          end
                          if loc(i,section_number)~=0
                              loc(i,section_number+4)=0;
                          end
                      end
                end
                
              %== Find all permutations which satisfy all VT conditions ==%
                if algo_num>confined_iters
                else
                    loc1=loc(:,1:section_number);
                    loc=loc1;
                   
                    satisfied_cases=[];
                    ca=0;
                    for i=1:size(loc,1)
                        if loc(i,section_number)~=0
                            ca=ca+1;
                           satisfied_cases(ca,:)=loc(i,:);
                        end
                    end
                    algo_num_1=algo_num_1+algo_num;
                    Data_line(DN,:)=data_bits_total;
                    VT_code(DN,:)=C;
                    if DN==1
                        satisfied_seq_1=satisfied_cases;
                        C_sep_1=C_seperated;
                        size_satisfied_seq_1=[];
                        for s1_q1=1:size(satisfied_seq_1,1)
                            for s1_q2=1:size(satisfied_seq_1,2)
                                size_satisfied_seq_1(s1_q1,s1_q2)=C_seperated(satisfied_seq_1(s1_q1,s1_q2),1);
                            end
                        end
                    else
                        satisfied_seq_2=satisfied_cases;
                        C_sep_2=C_seperated;
                        size_satisfied_seq_2=[];
                        for s2_q1=1:size(satisfied_seq_2,1)
                            for s2_q2=1:size(satisfied_seq_2,2)
                                size_satisfied_seq_2(s2_q1,s2_q2)=C_seperated(satisfied_seq_2(s2_q1,s2_q2),1);
                            end
                        end
                    end
                end
            end
    end
    
    %===== Error Calculation =====%
    if algo_num>confined_iters
    else
        if ran_vec==0
        else
            for s_q1=1:size(size_satisfied_seq_1,1)
                for s_q2=1:size(size_satisfied_seq_2,1)
                    if norm(satisfied_seq_1(s_q1,:)-satisfied_seq_2(s_q2,:),2)==0
                        % 1st reconstruced VT codeword sub-string x^(1) 
                        first_sec1=[];
                        for j_c=1:section_number
                            i_num=satisfied_seq_1(s_q1,j_c);
                            num=C_sep_1(i_num,1);
                            first_sec1=[first_sec1 C_sep_1(i_num,2:(num+1))];
                        end 
                        % 2nd reconstruced VT codeword sub-string x^(2)
                        first_sec2=[];
                        for j_c=1:section_number
                            i_num=satisfied_seq_2(s_q2,j_c);
                            num=C_sep_2(i_num,1);
                            first_sec2=[first_sec2 C_sep_2(i_num,2:(num+1))];
                        end
                        
                        %== Error Checking ==%
                        if error_happen==0
                            if norm(first_sec1-VT_code(1,:),2)~=0 | norm(first_sec2-VT_code(2,:),2)~=0
                                error_new=error_new+1;
                                error_happen=1;
                            error_seq_1=satisfied_seq_1(s_q1,:);
                            error_seq_2=satisfied_seq_2(s_q2,:);
                            end
                        end
                    end
                end
            end

        end
    end
end

%% Average Decoding Error and Average iterations/factorial(DNA fragments)
final_error_new=error_new/data_num;
disp('Average Error Rate:');
disp(final_error_new);
algo_average=algo_num_1/(factorial(log2(n))*(data_num));
disp('Average Decoding Iterations:');
disp(algo_average);

