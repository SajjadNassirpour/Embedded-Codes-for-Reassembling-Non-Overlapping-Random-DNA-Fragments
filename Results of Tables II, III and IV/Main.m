clc
clear all
close all

%% Parameters
k1=8;  % 2^k1: Length of each VT codeword section
disp('k1:');
disp(k1);
k2=0;  % 2^k2: Number of VT codeword sections
disp('k2:');
disp(k2);
p=0.5; % Bernoulli probability of each element in data-bit string d
data_num=10000;  % Number of different data-bit strings 
chs_comb=100000000;  % Number of permutations (less than factorial(DNA fragments))
n=(2^k2)*(2^k1); % DNA Length
disp('DNA Length:');
disp(n);
cut_rule=log2(n); % Breaking Pattern
error=0; % Number of cases that error happens

%% Main-procedure
for loop1=1:data_num
    
    data_DNA=[]; %data-bit strings
    C_DNA=[]; % VT codeword 
    
    for DN=1:2
        
        %=== Generate data bits and VT encoder ===%
        data_bits_total=[];
        C=[];
        for kk=1:(2^k2)
            data_bits=rand_gen(k1,p);
            data_bits_total=[data_bits_total data_bits];
            r_desired=kk-1;  %r_i=i-1
            C_section=VT_encoder(data_bits,r_desired,k1);  % Generate VT codeword section
            C=[C C_section];  % Merge all VT codeword sections
        end
        C_DNA(DN,:)=C;
        data_DNA(DN,:)=data_bits_total;
        
        %=== Break DNA string into some DNA fragments ===%
        if DN==1
            ran_vec=break_points(cut_rule,n); % Generate breaking points
        end

        if ran_vec==0
        else
            C_seperated=DNA_frag(C,ran_vec,cut_rule); % Set of DNA fragments
            if DN==1
                C_sep_1=C_seperated;
            else
                C_sep_2=C_seperated;
            end
        end
    end
    diff=zeros(2,1);
    section_number=length(ran_vec)+1;  % Number of DNA fragments
    
    %=== VT decoder ===%
    if ran_vec==0
    else
        m1=min(factorial(section_number),chs_comb);  % Check min(factorial(fragments),chs_comb)
        [A,B,diff]=VT_decoder(C_DNA,C_sep_1,C_sep_2,k1,k2,m1,section_number); % VT decoder
    end  
    
    %=== Error Calculation ===%
    if ran_vec==0 | sum(diff)<2
    else
        no_error=section_number:-1:1;
        error_finder=norm(no_error-A);

        if error_finder~=0 
            error=error+1;
        end
    end        
end
err=error/data_num;
disp('Average Error Rate:');
disp(err);


