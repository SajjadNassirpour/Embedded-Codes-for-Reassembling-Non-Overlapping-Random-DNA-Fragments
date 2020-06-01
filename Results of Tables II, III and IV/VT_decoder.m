function [A,B,diff]=VT_decoder(C_DNA,C_sep_1,C_sep_2,k1,k2,m1,section_number)
    %===== This function performs VT-decoder =====%
    
    A=0;  %1st Reconstructed VT codeword substring as x^(1)
    B=0;  %2nd Reconstructed VT codeword substring as x^(2)
    n=size(C_DNA,2);
    n_sec=2^k1;  % Length of VT codeword section
    diff=zeros(2,1);
    q=1;
    %=== VT-decoder procedure ===%
    while q<m1+1 & sum(diff)<2
        diff=zeros(2,1);
        comb=randperm(section_number);
        n1=section_number;
        for DN=1:2
            if DN==1
                C_seperated=C_sep_1;
            else
                C_seperated=C_sep_2;
            end
            C=C_DNA(DN,:);
            po_0=0;
            
            %=== Reconstruct VT codeword from DNA fragments ===%
            for qq=1:n1
                section_indicator=comb(qq);
                section_data_size=C_seperated(section_indicator,1);
                po_1=po_0+section_data_size;
                C_combined((po_0+1):po_1)=C_seperated(section_indicator,2:(section_data_size+1));
                po_0=po_1;
            end
            sec_ind=1;
            r_break=0;
            
            %=== Check VT conditions ===%
            while sec_ind<(2^k2)+1 && r_break==0
                r_desired=sec_ind-1;
                C_combined_sec=C_combined(((sec_ind-1)*n_sec+1):sec_ind*n_sec);
                VT_sum=0;
                for i=1:n_sec
                    VT_sum=VT_sum+(i*C_combined_sec(1,i));
                end
                r_combined=mod(VT_sum,(n_sec+1));
                if r_combined~=r_desired
                    r_break=1;
                end
                sec_ind=sec_ind+1;
            end
            
            %= If all VT conditions are satisfied, then we decode data =%
            if r_break==0
                ii=1;
                jj=0;
                while ii<n+1 && diff(DN)==0
                    if ii==2^jj
                       ii=ii+1;
                       jj=jj+1;
                    else
                        if C_combined(1,ii)~=C(1,ii) 

                            if DN==1
                               A =comb;
                            else
                                B=comb;
                            end
                            diff(DN)=1;
                        else
                            ii=ii+1;
                        end
                    end
                end
            end
        end
        q=q+1;   
    end
end