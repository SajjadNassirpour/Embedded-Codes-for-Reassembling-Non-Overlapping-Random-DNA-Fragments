function C=VT_encoder(data_bits,r_desired,k)
    %===== This function performs VT encoder =====%
    
    n=2^k;  % length of each VT codeword section
    j=0;
    jj=1;
    C=zeros(1,n);
    for i=1:n
        if i==2^j
           C(1,i)=0;
           j=j+1;
        else
           C(1,i)=data_bits(1,jj); %Initial version of VT codeword section
           jj=jj+1;
        end
    end      
    
    %=== Calculate Varshamov-Tenengolts summation ===%
    VT_sum=0;
    for i=1:n
        VT_sum=VT_sum+(i*C(1,i));
    end
    r=mod(VT_sum,(n+1)); % Residue of intial VT codeword
    
    %=== Modify remainder value r to achieve desired residue ===%
    C_modified=C;
    if r~=r_desired
        if (n+1-r+r_desired)<2^(k+1)
            modify_remainder=de2bi((n+1-r+r_desired),(k+1));
        else
            modify_remainder=de2bi((n-1),(k+1));
        end
        j=0;
        for i=1:n
            if i==2^j
               C_modified(1,i)=modify_remainder(1,(j+1));
               j=j+1;      
            end
        end  
    end
    C=C_modified; % VT codeword section
end