function C_seperated=DNA_frag(C,ran_vec,cut_rule)
    %===== This function generates set of DNA fragments =====%
    
    n=length(C); % Length of VT codeword
    section_number=length(ran_vec)+1;  % Number of DNA fragments
    C_seperated=zeros(section_number,floor(cut_rule)); 
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
end