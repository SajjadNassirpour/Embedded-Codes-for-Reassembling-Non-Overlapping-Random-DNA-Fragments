function ran_vec=break_points(cut_rule,n)
    %===== This function generates Breaking points =====%
    %===== of DNA string based on breaking pattern =====%

    N=floor(cut_rule)-1;  %Number of points
    ran_vec=geornd(cut_rule/n,1,2*(N+1)); % Generate  Geo distribution
    if size(find(ran_vec==0),2)~=0
        ran_vec=ran_vec+ones(1,2*(N+1));
    end
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
    ran_vec=rran_vec; % Breaking points
end