function rand_seq = rand_gen(k,p)
    %===== This function generates data-bit string =====%
        %=====  based on Bernoulli distribution =====%
    data_bits_length=2^k-k-1;
    rand_seq =(rand(1,data_bits_length)<=p);
end