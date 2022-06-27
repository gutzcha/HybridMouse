function data_b = get_power_vector(data,frequency_range)


if ~exist('range','var')||isempty(frequency_range)
    frequency_range = [1, 125];
end

range_size= diff(frequency_range)+1;

data_b=zeros(1,range_size);
if isempty(data)
    return
end
if iscell(data)
    data = data{:};
end

%determine the size of the signal array
L= size(data,1);

NFFT = 2^nextpow2(L(1,1));

%Get fft magbitude
data_f = fft(data,NFFT)/L(1,1);
data_s=smooth(2*abs(data_f(1:NFFT/2+1)));

binsum=size(data_s)/range_size;

for j=frequency_range(1):frequency_range(end)
    data_b(j)= sum(data_s((j-1)*floor(binsum(1,1))+1:j*floor(binsum(1,1))));
end

end
