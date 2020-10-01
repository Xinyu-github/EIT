function [Mcomplex] = AixTOM_readFile(filename)

fid = fopen(filename, 'r');
% if (fid == -1)
%     mag = 0;
%     phase = 0;
% else
%tline = fgetl(fid);
%fgetl(fid);
%mag = fscanf(fid, '%d', 1);
%phase = fscanf(fid, '%d', 1);

header = fscanf(fid, '%d %d %d %d %f',5); %7 in prior version
patternid = header(1);
measurements_per_frame = header(2);
%frame = header(3);
samples = header(3);
samples_to_read  = header(4);

M=fscanf(fid,'%li,',[(samples_to_read * 4) * measurements_per_frame + 8+1,inf]);
fclose(fid);

Mcomplex = [];
McomplexU = [];
McomplexI = [];

for i=1:length(M(1,:))
    meas_data = M(2:end,i);
    [Mcomplex(i,:,:), McomplexU(i,:,:), McomplexI(i,:,:)] = AixTOM_extractData(meas_data, samples_to_read, measurements_per_frame); 
end

end