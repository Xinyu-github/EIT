function [Mcomplex, McomplexU, McomplexI] = AixTOM_extractData(meas_data, samples_to_read, measurements_per_frame)

meas_data_double = double(meas_data);

Mcomplex = [];
McomplexU = [];
McomplexI = [];
k=1; 
l=1; %measurement
%start_sample = meas_data_double(k+samples_to_read*2*2+0);%4
while k<(samples_to_read * 2 * 2 ) * measurements_per_frame
    %start = meas_data_double(k+samples_to_read*2*2+3)

    for j=1:1:samples_to_read
        m=k+(j-1)*4;
        %McomplexU(l,2*j) = (meas_data_double(m)+1i*meas_data_double(m+1))*meas_data_double(k+samples_to_read*2*2+2); %higher frequency
        %McomplexU(l,j*2-1) = (meas_data_double(m+2)+1i*meas_data_double(m+3))*meas_data_double(k+samples_to_read*2*2+2); %lower frequency
        
        %Data arrive in 4*16bit packages, each sample contains: xfft_meas1_tdata & xfft_meas1_I_tdata;
        % lowest bit will be transmited first: m: RE(meas_I), m+1=1:
        % IMAG(meas_I), RE(meas) m+1=1: IMAG(meas)
        McomplexU(l,j) = (meas_data_double(m+3)+1i*meas_data_double(m+2));%*meas_data_double(k+samples_to_read*2*2+2); %higher frequency

        %McomplexI(l,2*j) = (meas_data_double(m+4)+1i*meas_data_double(m+5))*meas_data_double(k+samples_to_read*2*2+2); %higher frequency
        %McomplexI(l,2*j-1) = (meas_data_double(m+6)+1i*meas_data_double(m+7))*meas_data_double(k+samples_to_read*2*2+2); %lower frequency
        McomplexI(l,j) = (meas_data_double(m+1)+1i*meas_data_double(m));%*meas_data_double(k+samples_to_read*2*2+2); %higher frequency

        %Mcomplex(l,j) = (meas_data_double(m)+1i*meas_data_double(m+1));%*M(k+6,1);
            
    end
    Mcomplex(l,:) = McomplexU(l,:)./McomplexI(l,:);    
    
%     current_sample = meas_data_double(m+4);%8
%     if current_sample-start_sample+1 ~= l
%         display 'sample missmatch';
%        break 
%     end
    l=l+1;
    k = k+samples_to_read*4;
end

end 
% if M((j*samples)+2*(j-1)+1,1)~=j
%     display('measurement missmatch')
%     return;
% end



% for k=1:5%length(Mcomplex(:,1))
% figure;
% plot(real(Mcomplex(k,:)));
% end

% % figure;
% % plot(abs(real(Mcomplex(13,:))));
% % 
% % figure;
% % plot(abs(real(Mcomplex(:,1))));
% % hold on;
% % plot(abs(real(Mcomplex(:,2))));

% figure;
% plot(abs(real(Mcomplex(3,:,1))));
% hold on;
% plot(abs(real(Mcomplex(3,:,2))));
% 
% figure;
% plot(abs(imag(Mcomplex(1,:,1))));
% hold on;
% plot(abs(imag(Mcomplex(1,:,2))));
% 
% figure;
% plot(abs(imag(Mcomplex(10,:,1))));
% hold on;
% plot(abs(imag(Mcomplex(10,:,2))));
% 

%end