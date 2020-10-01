
%n = [10,10,10,10,10,10,11,11,11,12,12];
TYPE = {'withbeta','diff'};
% index of best hp for various position in mode1 
% selected by L-curve
n = [7,7,4,4,5,6,6,6,11,10,9]; 

% index of best hp for various position in mode2
n1 = [10,10,10,10,14,12,9,8,12,19,14];
position = 0:-0.1:-1;
for i = 1:11
%     [res,reg]=generate_lcurve(img_beta{i});
    for j = 1:length(TYPE)
        load(strcat('circle_',TYPE{j},'_case',num2str(i),'.mat'));
    end
    % figure of merit for mode1
    FOM1(i).param= eval_GREIT_fig_merit_fdEIT(img_beta{1}(n(i)), [position(i);0;0.15]);
    IA1(i) =FOM1(i).param(1,:); %image amplitude
    PE1(i) = FOM1(i).param(2,:);%position error
    RES1(i) = FOM1(i).param(3,:);% resolution
    SD1(i) = FOM1(i).param(4,:); % shape deformation
    
    img = beta_selection(img_beta{1}(n(i)),img_mode2(n1(i)),0);% combine mode1 and mode2
    % figure of merit for combined image
    FOM2(i).param= eval_GREIT_fig_merit_fdEIT(img, [position(i);0;0.15]);
    IA2(i) =FOM2(i).param(1,:);
    PE2(i) = FOM2(i).param(2,:);
    RES2(i) = FOM2(i).param(3,:);
    SD2(i) = FOM2(i).param(4,:);

end

figure;
subplot(4,1,1)
plot(abs(position),real(IA1))
hold on
plot(abs(position),real(IA2))

legend('cost func1','cost func2')    
title('Image Amplitude')
%  params(1,:) = Image Amplitude
%  params(2,:) = Position Error => + toward centre, - toward edge
%  params(3,:) = Resolution
%  params(4,:) = Shape Deformation

subplot(4,1,2)
plot(abs(position),PE1)
hold on
plot(abs(position),PE2)
legend('cost func1','cost func2')   
title('Position Error')

subplot(4,1,3)
plot(abs(position),RES1)
hold on
plot(abs(position),RES2)
legend('cost func1','cost func2')   
title('Resolution')

subplot(4,1,4)
plot(abs(position),SD1)
hold on
plot(abs(position),SD2)
legend('cost func1','cost func2') 
title('Shape Deformation')
