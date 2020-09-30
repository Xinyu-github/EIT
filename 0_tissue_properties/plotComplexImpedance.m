% 1kHZ - 1MHz
%Heart_complex = Heart_cond(41:101) + 1i*eps0.*Freq(41:101)*2*pi.*Heart_perm(41:101);
%Heart_complex = Heart_cond + 1i*eps0.*Freq*2*pi.*Heart_perm;

% Freq(55) = 5kHz
% 61 = 100kHz
% 98 = 700kHz
%Heart_complex = Heart_cond + 1i*eps0.*Freq*2*pi.*Heart_perm;
clear

load('All_f_set_1.mat');

clear title xlabel ylabel legend
figure
hold all;

load('All_f_set_1.mat');


plot(real(Heart_pho(41:101)), imag(Heart_pho(41:101)), '-.g',...
    'LineWidth',2)
plot(real(Muscle_pho(41:101)), imag(Muscle_pho(41:101)), '-.b',...
    'LineWidth',2)
plot(real(LungDef_pho(41:101)), imag(LungDef_pho(41:101)), '-.m',...
    'LineWidth',2)
plot(real(LungInf_pho(41:101)), imag(LungInf_pho(41:101)), '-.k',...
    'LineWidth',2)

plot(real(Heart_pho(61)), imag(Heart_pho(61)), 'o',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',6)%'MarkerFaceColor',[.49 1 .63],...
plot(real(Heart_pho(98)), imag(Heart_pho(98)), 'x',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',8)%'MarkerFaceColor',[.49 1 .63],...

plot(real(Muscle_pho(61)), imag(Muscle_pho(61)), 'o',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',6)%'MarkerFaceColor',[.49 1 .63],...
plot(real(Muscle_pho(98)), imag(Muscle_pho(98)), 'x',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',8)%'MarkerFaceColor',[.49 1 .63],...

plot(real(LungDef_pho(61)), imag(LungDef_pho(61)), 'o',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',6)%'MarkerFaceColor',[.49 1 .63],...
plot(real(LungDef_pho(98)), imag(LungDef_pho(98)), 'x',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',8)%'MarkerFaceColor',[.49 1 .63],...

plot(real(LungInf_pho(61)), imag(LungInf_pho(61)), 'o',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',6)%'MarkerFaceColor',[.49 1 .63],...
plot(real(LungInf_pho(98)), imag(LungInf_pho(98)), 'x',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',8)%'MarkerFaceColor',[.49 1 .63],...



p1 = [real(Heart_pho(61)) imag(Heart_pho(61))];                      
p2 = [real(Heart_pho(98)) imag(Heart_pho(98))];
arrow(p1, p2,'LineWidth',1.5,'Color','g','length',10);
%dp = p2-p1;                         % Difference
%quiver(p1(1),p1(2),dp(1),dp(2),0,'MaxHeadSize',0.5,'LineWidth',2)

p1 = [real(Muscle_pho(61)) imag(Muscle_pho(61))];                      
p2 = [real(Muscle_pho(98)) imag(Muscle_pho(98))];
arrow(p1, p2,'LineWidth',1.5,'Color','b','length',10);

p1 = [0 0];                      
p2 = [real(Muscle_pho(61)) imag(Muscle_pho(61))]; 
arrow(p1, p2,'LineWidth',1.5,'Color','b','length',10,'LineStyle',':');
p1 = [0 0];                
p2 = [real(Muscle_pho(98)) imag(Muscle_pho(98))];
arrow(p1, p2,'LineWidth',1.5,'Color','b','length',10,'LineStyle',':');


p1 = [real(LungDef_pho(61)) imag(LungDef_pho(61))];                      
p2 = [real(LungDef_pho(98)) imag(LungDef_pho(98))];
arrow(p1, p2,'LineWidth',1.5,'Color','m','length',10);

p1 = [real(LungInf_pho(61)) imag(LungInf_pho(61))];                      
p2 = [real(LungInf_pho(98)) imag(LungInf_pho(98))];
arrow(p1, p2,'LineWidth',1.5,'Color','k','length',10);

grid
%axis([0  10    0  10])
%text(p1(1),p1(2), sprintf('(%.0f,%.0f)',p1))
%text(p2(1),p2(2), sprintf('(%.0f,%.0f)',p2))
xlabel('real-component [\Omega/m]');
ylabel('imag-component [\Omega/m]');

lgd = legend('Heart', 'Muscle', 'Lung deflated', 'Lung inflated', '100 kHz', '700 kHz');
set(lgd,'Interpreter','latex');
%set(lgd,'FontSize',17);

%matlab2tikz('tikz/tissueProperties.tikz');

figure
p1 = [real(Heart_pho(61)) imag(Heart_pho(61))];                      
p2 = [real(Heart_pho(98)) imag(Heart_pho(98))];
arrow([0 0], p2-p1,'LineWidth',1.5,'Color','g','length',10);
%dp = p2-p1;                         % Difference
%quiver(p1(1),p1(2),dp(1),dp(2),0,'MaxHeadSize',0.5,'LineWidth',2)

p1 = [real(Muscle_pho(61)) imag(Muscle_pho(61))];                      
p2 = [real(Muscle_pho(98)) imag(Muscle_pho(98))];
arrow([0 0], p2-p1,'LineWidth',1.5,'Color','b','length',10);

p1 = [real(LungDef_pho(61)) imag(LungDef_pho(61))];                      
p2 = [real(LungDef_pho(98)) imag(LungDef_pho(98))];
arrow([0 0], p2-p1,'LineWidth',1.5,'Color','m','length',10);

p1 = [real(LungInf_pho(61)) imag(LungInf_pho(61))];                      
p2 = [real(LungInf_pho(98)) imag(LungInf_pho(98))];
arrow([0 0], p2-p1,'LineWidth',1.5,'Color','k','length',10);

grid
