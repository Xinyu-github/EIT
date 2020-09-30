clear
load('All_f_set_1.mat');
% 1kHZ - 1MHz
%Heart_complex = Heart_cond(41:101) + 1i*eps0.*Freq(41:101)*2*pi.*Heart_perm(41:101);
%Heart_complex = Heart_cond + 1i*eps0.*Freq*2*pi.*Heart_perm;
Heart_complex = Heart_cond + 1i*eps0.*Freq*2*pi.*Heart_perm;
Muscle_complex = Muscle_cond + 1i*eps0.*Freq*2*pi.*Muscle_perm;
LungDef_complex = LungDef_cond + 1i*eps0.*Freq*2*pi.*LungDef_perm;
BoneCortical_complex = BoneCortical_cond + 1i*eps0.*Freq*2*pi.*BoneCortical_perm;
LungInf_complex = LungInf_cond + 1i*eps0.*Freq*2*pi.*LungInf_perm;
%BodyFluid_complex = BodyFluid_cond + 1i*eps0.*Freq*2*pi.*BodyFluid_perm;

% Freq(55) = 5kHz
% 61 = 100kHz
% 98 = 700kHz
%Heart_complex = Heart_cond + 1i*eps0.*Freq*2*pi.*Heart_perm;

for i=1:length(Freq)
    all_tissues = getTissues_mod(Freq(i));
    Heart_complex_gabriel(i,1) = all_tissues.Heart_combined;
    Muscle_complex_gabriel(i,1) = all_tissues.Muscle_combined;
    LungDef_complex_gabriel(i,1) = all_tissues.LungDeflated_combined;
    LungInf_complex_gabriel(i,1) = all_tissues.LungInflated_combined;
    BoneCortical_complex_gabriel(i,1) = all_tissues.BoneAussen_combined;
end
    Heart_complex_gabriel = Heart_complex_gabriel;
    Muscle_complex= Muscle_complex_gabriel;
    LungDef_complex= LungDef_complex_gabriel;
    LungInf_complex= LungInf_complex_gabriel;
    BoneCortical_complex= BoneCortical_complex_gabriel;


figure
clear title xlabel ylabel legend
set(gcf, 'Position', [100, 100, 1440, 700],'Color','white');
hold all;

plot(real(Heart_complex(41:101)), 100*imag(Heart_complex(41:101)), '-.g',...
    'LineWidth',2)
plot(real(Muscle_complex(41:101)), 100*imag(Muscle_complex(41:101)), '-.b',...
    'LineWidth',2)
plot(real(LungDef_complex(41:101)), 100*imag(LungDef_complex(41:101)), '-.m',...
    'LineWidth',2)
plot(real(LungInf_complex(41:101)), 100*imag(LungInf_complex(41:101)), '-.k',...
    'LineWidth',2)
plot(real(BoneCortical_complex(41:101)), 100*imag(BoneCortical_complex(41:101)), '-.c',...
    'LineWidth',2)
%plot(real(BodyFluid_complex(41:101)), 100*imag(BodyFluid_complex(41:101)), '-.m',...
%    'LineWidth',2)


plot(real(Heart_complex(81)), 100*imag(Heart_complex(81)), 'o',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',6)%'MarkerFaceColor',[.49 1 .63],...
plot(real(Heart_complex(98)), 100*imag(Heart_complex(98)), 'x',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',8)%'MarkerFaceColor',[.49 1 .63],...

plot(real(Muscle_complex(81)), 100*imag(Muscle_complex(81)), 'o',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',6)%'MarkerFaceColor',[.49 1 .63],...
plot(real(Muscle_complex(98)), 100*imag(Muscle_complex(98)), 'x',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',8)%'MarkerFaceColor',[.49 1 .63],...

plot(real(LungDef_complex(81)), 100*imag(LungDef_complex(81)), 'o',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',6)%'MarkerFaceColor',[.49 1 .63],...
plot(real(LungDef_complex(98)), 100*imag(LungDef_complex(98)), 'x',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',8)%'MarkerFaceColor',[.49 1 .63],...

plot(real(LungInf_complex(81)), 100*imag(LungInf_complex(81)), 'o',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',6)%'MarkerFaceColor',[.49 1 .63],...
plot(real(LungInf_complex(98)), 100*imag(LungInf_complex(98)), 'x',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',8)%'MarkerFaceColor',[.49 1 .63],...

plot(real(BoneCortical_complex(81)), 100*imag(BoneCortical_complex(81)), 'o',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',6)%'MarkerFaceColor',[.49 1 .63],...
plot(real(BoneCortical_complex(98)), 100*imag(BoneCortical_complex(98)), 'x',...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerSize',8)%'MarkerFaceColor',[.49 1 .63],...

p1 = [real(Heart_complex(81)) 100*imag(Heart_complex(81))];                      
p2 = [real(Heart_complex(98)) 100*imag(Heart_complex(98))];
arrow(p1, p2,'LineWidth',1.5,'Color','g','length',10);
%dp = p2-p1;                         % Difference
%quiver(p1(1),p1(2),dp(1),dp(2),0,'MaxHeadSize',0.5,'LineWidth',2)

p1 = [real(Muscle_complex(81)) 100*imag(Muscle_complex(81))];                      
p2 = [real(Muscle_complex(98)) 100*imag(Muscle_complex(98))];
arrow(p1, p2,'LineWidth',1.5,'Color','b','length',10);

% p1 = [0 0];                      
% p2 = [real(Muscle_complex(61)) imag(Muscle_complex(61))]; 
% arrow(p1, p2,'LineWidth',1.5,'Color','b','length',10,'LineStyle',':');
% p1 = [0 0];                
% p2 = [real(Muscle_complex(98)) imag(Muscle_complex(98))];
% arrow(p1, p2,'LineWidth',1.5,'Color','b','length',10,'LineStyle',':');


p1 = [real(LungDef_complex(81)) 100*imag(LungDef_complex(81))];                      
p2 = [real(LungDef_complex(98)) 100*imag(LungDef_complex(98))];
arrow(p1, p2,'LineWidth',1.5,'Color','m','length',10);


%p1 = [0 0];                      
%p2 = [real(LungInf_complex(61)) 100*imag(LungInf_complex(61))]; 
%arrow(p1, p2,'LineWidth',1.5,'Color','k','length',10,'LineStyle',':');
%p1 = [0 0];                
%p2 = [real(LungInf_complex(98)) 100*imag(LungInf_complex(98))];
%arrow(p1, p2,'LineWidth',1.5,'Color','k','length',10,'LineStyle',':');

p1 = [real(LungInf_complex(81)) 100*imag(LungInf_complex(81))];                      
p2 = [real(LungInf_complex(98)) 100*imag(LungInf_complex(98))];
arrow(p1, p2,'LineWidth',1.5,'Color','k','length',10);

p1 = [real(BoneCortical_complex(81)) 100*imag(BoneCortical_complex(81))];                      
p2 = [real(BoneCortical_complex(98)) 100*imag(BoneCortical_complex(98))];
arrow(p1, p2,'LineWidth',1.5,'Color','c','length',10);


grid
%axis([0  10    0  10])
%text(p1(1),p1(2), sprintf('(%.0f,%.0f)',p1))
%text(p2(1),p2(2), sprintf('(%.0f,%.0f)',p2))
axis([0 0.52 0 12]);
xlabel('$\Re{\{\gamma(\omega)\}}$ [S\,m]','Interpreter','latex');
ylabel('$\Im{\{\gamma(\omega)\}}\cdot 10^{-2}$ [S\,m]','Interpreter','latex');

lgd = legend('Heart', 'Muscle', 'Lung deflated', 'Lung inflated', 'Bone cortical', '100 kHz', '700 kHz');
set(lgd,'Interpreter','latex','location','northwest');
fSize=18;
set(lgd,'FontSize',fSize);
set(gca,'FontSize',fSize)

matlab2tikz('tikz/tissueProperties_conductivity.tikz');

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
