function tissues = my_getTissues(freq)
%GETTISSUES tissues = gettissues(freq)
%�bergibt alle verf�gbaren Tissuewerte f�r eine bestimmte frequenz in einer
%Struktur

% Formel: http://en.wikipedia.org/wiki/Permittivity#Lossy_medium

% Modified by Sebastien Dambrun to give the combined conductivities of
% somme tissues that were not already given.

epsilon0_const = 8.85418781762*10^(-12);

tissues.Aorta_epsilonr = real(4+(40/(1+(1i*2*pi*freq*0.000000000008842)^(1-0.1)))+(50/(1+(1i*2*pi*freq*0.000000003183)^(1-0.1)))+(100000/(1+(1i*2*pi*freq*0.000159155)^(1-0.2)))+(10000000/(1+(1i*2*pi*freq*0.001592)^(1-0)))+(0.25/(1i*2*pi*freq*epsilon0_const)));
tissues.Aorta_sigma = -2*pi*freq*epsilon0_const*imag(4+(40/(1+(1i*2*pi*freq*0.000000000008842)^(1-0.1)))+(50/(1+(1i*2*pi*freq*0.000000003183)^(1-0.1)))+(100000/(1+(1i*2*pi*freq*0.000159155)^(1-0.2)))+(10000000/(1+(1i*2*pi*freq*0.001592)^(1-0)))+(0.25/(1i*2*pi*freq*epsilon0_const)));

tissues.Bladder_epsilonr = real(2.5+(16/(1+(1i*2*pi*freq*0.000000000008842)^(1-0.1)))+(400/(1+(1i*2*pi*freq*0.000000159155)^(1-0.1)))+(100000/(1+(1i*2*pi*freq*0.000159155)^(1-0.2)))+(10000000/(1+(1i*2*pi*freq*0.015915)^(1-0)))+(0.2/(1i*2*pi*freq*epsilon0_const)));
tissues.Bladder_sigma = -2*pi*freq*epsilon0_const*imag(2.5+(16/(1+(1i*2*pi*freq*0.000000000008842)^(1-0.1)))+(400/(1+(1i*2*pi*freq*0.000000159155)^(1-0.1)))+(100000/(1+(1i*2*pi*freq*0.000159155)^(1-0.2)))+(10000000/(1+(1i*2*pi*freq*0.015915)^(1-0)))+(0.2/(1i*2*pi*freq*epsilon0_const)));

tissues.Blood_epsilonr = real(4+(56/(1+(1i*2*pi*freq*0.000000000008377)^(1-0.1)))+(5200/(1+(1i*2*pi*freq*0.000000132629)^(1-0.1)))+(0/(1+(1i*2*pi*freq*0.000159155)^(1-0.2)))+(0/(1+(1i*2*pi*freq*0.015915)^(1-0)))+(0.7/(1i*2*pi*freq*epsilon0_const)));
tissues.Blood_sigma = -2*pi*freq*epsilon0_const*imag(4+(56/(1+(1i*2*pi*freq*0.000000000008377)^(1-0.1)))+(5200/(1+(1i*2*pi*freq*0.000000132629)^(1-0.1)))+(0/(1+(1i*2*pi*freq*0.000159155)^(1-0.2)))+(0/(1+(1i*2*pi*freq*0.015915)^(1-0)))+(0.7/(1i*2*pi*freq*epsilon0_const)));
tissues.Blood_combined = tissues.Blood_sigma + 1i*tissues.Blood_epsilonr*(2*pi*freq*epsilon0_const);
tissues.Blood_epsilon = tissues.Blood_epsilonr - 1i*tissues.Blood_sigma/(2*pi*freq*epsilon0_const);

tissues.BoneSchwamm_epsilonr = real(2.5+(18/(1+(1i*2*pi*freq*0.000000000013263)^(1-0.22)))+(300/(1+(1i*2*pi*freq*0.000000079577)^(1-0.25)))+(20000/(1+(1i*2*pi*freq*0.000159155)^(1-0.2)))+(20000000/(1+(1i*2*pi*freq*0.015915)^(1-0)))+(0.07/(1i*2*pi*freq*epsilon0_const)));
tissues.BoneSchwamm_sigma = -2*pi*freq*epsilon0_const*imag(2.5+(18/(1+(1i*2*pi*freq*0.000000000013263)^(1-0.22)))+(300/(1+(1i*2*pi*freq*0.000000079577)^(1-0.25)))+(20000/(1+(1i*2*pi*freq*0.000159155)^(1-0.2)))+(20000000/(1+(1i*2*pi*freq*0.015915)^(1-0)))+(0.07/(1i*2*pi*freq*epsilon0_const)));
tissues.Bone_Schwamm_combined = tissues.BoneSchwamm_sigma + 1i*tissues.BoneSchwamm_epsilonr*(2*pi*freq*epsilon0_const);
tissues.Blood_epsilon = tissues.BoneSchwamm_epsilonr - 1i*tissues.BoneSchwamm_sigma/(2*pi*freq*epsilon0_const);

tissues.BoneAussen_epsilonr = real(2.5+(10/(1+(1i*2*pi*freq*0.000000000013263)^(1-0.2)))+(180/(1+(1i*2*pi*freq*0.000000079577)^(1-0.2)))+(5000/(1+(1i*2*pi*freq*0.000159155)^(1-0.2)))+(100000/(1+(1i*2*pi*freq*0.015915)^(1-0)))+(0.02/(1i*2*pi*freq*epsilon0_const)));
tissues.BoneAussen_sigma = -2*pi*freq*epsilon0_const*imag(2.5+(10/(1+(1i*2*pi*freq*0.000000000013263)^(1-0.2)))+(180/(1+(1i*2*pi*freq*0.000000079577)^(1-0.2)))+(5000/(1+(1i*2*pi*freq*0.000159155)^(1-0.2)))+(100000/(1+(1i*2*pi*freq*0.015915)^(1-0)))+(0.02/(1i*2*pi*freq*epsilon0_const)));
tissues.BoneAussen_combined = tissues.BoneAussen_sigma + 1i*tissues.BoneAussen_epsilonr*(2*pi*freq*epsilon0_const);
tissues.BoneAussen_epsilon = tissues.BoneAussen_epsilonr - 1i*tissues.BoneAussen_sigma/(2*pi*freq*epsilon0_const);

tissues.BoneMarkInfiltrated_epsilonr = real(2.5+(9/(1+(1i*2*pi*freq*0.000000000014469)^(1-0.2)))+(80/(1+(1i*2*pi*freq*0.000000015915)^(1-0.1)))+(10000/(1+(1i*2*pi*freq*0.001591549)^(1-0.1)))+(2000000/(1+(1i*2*pi*freq*0.015915)^(1-0.1)))+(0.1/(1i*2*pi*freq*epsilon0_const)));
tissues.BoneMarkInfiltrated_sigma = -2*pi*freq*epsilon0_const*imag(2.5+(9/(1+(1i*2*pi*freq*0.000000000014469)^(1-0.2)))+(80/(1+(1i*2*pi*freq*0.000000015915)^(1-0.1)))+(10000/(1+(1i*2*pi*freq*0.001591549)^(1-0.1)))+(2000000/(1+(1i*2*pi*freq*0.015915)^(1-0.1)))+(0.1/(1i*2*pi*freq*epsilon0_const)));

tissues.BoneMarkNotInfiltrated_epsilonr = real(2.5+(3/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.2)))+(25/(1+(1i*2*pi*freq*0.000000015915)^(1-0.1)))+(5000/(1+(1i*2*pi*freq*0.001591549)^(1-0.1)))+(2000000/(1+(1i*2*pi*freq*0.015915)^(1-0.1)))+(0.001/(1i*2*pi*freq*epsilon0_const)));
tissues.BoneMarkNotInfiltrated_sigma = -2*pi*freq*epsilon0_const*imag(2.5+(3/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.2)))+(25/(1+(1i*2*pi*freq*0.000000015915)^(1-0.1)))+(5000/(1+(1i*2*pi*freq*0.001591549)^(1-0.1)))+(2000000/(1+(1i*2*pi*freq*0.015915)^(1-0.1)))+(0.001/(1i*2*pi*freq*epsilon0_const)));

tissues.Colon_epsilonr = real(4+(50/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.1)))+(3000/(1+(1i*2*pi*freq*0.000000159155)^(1-0.2)))+(100000/(1+(1i*2*pi*freq*0.000159155)^(1-0.2)))+(40000000/(1+(1i*2*pi*freq*0.001592)^(1-0)))+(0.01/(1i*2*pi*freq*epsilon0_const)));
tissues.Colon_sigma = -2*pi*freq*epsilon0_const*imag(4+(50/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.1)))+(3000/(1+(1i*2*pi*freq*0.000000159155)^(1-0.2)))+(100000/(1+(1i*2*pi*freq*0.000159155)^(1-0.2)))+(40000000/(1+(1i*2*pi*freq*0.001592)^(1-0)))+(0.01/(1i*2*pi*freq*epsilon0_const)));

tissues.Dunndarm_epsilonr = real(4+(50/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.1)))+(10000/(1+(1i*2*pi*freq*0.000000159155)^(1-0.1)))+(500000/(1+(1i*2*pi*freq*0.000159155)^(1-0.2)))+(40000000/(1+(1i*2*pi*freq*0.015915)^(1-0)))+(0.5/(1i*2*pi*freq*epsilon0_const)));
tissues.Dunndarm_sigma = -2*pi*freq*epsilon0_const*imag(4+(50/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.1)))+(10000/(1+(1i*2*pi*freq*0.000000159155)^(1-0.1)))+(500000/(1+(1i*2*pi*freq*0.000159155)^(1-0.2)))+(40000000/(1+(1i*2*pi*freq*0.015915)^(1-0)))+(0.5/(1i*2*pi*freq*epsilon0_const)));

tissues.FatAverageInfiltrated_epsilonr = real(2.5+(9/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.2)))+(35/(1+(1i*2*pi*freq*0.000000015915)^(1-0.1)))+(33000/(1+(1i*2*pi*freq*0.000159155)^(1-0.05)))+(10000000/(1+(1i*2*pi*freq*0.015915)^(1-0.01)))+(0.035/(1i*2*pi*freq*epsilon0_const)));
tissues.FatAverageInfiltrated_sigma = -2*pi*freq*epsilon0_const*imag(2.5+(9/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.2)))+(35/(1+(1i*2*pi*freq*0.000000015915)^(1-0.1)))+(33000/(1+(1i*2*pi*freq*0.000159155)^(1-0.05)))+(10000000/(1+(1i*2*pi*freq*0.015915)^(1-0.01)))+(0.035/(1i*2*pi*freq*epsilon0_const)));

tissues.FatNotInfiltrated_epsilonr = real(2.5+(3/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.2)))+(15/(1+(1i*2*pi*freq*0.000000015915)^(1-0.1)))+(33000/(1+(1i*2*pi*freq*0.000159155)^(1-0.05)))+(10000000/(1+(1i*2*pi*freq*0.007958)^(1-0.01)))+(0.01/(1i*2*pi*freq*epsilon0_const)));
tissues.FatNotInfiltrated_sigma = -2*pi*freq*epsilon0_const*imag(2.5+(3/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.2)))+(15/(1+(1i*2*pi*freq*0.000000015915)^(1-0.1)))+(33000/(1+(1i*2*pi*freq*0.000159155)^(1-0.05)))+(10000000/(1+(1i*2*pi*freq*0.007958)^(1-0.01)))+(0.01/(1i*2*pi*freq*epsilon0_const)));

tissues.Heart_epsilonr = real(4+(50/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.1)))+(1200/(1+(1i*2*pi*freq*0.0000001592)^(1-0.05)))+(450000/(1+(1i*2*pi*freq*0.00007234)^(1-0.22)))+(25000000/(1+(1i*2*pi*freq*0.004547)^(1-0)))+(0.05/(1i*2*pi*freq*epsilon0_const)));
tissues.Heart_sigma = -2*pi*freq*epsilon0_const*imag(4+(50/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.1)))+(1200/(1+(1i*2*pi*freq*0.0000001592)^(1-0.05)))+(450000/(1+(1i*2*pi*freq*0.00007234)^(1-0.22)))+(25000000/(1+(1i*2*pi*freq*0.004547)^(1-0)))+(0.05/(1i*2*pi*freq*epsilon0_const)));
tissues.Heart_combined = tissues.Heart_sigma + 1i*tissues.Heart_epsilonr*(2*pi*freq*epsilon0_const);
tissues.Heart_epsilon = tissues.Heart_epsilonr - 1i*tissues.Heart_sigma/(2*pi*freq*epsilon0_const);

tissues.Knorpel_epsilonr = real(4+(38/(1+(1i*2*pi*freq*0.000000000013263)^(1-0.15)))+(2500/(1+(1i*2*pi*freq*0.000000144686)^(1-0.15)))+(100000/(1+(1i*2*pi*freq*0.00031831)^(1-0.1)))+(40000000/(1+(1i*2*pi*freq*0.015915)^(1-0)))+(0.15/(1i*2*pi*freq*epsilon0_const)));
tissues.Knorpel_sigma = -2*pi*freq*epsilon0_const*imag(4+(38/(1+(1i*2*pi*freq*0.000000000013263)^(1-0.15)))+(2500/(1+(1i*2*pi*freq*0.000000144686)^(1-0.15)))+(100000/(1+(1i*2*pi*freq*0.00031831)^(1-0.1)))+(40000000/(1+(1i*2*pi*freq*0.015915)^(1-0)))+(0.15/(1i*2*pi*freq*epsilon0_const)));

tissues.LungDeflated_epsilonr = real(4+(45/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.1)))+(1000/(1+(1i*2*pi*freq*0.0000001592)^(1-0.1)))+(500000/(1+(1i*2*pi*freq*0.0001592)^(1-0.2)))+(10000000/(1+(1i*2*pi*freq*0.01592)^(1-0)))+(0.2/(1i*2*pi*freq*epsilon0_const)));
tissues.LungDeflated_sigma = -2*pi*freq*epsilon0_const*imag(4+(45/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.1)))+(1000/(1+(1i*2*pi*freq*0.0000001592)^(1-0.1)))+(500000/(1+(1i*2*pi*freq*0.0001592)^(1-0.2)))+(10000000/(1+(1i*2*pi*freq*0.01592)^(1-0)))+(0.2/(1i*2*pi*freq*epsilon0_const)));
tissues.LungDeflated_combined = tissues.LungDeflated_sigma + 1i*tissues.LungDeflated_epsilonr*(2*pi*freq*epsilon0_const);
tissues.LungDeflated_epsilon = tissues.LungDeflated_epsilonr - 1i*tissues.LungDeflated_sigma/(2*pi*freq*epsilon0_const);

tissues.LungInflated_epsilonr = real(2.5+(18/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.1)))+(500/(1+(2*pi*1i*freq*0.00000006366)^(1-0.1)))+(250000/(1+(2*pi*1i*freq*0.0001592)^(1-0.2)))+(40000000/(1+(2*pi*1i*freq*0.007958)^(1-0)))+(0.03/(1i*2*pi*freq*epsilon0_const)));
tissues.LungInflated_sigma = -2*pi*freq*epsilon0_const*imag(2.5+(18/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.1)))+(500/(1+(2*pi*1i*freq*0.00000006366)^(1-0.1)))+(250000/(1+(2*pi*1i*freq*0.0001592)^(1-0.2)))+(40000000/(1+(2*pi*1i*freq*0.007958)^(1-0)))+(0.03/(1i*2*pi*freq*epsilon0_const)));
tissues.LungInflated_combined = tissues.LungInflated_sigma + 1i*tissues.LungInflated_epsilonr*(2*pi*freq*epsilon0_const);
tissues.LungInflated_epsilon = tissues.LungInflated_epsilonr-1i*tissues.LungInflated_sigma/(2*pi*freq*epsilon0_const);

tissues.Muscle_epsilonr = real(4+(50/(1+(1i*2*pi*freq*0.000000000007234)^(1-0.1)))+(7000/(1+(1i*2*pi*freq*0.000000353678)^(1-0.1)))+(1200000/(1+(1i*2*pi*freq*0.00031831)^(1-0.1)))+(25000000/(1+(1i*2*pi*freq*0.002274)^(1-0)))+(0.2/(1i*2*pi*freq*epsilon0_const)));
tissues.Muscle_sigma = -2*pi*freq*epsilon0_const*imag(4+(50/(1+(1i*2*pi*freq*0.000000000007234)^(1-0.1)))+(7000/(1+(1i*2*pi*freq*0.000000353678)^(1-0.1)))+(1200000/(1+(1i*2*pi*freq*0.00031831)^(1-0.1)))+(25000000/(1+(1i*2*pi*freq*0.002274)^(1-0)))+(0.2/(1i*2*pi*freq*epsilon0_const)));
tissues.Muscle_combined = tissues.Muscle_sigma + 1i*tissues.Muscle_epsilonr*(2*pi*freq*epsilon0_const);
tissues.Muscle_epsilon = tissues.Muscle_epsilonr*epsilon0_const-1i*tissues.Muscle_sigma/(2*pi*freq) ;
tissues.Muscle_sigmaTest = -2*pi*freq*epsilon0_const*(4+(50/(1+(1i*2*pi*freq*0.000000000007234)^(1-0.1)))+(7000/(1+(1i*2*pi*freq*0.000000353678)^(1-0.1)))+(1200000/(1+(1i*2*pi*freq*0.00031831)^(1-0.1)))+(25000000/(1+(1i*2*pi*freq*0.002274)^(1-0)))+(0.2/(1i*2*pi*freq*epsilon0_const)));

tissues.SkinDry_epsilonr = real(4+(32/(1+(1i*2*pi*freq*0.000000000007234)^(1-0)))+(1100/(1+(1i*2*pi*freq*0.000000032481)^(1-0.2)))+(0/(1+(1i*2*pi*freq*0.000159155)^(1-0.2)))+(0/(1+(1i*2*pi*freq*0.015915)^(1-0.2)))+(0.0002/(1i*2*pi*freq*epsilon0_const)));
tissues.SkinDry_sigma = -2*pi*freq*epsilon0_const*imag(4+(32/(1+(1i*2*pi*freq*0.000000000007234)^(1-0)))+(1100/(1+(1i*2*pi*freq*0.000000032481)^(1-0.2)))+(0/(1+(1i*2*pi*freq*0.000159155)^(1-0.2)))+(0/(1+(1i*2*pi*freq*0.015915)^(1-0.2)))+(0.0002/(1i*2*pi*freq*epsilon0_const)));

tissues.SkinWet_epsilonr = real(4+(39/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.1)))+(280/(1+(1i*2*pi*freq*0.000000079577)^(1-0)))+(30000/(1+(1i*2*pi*freq*0.000001592)^(1-0.16)))+(30000/(1+(1i*2*pi*freq*0.001592)^(1-0.2)))+(0.0004/(1i*2*pi*freq*epsilon0_const)));
tissues.SkinWet_sigma = -2*pi*freq*epsilon0_const*imag(4+(39/(1+(1i*2*pi*freq*0.000000000007958)^(1-0.1)))+(280/(1+(1i*2*pi*freq*0.000000079577)^(1-0)))+(30000/(1+(1i*2*pi*freq*0.000001592)^(1-0.16)))+(30000/(1+(1i*2*pi*freq*0.001592)^(1-0.2)))+(0.0004/(1i*2*pi*freq*epsilon0_const)));

tissues.Average_epsilonr = real(3.63636363636364+(39.0681818181818/(1+(1i*2*pi*freq*8.8294090909091E-12)^(1-0.110681818181818)))+(1634.09090909091/(1+(1i*2*pi*freq*9.32335227272727E-08)^(1-0.115227272727273)))+(135636.363636364/(1+(1i*2*pi*freq*0.000214730227272727)^(1-0.176363636363636)))+(21844204.5454545/(1+(1i*2*pi*freq*0.0118148409090909)^(1-0.0202272727272727)))+(0.305045454545454/(1i*2*pi*freq*epsilon0_const)));
tissues.Average_sigma = -2*pi*freq*epsilon0_const*imag(3.63636363636364+(39.0681818181818/(1+(1i*2*pi*freq*8.8294090909091E-12)^(1-0.110681818181818)))+(1634.09090909091/(1+(1i*2*pi*freq*9.32335227272727E-08)^(1-0.115227272727273)))+(135636.363636364/(1+(1i*2*pi*freq*0.000214730227272727)^(1-0.176363636363636)))+(21844204.5454545/(1+(1i*2*pi*freq*0.0118148409090909)^(1-0.0202272727272727)))+(0.305045454545454/(1i*2*pi*freq*epsilon0_const)));
tissues.Average_combined = tissues.Average_sigma + 1i*2*pi*freq*epsilon0_const*tissues.Average_epsilonr;
tissues.Average_epsilon = tissues.Average_epsilonr - 1i*tissues.Average_sigma/(2*pi*freq*epsilon0_const); 
tissues.Average_sigma2 = -tissues.Average_sigma/(2*pi*freq*epsilon0_const);

if freq == 100e3
tissues.fluid_epsilonr = 97.99;
tissues.fluid_sigma = 1.5;
tissues.fluid_combined = tissues.fluid_sigma + 1i*tissues.fluid_epsilonr*(2*pi*freq*epsilon0_const);
end

if freq == 700e3
tissues.fluid_epsilonr = 88.05;
tissues.fluid_sigma = 1.5;
tissues.fluid_combined = tissues.fluid_sigma + 1i*tissues.fluid_epsilonr*(2*pi*freq*epsilon0_const);
end

if freq ~= 100e3 || freq ~= 700e3
    
   disp('this provides fluid conductivity just for 100e3,700e3 Hz'); 
end
