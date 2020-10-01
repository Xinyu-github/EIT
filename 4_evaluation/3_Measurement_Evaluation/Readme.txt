IMPORTANT: See further below for info on the indexing.

In this folder all the measurements taken for the master's thesis of Simon 
Stümke in the winter semester of 2019/2020 are stored. The folder also
includes scripts to reconstruct all measurements in a time differential as
well as frequency differential way. The scripts are stores in the 'scripts'
subfolder. In the top section of each of these scripts the parameters that
should be used for the reconstruction can be specified.

The folder '2020_01_13_aixtom_dummy','2020_03_05_swisstom', 
'2020_03_05_veggie_conductivities' and 'Older' can be found under
U:\EIT_Hardware_II\Messungen\2020_MeasurementStuemke
and should not be stored in the git.

The folder '2020_01_13_aixtom_dummy' contains measurements that were taken
with the EIT-Phantom-001 using the aixtom dummy. The measurments were taken
using the adjacent measurement pattern and 16 electrodes. The
measurements themselves are stored in the .txt-files. The measurement
specifications are stored in the .mat-file of the same name.

The folder '2020_03_05_swisstom' contains measurements that were taken on the
3rd and 5th of March 2020 using the Swisstom EIT Pioneer Set. The
measurements were taken using the mobile EIT water tank at MEDIT. The
tank was filled with a salt-water-solution with a background condictivity
of 11,2 S/m in which a potato and one half of a butternut pumpking were
submerged at different poisitions. The measurements were taken using the
adjacent and Skip-4-4 pattern using 32 electrodes. The folder 'simulation
equivalent' contains a script that creates an FE-Model of the tank in which
the pumpking and potato are submerges. This model can be used to gerenate
reference voltages to compare to the measured voltages.

The folder '2020_03_05_veggie_conductivities' constains measurements that 
were taken on the 5th of March 2020 with the Agilent E4980A precision lcr 
meter. A piece of pumpkin and a piece of potato were measured in order to
determine ther dielectric properties. The pieces measured 3,7cmX3.4cmX2,8cm.
They were measured using the 'short' and 'long' measurement mode between
10e4 Hz and 10e6 Hz with either 30 or 60 individual measurements per
measurement cycle. The result are stored in the .mat-file that can be 
found in the folder.

The folder 'Older' contains other miscellaneous older measurements. 

The folder 'Scripts' contains all scripts that are needed to reconstruct
and/or analyse the measured data.
The functions 'AixTOM_extractData.m' and 'AixTOM_readFile.m' were developed
by Tobias Menden and are taken from the AixCom program files in order to
use data measured with the AixTOM system for reconstruction.

IMPORTANT: All measurements are stored in a raw and an indexed style. This
helps to find the different vegetable positions especially for the swisstom
measurments since these were taken in a continuous measurement mode. The 
single measurements for the different vegetable positions were thus
separated and extracted from the continuous file. The resulting files
were indexed using the convention explained below:

11 - Potato
12 - Pumpkin
13 - Pumpking and Potato

20 - Homogenous Tank without vegetables
21 - Close to the tank edge at electrode 1 (oriented towards electrode 1 for 'pumpkin and potato')
22 - Close to the tank edge at electrode 5
23 - Close to the tank edge at electrode 9
24 - Close to the tank edge at electrode 13
25 - Tank middle