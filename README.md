# SplitbeltMechWork
Calculating the mechanical work performed during split-belt walking.

SplitBelt_ILM.m is the main code that estimates the mechanical work done when walking on a split-belt treadmill. It uses an extension of the Individual Limbs Method described in Donelan and Kuo 2002. The folder TestCodeData contains data collected from two participants in our study Sanchaz et. al. 2019 (https://physoc.onlinelibrary.wiley.com/doi/abs/10.1113/JP277725) that serves as sample data on which to run this code. The code is set up to loop through all the participants in the folder. This repository also contains three functions that can be used to debug or verify that the outputs make sense. All code are commented to explain the calculations, and the input and output variables. But there are a number of variables related to work whose names may be confusing. So, here is a list of only the work related variable names and what they represent. They may represent work or work rate depending on the structure in which they are stored. Note that the work (rate) done on the treadmill belts is always by the legs. But the work (rate) done on the centre of mass is through the legs but may be by the legs or the belt.

net: Total net work (rate) done by the legs\
totpos: Total positive work (rate) done by the legs\
totneg: Total negative work/work rate done by the legs\
lpos: Total positive work /work rate done by the left leg\
lneg: Total negative work/work rate done by the left leg\
rpos: Total positive work/work rate done by the right leg\
rneg: Total negative work/work rate done by the right leg

nettread: Total net work/work rate done on the belts\
treadpos: Total positive work/work rate done on the belts\
treadneg: Total negative work/work rate done on the belts\
ltreadpos: Total positive work/work rate done on the left belt\
ltreadneg: Total negative work/work rate done on the left belt\
rtreadpos: Total positive work/work rate done on the right belt\
rtreadneg: Total negative work/work rate done on the right belt

netcom: Total net work/work rate done on the centre of mass\
compos: Total positive work /work rate done on the centre of mass\
comneg: Total negative work/work rate done on the centre of mass\
lcompos: Total positive work/work rate done on the centre of mass through the left leg\
lcomneg: Total negative work/work rate done on the centre of mass through the left leg\
rcompos: Total positive work/work rate done on the centre of mass through the right leg\
rcomneg: Total negative work/work rate done on the centre of mass through the right leg
