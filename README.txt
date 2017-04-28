                                       blOCh 
                                       
                                    Version 5.0

                                   2017 April 26


1.0 CONTENTS

	2.0 INTRODUCTION

	3.0 CITATION
	
	4.0 LICENSE
	
	5.0 CONTRIBUTORS
	
	

	
2.0 INTRODUCTION

	This is the newest public version of blOCh, a MATLAB library for computing and     simulating RF pulses using various optimization methods from optimal control theory and by means of the Bloch equations. It extends to three spatial dimensions and a spectral dimension too. Hence, it can optimize RF pulses for a target defined in 4D.
	blOCh includes four main script files (blOCh__spc, blOCh__khr, blOCh__opt, and blOCh__sim) that can do the basics. And there is a package folder with supplementary scripts for advanced stuff. The four main scripts involve spatial dependencies, temporal dependencies, optimization parameteres, and simulation facilities, respectively. And each of them can be replaced by alternative functions. 
	It is also possible to add more functionality to blOCh, e.g., by including frameworks from other researchers through the packages folder. I would need to consider a request and I would need the corresponding authors acceptance of including the framework in blOCh.
	As many opensource softwares, blOCh does not yet come with a manual and as in an ideal world I should make one. Until then you are welcome to contact me for help and collaborations, and I will setup blOCh for your specific need if I can.
	
	
	
	Best wishes
	
	Mads Vinding






3.0 CITATION

	If you use this software please cite the corresponding papers below.

	For your convenience, look in the folder “citations” for RIS and BIB files.


Application of the limited-memory quasi-Newton algorithm for multi-dimensional, large flip-angle RF pulses at 7T.
Vinding, M.S., et al.,
Magn Reson Mater Phy. 2017, 30, 29-39.
doi:10.1007/s10334-016-0580-1

Local SAR, global SAR, and power-constrained large-flip-angle pulses with optimal control and virtual observation points.
Vinding, M.S., et al.,
Magn Reson Med. 2017, 77, 1, 374-384.
doi: 10.1002/mrm.26086

Real-time 2D spatially selective MRI experiments: Comparative analysis of optimal control design methods.
Maximov, I. et al., 
J Magn Reson. 2015, 254, 110-20. 
doi: 10.1016/j.jmr.2015.03.003

Fast numerical design of spatial-selective rf pulses in MRI using Krotov and quasi-Newton based optimal control methods.
Vinding, M.S., et al., 
J Chem Phys. 2012, 137, 054203.
doi: 10.1063/1.4739755




4.0 LICENSE

	Please see the license file, LICENSE.txt
	
	
5.0 CONTRIBUTORS

	Besides myself blOCh receives invaluable support development from
	
	Dr. Ivan I. Maximov
	MSc. David Goodwin
		



