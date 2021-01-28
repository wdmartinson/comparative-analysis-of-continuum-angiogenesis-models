# comparative-analysis-of-continuum-angiogenesis-models

Please note that the functions "Independent_Autonomous_Plots.m" and "Independent_Autonomous_Plots_a_e1.m" require the use of the MATLAB Chebfun library, which can be downloaded at https://www.chebfun.org.

The source code used to generate multiple realizations of the agent-based model outlined in Pillay et al. (2017, Phys. Rev. E) is located in the file "ABM_2D_Model.m". Please note that users must change parameter values within the source code if they wish to generate their own files for investigation.

To generate numerical solutions of the snail-trail and Pillay et al. PDE models described in Martinson et al. (2021, J. Math. Biol. using initial conditions supplied by multiple realizations of the Pillay et al. agent-based model, please run the function "ABM_2D_Model.m" and then "Larger_Domain_Angiogenesis_Models.m". Please note that users must edit the source code directly if they wish to save/load files under a different name than the default. The code makes use of the functions "SnailTrail_1D_PDE.m" and "Pillay_1D_Model.m" to generate solutions. If users wish to simulate the leading order snail-trail/Pillay et al. PDE model described in Martinson et al. (2021, J. Math. Biol.), they must change the boolean variables use_original_STPDE and use_original_PPDE, which are respectively located in the functions "SnailTrail_1D_PDE.m" and "Pillay_1D_Model.m", to TRUE. Please note that the default value of these boolean variables is FALSE.

The function "SnailTrail_1D_PDE_OuterSoln.m" is deprecated and will be removed in a future release.
