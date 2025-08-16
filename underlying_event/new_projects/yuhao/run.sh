#./analyse_decayer ../../filelist/pp_13TeV_HIForced_minbias_0.15mb_rbmSQ0.9_rbmHQ1.4_geoNC3.list hist_outputallFSI_liang_Transmin.root >& log1 &                                           
#./analyse_decayer ../../filelist/pp_13TeV_HIForced_minbias_0.15mb_rbmSQ0.9_rbmHQ1.4_ntmax2_geoNC3.list hist_outputnohFSI_liang_Transmin.root >& log2 &                                    
#./analyse_decayer ../../filelist/pp_13TeV_HIForced_minbias_0mb_ntmax2_geoNC3.list hist_outputnoFSI_liang_Transmin.root >& log3 &                                                          
#./analyse_decayer ../../filelist/pp_13TeV_HIForced_minbias_0.15mb_rbmSQ0.9_rbmHQ1.4_geoNC3.list hist_outputallFSI_liang_NoEqua.root >& log1 &                                             
#./analyse_decayer ../../filelist/pp_13TeV_HIForced_minbias_0.15mb_rbmSQ0.9_rbmHQ1.4_ntmax2_geoNC3.list hist_outputnohFSI_liang_NoEqua.root >& log2 &                                      
#./analyse_decayer ../../filelist/pp_13TeV_HIForced_minbias_0mb_ntmax2_geoNC3.list hist_outputnoFSI_liang_NoEqua.root >& log3 &                                                            
./analyse_decayer /home/adminuser/Working/pythia8_ampt/filelist/pp_13TeV_HIForced_minbias_0mb_rbmCQ1.4_rbmBQ1.2_ntmax2_geoNC3_m0_6.6.list hist_outputnoFSI_new.root >& log1 &
./analyse_decayer /home/adminuser/Working/pythia8_ampt/filelist/pp_13TeV_HIForced_minbias_0.15mb_rbmSQ0.9_rbmCQ1.4_rbmBQ1.2_ntmax2_geoNC3_m0_6.6.list hist_outputnohFSI_new.root >& log2 &
./analyse_decayer /home/adminuser/Working/pythia8_ampt/filelist/pp_13TeV_HIForced_minbias_0mb_ntmax2_geoNC3.list hist_outputnoFSI.root >& logfilenoFSI &
./analyse_decayer /home/adminuser/Working/pythia8_ampt/filelist/pp_13TeV_HIForced_minbias_0.15mb_rbmSQ0.9_rbmHQ1.4_ntmax2_geoNC3.list hist_outputnohFSI.root >& logfilenohFSI &


