sarter = ['!../BinV201/mwsarta_atms fin=NH3_big0_cris_op_sergio_1014mb.rtp fout=junk.rtp'];
eval(sartaer)

[h3,ha3,p3,pa3] = rtpread('junk.rtp');
p3x = p3;
p3x = rmfield(p3x,'gxover');
p3x = rmfield(p3x,'txover');
p3x = rmfield(p3x,'pnote');
p3x = rmfield(p3x,'udef');
p3x = rmfield(p3x,'iudef');
p3x = rmfield(p3x,'freqcal');
p3x = rmfield(p3x,'robsqual');
p3x = rmfield(p3x,'ifov');
p3x = rmfield(p3x,'xtrack');
p3x = rmfield(p3x,'atrack');
p3x = rmfield(p3x,'sundist');
p3x = rmfield(p3x,'glint');
p3x = rmfield(p3x,'cemis');
p3x = rmfield(p3x,'cemis2');
p3x = rmfield(p3x,'cprtop');
p3x = rmfield(p3x,'cprtop2');
p3x = rmfield(p3x,'cprbot');
p3x = rmfield(p3x,'cprbot2');
p3x = rmfield(p3x,'cfrac');
p3x = rmfield(p3x,'cfrac2');
p3x = rmfield(p3x,'cfrac12');
p3x = rmfield(p3x,'cpsize');
p3x = rmfield(p3x,'cpsize2');
p3x = rmfield(p3x,'cngwat');
p3x = rmfield(p3x,'cngwat2');
p3x = rmfield(p3x,'ctype');
p3x = rmfield(p3x,'ctype2');
p3x = rmfield(p3x,'crho');
p3x = rmfield(p3x,'crho2');
p3x = rmfield(p3x,'cstemp');
p3x = rmfield(p3x,'cstemp2');

save xavier_atms0.mat p3x
