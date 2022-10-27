% Figure 2d

load('/Users/julia/Documents/code/ACC-sex-diffs/to_share/Processed_data/InverseGaussian_ACC_DMS_nphr_LAS_zscore0_grouping_none_ver1_laserType1.mat');
fig2_table = table([],{},[],[],[],[],[],[],'VariableNames',{'Animal';'Sex';'scale_ctrl';'shape_ctrl';'shift_ctrl';'scale_laser';'shape_laser';'shift_laser'});
for n = 1:size(fits.fits_f,2)
   fig2_table.Animal(n) = n;
   fig2_table.Sex(n) = {'female'};
   fig2_table.scale_ctrl(n) =fits.fits_f{n}{1}(1);
   fig2_table.shape_ctrl(n)= fits.fits_f{n}{1}(3);
   fig2_table.shift_ctrl(n) =fits.fits_f{n}{1}(2);
   fig2_table.scale_laser(n) =fits.fits_f{n}{2}(1);
   fig2_table.shape_laser(n)= fits.fits_f{n}{2}(3);
   fig2_table.shift_laser(n) =fits.fits_f{n}{2}(2);
end

for n = 1:size(fits.fits_m,2)
   fig2_table.Animal(n+size(fits.fits_f,2)) = n+size(fits.fits_f,2);
   fig2_table.Sex(n+size(fits.fits_f,2)) = {'male'};
   fig2_table.scale_ctrl(n+size(fits.fits_f,2)) =fits.fits_m{n}{1}(1);
   fig2_table.shape_ctrl(n+size(fits.fits_f,2))= fits.fits_m{n}{1}(3);
   fig2_table.shift_ctrl(n+size(fits.fits_f,2)) =fits.fits_m{n}{1}(2);
   fig2_table.scale_laser(n+size(fits.fits_f,2)) =fits.fits_m{n}{2}(1);
   fig2_table.shape_laser(n+size(fits.fits_f,2))= fits.fits_m{n}{2}(3);
   fig2_table.shift_laser(n+size(fits.fits_f,2)) =fits.fits_m{n}{2}(2);
end