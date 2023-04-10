fname = 'G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\test_extracted\NA@NA\T_Bl_Sq_0to10_0\NA@NA_T_Bl_Sq_0to10_0_202302241751.bin';
minPtsToLoad = 50;
load('G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\test_extracted\NA@NA\T_Bl_Sq_0to10_0\NA@NA_T_Bl_Sq_0to10_0_202302241751 sup data dir\camcalinfo.mat');
cc = camcalinfo;
eset = ExperimentSet.fromFiles(fname, 'minpts', minPtsToLoad, 'camcalinfo', cc, 'parallel', false);
uiopen('G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\test_extracted\NA@NA\T_Bl_Sq_0to10_0\NA@NA_T_Bl_Sq_0to10_0_20230224175_foreground.bmp', 1)
im = imread('G:\AS-Filer\PHY\mmihovil\Shared\Yiming Xu\data\test_extracted\NA@NA\T_Bl_Sq_0to10_0\NA@NA_T_Bl_Sq_0to10_0_20230224175_foreground.bmp');
pcolor(cc.camx,cc.camy, im)  % size doesn't match!!!!!!!!!!!!!!!!!