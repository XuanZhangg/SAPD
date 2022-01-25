/*-------------------------------------------------------------------------------------------------------*/
/*--------------------------------- DrivFace database ---------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------*/


The DrivFace database contains images sequences of subjects while driving in real scenarios. 
It is composed of 606 samples of 640×480 pixels each, acquired over different days from 4 drivers
(2 women and 2 men) with several facial features like glasses and beard.

The ground truth contains the annotation of the face bounding box and the facial key points
(eyes, nose and mouth). 
A set of labels assigning each image into 3 possible gaze direction classes are given. 
The first class is the “looking-right” class and contains the head angles between -45º and -30º. 
The second one is the “frontal” class and contains the head angles between -15º and 15º. 
The last one is the “looking-left” class and contains the head angles between 30º and 45º. 


ATTN: This database is free for academic usage.
      For other purposes, please contact PhD. Katerine Diaz (kdiaz@cvc.uab.es).


Files and scripts
• DrivImages.zip has the driver images. The imag's name has the format:
    * YearMonthDay_subject_Driv_imNum_HeadPose.jpg 

i.e. 20130529_01_Driv_011_f .jpg is a frame of the fisrts driver corresponding to the 11 sequence's image and the head pose is frontal.
subject = [1:4], imNum = [001:...], HeadPose = lr (looking-right), f (frontal) and lf (looking-left).
  

• drivPoints.txt contains the ground truth in table's format, where the columns have the follow information:
    * fileName is the imagen's name into DrivImages.zip
    * subject = [1:4]
    * imgNum  = int
    * label   = [1/2/3] (head pose class that corresponding to [lr/f/lf], respectively) 
    * ang     = [-45, -30/ -15 0 15/ 30 15] (head pose angle)
    * [xF yF wF hF] = face position 
    * [xRE yRE]     = rigth eye position 
    * [xLE yL]      = left eye position 
    * [xN yN]       = Nose position 
    * [xRM yRM]     = rigth corner of mouth
    * [xLM yLM]     = left corner of mouth   


• read_drivPoints.m is a Matlab function to read the drivPoints file. You can also use:
    * Table = readtable('drivPoints.txt');

• DrivFace.zip has all the previous files

Citations
Katerine Diaz-Chito, Aura Hernández-Sabaté, Antonio M. López, A reduced feature set for driver head pose estimation, 
Applied Soft Computing, Volume 45, August 2016, Pages 98-107, ISSN 1568-4946, http://dx.doi.org/10.1016/j.asoc.2016.04.027.  

Contact information 
Questions? Email Katerine Diaz (kdiaz@cvc.uab.es)
